function result = run_ddvpp_slp_admm_design(params)
%RUN_DDVPP_SLP_ADMM_DESIGN_V2
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% 对应论文框架
% -------------------------------------------------------------------------
% 2.2: 统一节点频率响应模型       -> 通过 evaluate_ddvpp_frequency_response
% 3.2: 特征值敏感度               -> 这里用有限差分实现工程可运行版本
% 3.3: RoCoF / QSS / nadir 约束    -> 外层构造 active constraints
% 4.1: SLP + trust-region          -> 外层 accept / reject + 半径更新
% 4.2: 内层 ADMM                   -> 节点自治更新 + 全局安全投影 + dual update
% 4.3: residual balancing          -> 动态调 rho
%
% -------------------------------------------------------------------------
% 使用方式
% -------------------------------------------------------------------------
% addpath(genpath(pwd));
% result = run_ddvpp_slp_admm_design_v2();
%
% 也可以手动覆盖：
% p = struct();
% p.case_function = 'modified_case118_slp_admm_active';
% p.disturbance_id = 1;
% p.nadir_limit_hz = 0.16;
% p.qss_limit_hz = 0.08;
% result = run_ddvpp_slp_admm_design_v2(p);
%
% -------------------------------------------------------------------------
% 说明
% -------------------------------------------------------------------------
% 1) 该脚本不依赖 quadprog / CVX，仅使用小规模 2x2 盒约束 QP 的闭式/枚举求解。
% 2) 当前仓库 evaluate 函数没有解析式灵敏度接口，因此这里对 sigma(M,D) 使用
%    有限差分近似，这样更稳，也更容易直接跑通。
% 3) 当前 case 中若 nadir / QSS 本来就不违规，则 SLP 仍会工作，但主要优化
%    重心会落在 RoCoF 上。这不是算法坏了，而是 case 设计下活跃约束不同。
%
% Zhongliang 现在最需要的是“能直接替换跑”的版本，所以这里优先保证接口兼容
% 和行为可解释。
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(params)
    params = struct();
end

cfg = local_default_config(params);
base_eval = local_eval_from_MD([], [], cfg);
model = base_eval.model;
gtab = model.ddvpp.gen_dynamic_table;
ctrl_rows = find(gtab.is_controllable);
if isempty(ctrl_rows)
    error('No controllable IBR rows are marked in ddvpp.gen_dynamic_table.is_controllable.');
end

% 读 case 内默认安全参数；允许用户手动覆写。
[cfg.nadir_limit_hz, cfg.qss_limit_hz, cfg.rocof_limit_hz_per_s, ...
 cfg.max_outer, cfg.max_inner, cfg.trust_init, cfg.trust_min, cfg.trust_max, ...
 cfg.rho_init, cfg.primal_tol, cfg.dual_tol] = local_fill_limits_from_case(cfg, base_eval);

% 当前工作点参数
M_curr = base_eval.model.M(:);
D_curr = base_eval.model.D(:);
trust_radius = cfg.trust_init;
outer_rho = cfg.rho_init;

% 初始模态跟踪种子
prev_modes = local_seed_modes(base_eval, cfg.n_tracked_modes);

history = repmat(local_empty_history_row(), 0, 1);

fprintf('===============================================================\n');
fprintf('DDVPP SLP-ADMM design started (v2)\n');
fprintf('Case function            : %s\n', char(string(cfg.case_function)));
fprintf('Disturbance mode         : %s\n', cfg.disturbance_mode);
if strcmpi(cfg.disturbance_mode, 'table_id')
    fprintf('Disturbance table row    : %d\n', cfg.disturbance_id);
end
fprintf('Controllable IBR count   : %d\n', numel(ctrl_rows));
fprintf('Tracked critical modes   : %d\n', cfg.n_tracked_modes);
fprintf('Security limits          : RoCoF<=%.3f, |nadir|<=%.3f, |QSS|<=%.3f\n', ...
    cfg.rocof_limit_hz_per_s, cfg.nadir_limit_hz, cfg.qss_limit_hz);
fprintf('Initial trust radius     : %.4f\n', trust_radius);
fprintf('Initial ADMM rho         : %.4f\n', outer_rho);
fprintf('===============================================================\n\n');

local_print_summary('Baseline', base_eval, cfg);

best_eval = base_eval;
accepted_eval = base_eval;
accepted_M = M_curr;
accepted_D = D_curr;
stop_reason = 'maximum outer iterations reached';

for outer = 1:cfg.max_outer
    fprintf('-------------------- Outer iteration %d --------------------\n', outer);

    curr_eval = local_eval_from_MD(M_curr, D_curr, cfg);
    tracked = local_track_modes(curr_eval, prev_modes, cfg.n_tracked_modes);
    active = local_build_active_constraints(curr_eval, tracked, ctrl_rows, cfg);

    local_print_active_constraints(active, tracked);

    % 若当前已满足全部约束，则可以提前结束。
    if active.total_merit <= cfg.outer_merit_tol
        accepted_eval = curr_eval;
        accepted_M = M_curr;
        accepted_D = D_curr;
        stop_reason = 'all active constraints satisfied';
        hist_row = local_empty_history_row();
        hist_row.outer_iter = outer;
        hist_row.accepted = true;
        hist_row.total_merit_before = active.total_merit;
        hist_row.total_merit_after = active.total_merit;
        hist_row.trust_radius = trust_radius;
        hist_row.rho_final = outer_rho;
        history(end+1) = hist_row; %#ok<AGROW>
        fprintf('All constraints already satisfied. Stop.\n');
        break;
    end

    % ---------------------------------------------------------------------
    % 论文 3.2 节：灵敏度矩阵 S_sigma
    % 这里用有限差分替代解析式，保证和当前仓库 evaluate 接口直接兼容。
    % ---------------------------------------------------------------------
    sens = local_build_sensitivities(curr_eval, tracked, ctrl_rows, cfg);

    % ---------------------------------------------------------------------
    % 论文 4.2 节：内层 ADMM
    % A_i * [dm_i; dd_i] 表示节点 i 对以下“安全贡献量”的线性化贡献：
    %   1) tracked modes 的 sigma 左移量
    %   2) 全局 sum(Delta D)
    %   3) 活跃负荷节点 RoCoF 降幅
    % ---------------------------------------------------------------------
    admm = local_inner_admm(curr_eval, active, sens, ctrl_rows, trust_radius, outer_rho, cfg);

    fprintf('  inner ADMM            : it=%d, primal=%.3e, dual=%.3e, rho_final=%.4f\n', ...
        admm.iters, admm.primal_res, admm.dual_res, admm.rho_final);

    % 形成候选点并重新仿真（论文 4.1 节：SLP 外层评估）
    cand_M = M_curr;
    cand_D = D_curr;
    cand_M(ctrl_rows) = cand_M(ctrl_rows) + admm.dx(:,1);
    cand_D(ctrl_rows) = cand_D(ctrl_rows) + admm.dx(:,2);
    [cand_M, cand_D] = local_clip_MD_to_bounds(cand_M, cand_D, gtab);
    cand_eval = local_eval_from_MD(cand_M, cand_D, cfg);

    cand_active = local_build_active_constraints(cand_eval, tracked, ctrl_rows, cfg);
    [rho_tr, pred_merit] = local_trust_ratio(active, admm, cand_eval, cfg);

    accept = (cand_active.total_merit <= active.total_merit + cfg.accept_tol_abs);
    if rho_tr < cfg.trust_reject_threshold
        accept = false;
    end

    hist_row = local_empty_history_row();
    hist_row.outer_iter = outer;
    hist_row.primal_res = admm.primal_res;
    hist_row.dual_res = admm.dual_res;
    hist_row.rho_final = admm.rho_final;
    hist_row.total_merit_before = active.total_merit;
    hist_row.total_merit_after = cand_active.total_merit;
    hist_row.predicted_merit_after = pred_merit;
    hist_row.rho_tr = rho_tr;
    hist_row.trust_radius = trust_radius;
    hist_row.accepted = accept;
    history(end+1) = hist_row; %#ok<AGROW>

    if accept
        M_curr = cand_M;
        D_curr = cand_D;
        accepted_M = cand_M;
        accepted_D = cand_D;
        accepted_eval = cand_eval;
        prev_modes = local_seed_modes(cand_eval, cfg.n_tracked_modes);

        if rho_tr > 0.85
            trust_radius = min(cfg.trust_expand * trust_radius, cfg.trust_max);
        elseif rho_tr < 0.50
            trust_radius = max(cfg.trust_shrink * trust_radius, cfg.trust_min);
        end

        outer_rho = admm.rho_final;
        fprintf('Accepted step: rho_TR = %.4f, trust radius -> %.4f, rho -> %.4f\n\n', ...
            rho_tr, trust_radius, outer_rho);
        local_print_summary('Accepted iterate', cand_eval, cfg);
        fprintf('\n');
    else
        trust_radius = max(cfg.trust_shrink * trust_radius, cfg.trust_min);
        outer_rho = max(admm.rho_final / 2.0, cfg.rho_min);
        fprintf('Rejected step: rho_TR = %.4f, trust radius -> %.4f, rho -> %.4f\n\n', ...
            rho_tr, trust_radius, outer_rho);
    end

    if cand_active.total_merit <= cfg.outer_merit_tol && accept
        stop_reason = 'all active constraints satisfied';
        break;
    end
end

final_eval = local_eval_from_MD(accepted_M, accepted_D, cfg);
result = struct();
result.cfg = cfg;
result.before = base_eval;
result.after = final_eval;
result.history = history;
result.M_final = accepted_M;
result.D_final = accepted_D;
result.stop_reason = stop_reason;

fprintf('===============================================================\n');
fprintf('DDVPP SLP-ADMM design finished (v2)\n');
fprintf('Stop reason             : %s\n\n', stop_reason);
local_print_summary('Final', final_eval, cfg);

result.update_table = local_build_update_table(base_eval, final_eval, ctrl_rows);
disp('Top controllable nodes with the largest parameter updates:');
disp(result.update_table(1:min(10,height(result.update_table)), :));

result.figures = local_save_figures(base_eval, final_eval, cfg);
fprintf('Saved figures:\n');
fprintf('  %s\n', result.figures.before_all_nodes);
fprintf('  %s\n', result.figures.after_all_nodes);
fprintf('  %s\n', result.figures.compare_worst_buses);
fprintf('===============================================================\n');
end

% =========================================================================
function cfg = local_default_config(params)
cfg = struct();
cfg.case_function = local_get_field(params, 'case_function', 'modified_case118');
cfg.disturbance_mode = local_get_field(params, 'disturbance_mode', 'table_id');
cfg.disturbance_id = local_get_field(params, 'disturbance_id', 1);
cfg.disturbance_load_bus = local_get_field(params, 'disturbance_load_bus', []);
cfg.disturbance_mw = local_get_field(params, 'disturbance_mw', []);
cfg.disturbance_bus = local_get_field(params, 'disturbance_bus', []);
cfg.disturbance_pu = local_get_field(params, 'disturbance_pu', []);

cfg.enforce_local_bounds = true;

cfg.n_tracked_modes = local_get_field(params, 'n_tracked_modes', 3);
cfg.fd_rel_step_M = local_get_field(params, 'fd_rel_step_M', 5e-3);
cfg.fd_rel_step_D = local_get_field(params, 'fd_rel_step_D', 5e-3);

cfg.nadir_limit_hz = local_get_field(params, 'nadir_limit_hz', []);
cfg.qss_limit_hz = local_get_field(params, 'qss_limit_hz', []);
cfg.rocof_limit_hz_per_s = local_get_field(params, 'rocof_limit_hz_per_s', []);

cfg.max_outer = local_get_field(params, 'max_outer', []);
cfg.max_inner = local_get_field(params, 'max_inner', []);
cfg.trust_init = local_get_field(params, 'trust_init', []);
cfg.trust_min = local_get_field(params, 'trust_min', []);
cfg.trust_max = local_get_field(params, 'trust_max', []);
cfg.rho_init = local_get_field(params, 'rho_init', []);
cfg.primal_tol = local_get_field(params, 'primal_tol', []);
cfg.dual_tol = local_get_field(params, 'dual_tol', []);

cfg.cost_M_scale = local_get_field(params, 'cost_M_scale', 1.0);
cfg.cost_D_scale = local_get_field(params, 'cost_D_scale', 5.0);
cfg.cost_lin_scale = local_get_field(params, 'cost_lin_scale', 0.2);

cfg.rocof_active_topk = local_get_field(params, 'rocof_active_topk', 8);
cfg.rocof_gap_tol = local_get_field(params, 'rocof_gap_tol', 1e-5);
cfg.outer_merit_tol = local_get_field(params, 'outer_merit_tol', 1e-5);
cfg.accept_tol_abs = local_get_field(params, 'accept_tol_abs', 1e-8);
cfg.trust_expand = local_get_field(params, 'trust_expand', 1.5);
cfg.trust_shrink = local_get_field(params, 'trust_shrink', 0.5);
cfg.trust_reject_threshold = local_get_field(params, 'trust_reject_threshold', 0.10);

cfg.rho_min = local_get_field(params, 'rho_min', 1e-3);
cfg.rho_max = local_get_field(params, 'rho_max', 1e3);
cfg.rho_balance_mu = local_get_field(params, 'rho_balance_mu', 10.0);
cfg.rho_balance_tau = local_get_field(params, 'rho_balance_tau', 2.0);

cfg.save_dir = local_get_field(params, 'save_dir', 'ddvpp_outputs_v2');
end

% =========================================================================
function eval_out = local_eval_from_MD(M, D, cfg)
p = struct();
p.case_function = cfg.case_function;
p.disturbance_mode = cfg.disturbance_mode;
p.disturbance_id = cfg.disturbance_id;
p.disturbance_load_bus = cfg.disturbance_load_bus;
p.disturbance_mw = cfg.disturbance_mw;
p.disturbance_bus = cfg.disturbance_bus;
p.disturbance_pu = cfg.disturbance_pu;
p.enforce_local_bounds = cfg.enforce_local_bounds;
if ~isempty(M), p.override_M = M; end
if ~isempty(D), p.override_D = D; end
eval_out = evaluate_ddvpp_frequency_response(p);
end

% =========================================================================
function [nadir_lim, qss_lim, rocof_lim, max_outer, max_inner, trust_init, trust_min, trust_max, rho_init, primal_tol, dual_tol] = local_fill_limits_from_case(cfg, eval_out)
sec = eval_out.model.ddvpp.security_limits;

rocof_lim = sec.rocof_limit_hz_per_s;
if ~isempty(cfg.rocof_limit_hz_per_s), rocof_lim = cfg.rocof_limit_hz_per_s; end

nadir_lim = sec.nadir_limit_hz;
if ~isempty(cfg.nadir_limit_hz), nadir_lim = cfg.nadir_limit_hz; end

qss_lim = sec.qss_limit_hz;
if ~isempty(cfg.qss_limit_hz), qss_lim = cfg.qss_limit_hz; end

max_outer = sec.max_iterations_slp;
if ~isempty(cfg.max_outer), max_outer = cfg.max_outer; end

max_inner = sec.max_iterations_admm;
if ~isempty(cfg.max_inner), max_inner = cfg.max_inner; end

trust_init = sec.trust_region_init;
if ~isempty(cfg.trust_init), trust_init = cfg.trust_init; end
trust_min = sec.trust_region_min;
if ~isempty(cfg.trust_min), trust_min = cfg.trust_min; end
trust_max = sec.trust_region_max;
if ~isempty(cfg.trust_max), trust_max = cfg.trust_max; end

rho_init = sec.rho_init;
if ~isempty(cfg.rho_init), rho_init = cfg.rho_init; end

primal_tol = sec.primal_tol;
if ~isempty(cfg.primal_tol), primal_tol = cfg.primal_tol; end
dual_tol = sec.dual_tol;
if ~isempty(cfg.dual_tol), dual_tol = cfg.dual_tol; end
end

% =========================================================================
function modes = local_seed_modes(eval_out, n_track)
modal = eval_out.modal;
model = eval_out.model;
freq_idx = (model.n_gen + 1):(2 * model.n_gen);

cand = modal.dangerous_complex_modes;
if isempty(cand)
    modes = repmat(local_empty_mode_row(), 0, 1);
    return;
end
sel = cand.mode_index(1:min(n_track, height(cand)));
modes = repmat(local_empty_mode_row(), numel(sel), 1);
for ii = 1:numel(sel)
    k = sel(ii);
    modes(ii).mode_index = k;
    modes(ii).lambda = modal.eigvals(k);
    modes(ii).sigma = real(modal.eigvals(k));
    modes(ii).omega_d = imag(modal.eigvals(k));
    modes(ii).vr_freq = modal.Vr(freq_idx, k);
end
end

% =========================================================================
function tracked = local_track_modes(eval_out, prev_modes, n_track)
modal = eval_out.modal;
model = eval_out.model;
freq_idx = (model.n_gen + 1):(2 * model.n_gen);

cand_idx = find(imag(modal.eigvals) > 1e-8);
if isempty(cand_idx)
    tracked = repmat(local_empty_mode_row(), 0, 1);
    return;
end

cand_idx = cand_idx(:).';
cand_v = cell(numel(cand_idx), 1);
for jj = 1:numel(cand_idx)
    cand_v{jj} = modal.Vr(freq_idx, cand_idx(jj));
end

tracked = repmat(local_empty_mode_row(), min(n_track, numel(cand_idx)), 1);
used = false(size(cand_idx));

if isempty(prev_modes)
    prev_modes = repmat(local_empty_mode_row(), 0, 1);
end

n_from_prev = min(numel(prev_modes), numel(tracked));
for ii = 1:n_from_prev
    best_score = -inf;
    best_pos = 0;
    v_prev = prev_modes(ii).vr_freq;
    if isempty(v_prev)
        continue;
    end
    for jj = 1:numel(cand_idx)
        if used(jj), continue; end
        score = local_mac(v_prev, cand_v{jj});
        if score > best_score
            best_score = score;
            best_pos = jj;
        end
    end
    if best_pos > 0
        used(best_pos) = true;
        k = cand_idx(best_pos);
        tracked(ii).mode_index = k;
        tracked(ii).lambda = modal.eigvals(k);
        tracked(ii).sigma = real(modal.eigvals(k));
        tracked(ii).omega_d = imag(modal.eigvals(k));
        tracked(ii).vr_freq = cand_v{best_pos};
    end
end

fill_ptr = n_from_prev + 1;
if fill_ptr <= numel(tracked)
    score = zeros(numel(cand_idx), 1);
    for jj = 1:numel(cand_idx)
        gain = max(abs(modal.mode_step_gain_load(:, cand_idx(jj))));
        score(jj) = gain;
    end
    [~, ord] = sort(score, 'descend');
    for kk = 1:numel(ord)
        jj = ord(kk);
        if used(jj), continue; end
        k = cand_idx(jj);
        tracked(fill_ptr).mode_index = k;
        tracked(fill_ptr).lambda = modal.eigvals(k);
        tracked(fill_ptr).sigma = real(modal.eigvals(k));
        tracked(fill_ptr).omega_d = imag(modal.eigvals(k));
        tracked(fill_ptr).vr_freq = cand_v{jj};
        fill_ptr = fill_ptr + 1;
        if fill_ptr > numel(tracked)
            break;
        end
    end
end
end

% =========================================================================
function active = local_build_active_constraints(eval_out, tracked, ctrl_rows, cfg)
model = eval_out.model;
f0 = model.mpc.userdata.dynamic.base_frequency_hz;
T = eval_out.metrics.per_bus;
resp_load_hz = eval_out.response.resp_load_pu * f0;
qss_hz = resp_load_hz(:, end);
abs_qss = abs(qss_hz);

% --- RoCoF 约束：直接对活跃负荷节点施加线性化约束 ---
rocof = T.rocof0_hz_per_s;
rocof_gap = abs(rocof) - cfg.rocof_limit_hz_per_s;
active_rocof_rows = find(rocof_gap > cfg.rocof_gap_tol);
if ~isempty(active_rocof_rows)
    [~, ord] = sort(rocof_gap(active_rocof_rows), 'descend');
    active_rocof_rows = active_rocof_rows(ord(1:min(cfg.rocof_active_topk, numel(ord))));
end

% --- QSS 约束：使用当前全网最差 QSS gap 构造需求 ---
qss_gap = max(abs_qss) - cfg.qss_limit_hz;
qss_gap = max(qss_gap, 0.0);

% 这里把 QSS gap 折算成对 sum(Delta D) 的需求。
% 近似关系：|w_qss| ≈ |w_qss,0| * D_total0 / (D_total0 + DeltaD_total)
D_total_now = sum(model.D) + sum(model.mu);
worst_qss_now = max(abs_qss);
if qss_gap <= 0 || worst_qss_now <= 1e-12
    req_sum_delta_D = 0.0;
else
    req_sum_delta_D = D_total_now * (worst_qss_now / cfg.qss_limit_hz - 1.0);
    req_sum_delta_D = max(req_sum_delta_D, 0.0);
end

% --- nadir / sigma 约束：只有当前真的越限时才激活 ---
tracked_sigma_req = zeros(numel(tracked), 1);
tracked_bus = zeros(numel(tracked), 1);
tracked_margin = zeros(numel(tracked), 1);
worst_nadir_hz = eval_out.metrics.summary.worst_nadir_hz;
worst_nadir_bus = eval_out.metrics.summary.worst_nadir_bus;
coi_nadir_hz = eval_out.metrics.summary.coi_nadir_hz;
worst_nadir_gap = abs(worst_nadir_hz) - cfg.nadir_limit_hz;
if worst_nadir_gap > 0
    bus_idx = worst_nadir_bus;
    t_nadir = max(T.nadir_time_s(T.bus == bus_idx), 1e-3);
    margin = max(cfg.nadir_limit_hz - abs(coi_nadir_hz), 1e-6);
    for kk = 1:numel(tracked)
        k = tracked(kk).mode_index;
        if k <= 0, continue; end
        Cik = abs(eval_out.modal.residue_load(bus_idx, k)) * f0;
        if Cik <= 1e-12
            tracked_sigma_req(kk) = 0.0;
        else
            sigma_safe = log(margin / max(2.0 * Cik, 1e-12)) / t_nadir;
            tracked_sigma_req(kk) = max(sigma_safe - real(eval_out.modal.eigvals(k)), 0.0);
        end
        tracked_bus(kk) = bus_idx;
        tracked_margin(kk) = margin;
    end
end

% merit 用于外层 accept/reject
merit = sum(max(rocof_gap, 0)) + max(worst_nadir_gap, 0) + max(qss_gap, 0);

active = struct();
active.rocof = rocof;
active.rocof_gap = max(rocof_gap, 0.0);
active.active_rocof_rows = active_rocof_rows(:);
active.qss_hz = qss_hz;
active.req_sum_delta_D = req_sum_delta_D;
active.tracked_sigma_req = tracked_sigma_req(:);
active.tracked_bus = tracked_bus(:);
active.tracked_margin = tracked_margin(:);
active.worst_nadir_bus = worst_nadir_bus;
active.worst_nadir_hz = worst_nadir_hz;
active.worst_qss_bus = find(abs_qss == max(abs_qss), 1, 'first');
active.worst_qss_hz = qss_hz(active.worst_qss_bus);
active.total_merit = merit;
active.ctrl_rows = ctrl_rows(:);
end

% =========================================================================
function local_print_active_constraints(active, tracked)
fprintf('Tracked critical modes and required left shift:\n');
for kk = 1:numel(tracked)
    fprintf('  mode %d : idx=%d, lambda=%+.6f %+.6fj, req Delta sigma >= %.6f\n', ...
        kk, tracked(kk).mode_index, real(tracked(kk).lambda), imag(tracked(kk).lambda), ...
        active.tracked_sigma_req(kk));
end
fprintf('  worst-nadir bus        : %d\n', active.worst_nadir_bus);
fprintf('  current nadir          : %.6f Hz\n', active.worst_nadir_hz);
fprintf('  required sum(Delta D)  : %.6f\n', active.req_sum_delta_D);
if isempty(active.active_rocof_rows)
    fprintf('  active RoCoF buses      : none\n');
else
    fprintf('  active RoCoF buses      : ');
    fprintf('%d ', active.active_rocof_rows);
    fprintf('\n');
end
end

% =========================================================================
function sens = local_build_sensitivities(base_eval, tracked, ctrl_rows, cfg)
model = base_eval.model;
M = model.M(:);
D = model.D(:);

n_track = numel(tracked);
n_ctrl = numel(ctrl_rows);
Ssig_M = zeros(n_track, n_ctrl);
Ssig_D = zeros(n_track, n_ctrl);

for ii = 1:n_ctrl
    gidx = ctrl_rows(ii);

    % --- d sigma / d M_i ---
    dM = max(cfg.fd_rel_step_M * max(abs(M(gidx)), 1.0), 1e-4);
    Mp = M; Mp(gidx) = Mp(gidx) + dM;
    eval_M = local_eval_from_MD(Mp, D, cfg);
    tracked_M = local_track_modes(eval_M, tracked, n_track);
    for kk = 1:n_track
        if tracked(kk).mode_index <= 0 || tracked_M(kk).mode_index <= 0
            continue;
        end
        Ssig_M(kk, ii) = (real(tracked_M(kk).lambda) - real(tracked(kk).lambda)) / dM;
    end

    % --- d sigma / d D_i ---
    dD = max(cfg.fd_rel_step_D * max(abs(D(gidx)), 1.0), 1e-5);
    Dp = D; Dp(gidx) = Dp(gidx) + dD;
    eval_D = local_eval_from_MD(M, Dp, cfg);
    tracked_D = local_track_modes(eval_D, tracked, n_track);
    for kk = 1:n_track
        if tracked(kk).mode_index <= 0 || tracked_D(kk).mode_index <= 0
            continue;
        end
        Ssig_D(kk, ii) = (real(tracked_D(kk).lambda) - real(tracked(kk).lambda)) / dD;
    end
end

% --- 负荷节点 RoCoF 对 M 的灵敏度：解析线性化 ---
f0 = model.mpc.userdata.dynamic.base_frequency_hz;
F = model.F;
M_now = model.M(:);
u = model.M(:) .* base_eval.response.bstep((model.n_gen + 1):(2 * model.n_gen));
rocof_sens_load = zeros(model.n_load, n_ctrl);
for ii = 1:n_ctrl
    gidx = ctrl_rows(ii);
    rocof_sens_load(:, ii) = -f0 * F(:, gidx) * (u(gidx) / (M_now(gidx)^2));
end

sens = struct();
sens.Ssig_M = Ssig_M;
sens.Ssig_D = Ssig_D;
sens.rocof_sens_load = rocof_sens_load;
end

% =========================================================================
function admm = local_inner_admm(curr_eval, active, sens, ctrl_rows, trust_radius, rho_init, cfg)
model = curr_eval.model;
gtab = model.ddvpp.gen_dynamic_table;
f0 = model.mpc.userdata.dynamic.base_frequency_hz;
M = model.M(:);
D = model.D(:);

n_ctrl = numel(ctrl_rows);
n_track = size(sens.Ssig_M, 1);
act_rocof_rows = active.active_rocof_rows(:);
n_rocof = numel(act_rocof_rows);

% --- 组装“每节点”的局部贡献矩阵 A_i ---
% 行含义：
%   1..n_track         : sigma 左移贡献
%   n_track + 1        : Delta D 总和贡献
%   n_track + 2 .. end : 活跃 RoCoF 约束的降幅贡献
m = n_track + 1 + n_rocof;
A = cell(n_ctrl, 1);

rocof_sign = sign(active.rocof(act_rocof_rows));
rocof_sign(rocof_sign == 0) = 1;

for ii = 1:n_ctrl
    Ai = zeros(m, 2);
    if n_track > 0
        Ai(1:n_track, 1) = sens.Ssig_M(:, ii);
        Ai(1:n_track, 2) = sens.Ssig_D(:, ii);
    end
    Ai(n_track + 1, 2) = 1.0;
    for rr = 1:n_rocof
        bus_row = act_rocof_rows(rr);
        % 目标写成“正贡献 = 减小 |RoCoF|”
        Ai(n_track + 1 + rr, 1) = -rocof_sign(rr) * sens.rocof_sens_load(bus_row, ii);
    end
    A{ii} = Ai;
end

% --- 全局每节点共识下界 z_lb_per_node ---
% 采用 consensus ADMM：每个节点的贡献 A_i x_i 与同一个 z 对齐。
% 因而总需求先均分到每个节点的共识目标上。
z_lb_total = [active.tracked_sigma_req(:); active.req_sum_delta_D; active.rocof_gap(act_rocof_rows)];
z_lb = z_lb_total / max(n_ctrl, 1);

% --- 每个节点的本地盒约束 ---
ctrl_tbl_rows = ctrl_rows(:);
M_ctrl = M(ctrl_tbl_rows);
D_ctrl = D(ctrl_tbl_rows);
q_quad = gtab.local_cost_quad(ctrl_tbl_rows);
q_lin = gtab.local_cost_lin(ctrl_tbl_rows);

% generator-side 本地下界（论文 3.3 式(10) 的工程版）
% 这里用 u_i / M_i 对 generator node RoCoF 的硬下界，避免某些节点过小惯量。
u = model.M(:) .* curr_eval.response.bstep((model.n_gen + 1):(2 * model.n_gen));
local_m_req = f0 * abs(u(ctrl_tbl_rows)) / max(cfg.rocof_limit_hz_per_s, 1e-9);

lb = zeros(n_ctrl, 2);
ub = zeros(n_ctrl, 2);
for ii = 1:n_ctrl
    gidx = ctrl_tbl_rows(ii);
    lb(ii,1) = max([gtab.m_min(gidx) - M(gidx), local_m_req(ii) - M(gidx), -trust_radius]);
    ub(ii,1) = min([gtab.m_max(gidx) - M(gidx), +trust_radius]);
    lb(ii,2) = max([gtab.d_min(gidx) - D(gidx), -trust_radius]);
    ub(ii,2) = min([gtab.d_max(gidx) - D(gidx), +trust_radius]);
    if lb(ii,1) > ub(ii,1)
        lb(ii,1) = ub(ii,1);
    end
    if lb(ii,2) > ub(ii,2)
        lb(ii,2) = ub(ii,2);
    end
end

% --- ADMM 变量初始化 ---
z = z_lb;
y = zeros(m, n_ctrl);   % scaled dual per node
x = zeros(n_ctrl, 2);
rho = rho_init;
primal_res = inf;
dual_res = inf;

for it = 1:cfg.max_inner
    % Step 1: 节点自治更新 x-update
    for ii = 1:n_ctrl
        Ai = A{ii};
        H = diag([cfg.cost_M_scale * max(q_quad(ii), 1e-8), cfg.cost_D_scale * max(q_quad(ii), 1e-8)]) + rho * (Ai' * Ai);
        g = [cfg.cost_lin_scale * q_lin(ii); cfg.cost_lin_scale * q_lin(ii)] + rho * Ai' * (y(:,ii) - z);
        x(ii,:) = local_box_qp_2d(H, g, lb(ii,:).', ub(ii,:).').';
    end

    % Step 2: 全网安全投影 z-update
    z_prev = z;
    V = zeros(m, n_ctrl);
    for ii = 1:n_ctrl
        V(:,ii) = A{ii} * x(ii,:).' + y(:,ii);
    end
    z = mean(V, 2);
    z = max(z, z_lb);  % 对应论文 4.2 的投影到 Z_safe；本例是 lower-bound half-space

    % Step 3: dual update
    r_accum = 0.0;
    for ii = 1:n_ctrl
        ri = A{ii} * x(ii,:).' - z;
        y(:,ii) = y(:,ii) + ri;
        r_accum = r_accum + sum(ri.^2);
    end
    primal_res = sqrt(r_accum);
    dual_res = rho * sqrt(n_ctrl) * norm(z - z_prev, 2);

    % 论文 4.3：Residual Balancing
    if primal_res > cfg.rho_balance_mu * dual_res && rho < cfg.rho_max
        rho_new = min(cfg.rho_balance_tau * rho, cfg.rho_max);
        scale = rho / rho_new;
        y = y * scale;  % scaled dual 需要同步缩放
        rho = rho_new;
    elseif dual_res > cfg.rho_balance_mu * primal_res && rho > cfg.rho_min
        rho_new = max(rho / cfg.rho_balance_tau, cfg.rho_min);
        scale = rho / rho_new;
        y = y * scale;
        rho = rho_new;
    end

    if primal_res <= cfg.primal_tol && dual_res <= cfg.dual_tol
        break;
    end
end

admm = struct();
admm.dx = x;
admm.z = z;
admm.z_lb = z_lb;
admm.iters = it;
admm.rho_final = rho;
admm.primal_res = primal_res;
admm.dual_res = dual_res;

% 方便外层做 predicted improvement
admm.total_contrib = zeros(m, 1);
for ii = 1:n_ctrl
    admm.total_contrib = admm.total_contrib + A{ii} * x(ii,:).';
end
end

% =========================================================================
function [rho_tr, pred_merit] = local_trust_ratio(active_before, admm, cand_eval, cfg)
% 这里只构造一个工程上稳定的 trust-ratio：
% 1) 对 RoCoF gap 用线性预测与实际仿真对比
% 2) 没有 RoCoF 活跃时退化为整体 merit 改善率

pred_rocof_after = active_before.rocof_gap;
act_rows = active_before.active_rocof_rows(:);
if ~isempty(act_rows)
    contrib_rocof_total = admm.total_contrib((numel(active_before.tracked_sigma_req) + 2):end);
    pred_rocof_after(act_rows) = max(active_before.rocof_gap(act_rows) - contrib_rocof_total, 0.0);
end

f0 = cand_eval.model.mpc.userdata.dynamic.base_frequency_hz;
T_after = cand_eval.metrics.per_bus;
resp_load_hz_after = cand_eval.response.resp_load_pu * f0;
qss_after = max(abs(resp_load_hz_after(:, end))) - cfg.qss_limit_hz;
qss_after = max(qss_after, 0.0);
nadir_after = max(abs(cand_eval.metrics.summary.worst_nadir_hz) - cfg.nadir_limit_hz, 0.0);
rocof_after = max(abs(T_after.rocof0_hz_per_s) - cfg.rocof_limit_hz_per_s, 0.0);
actual_merit_after = sum(rocof_after) + qss_after + nadir_after;

pred_merit = sum(pred_rocof_after) + max(abs(active_before.worst_qss_hz) - cfg.qss_limit_hz, 0.0) ...
    + max(abs(active_before.worst_nadir_hz) - cfg.nadir_limit_hz, 0.0);

num = active_before.total_merit - actual_merit_after;
den = active_before.total_merit - pred_merit;
if den <= 1e-12
    rho_tr = 1.0;
else
    rho_tr = num / den;
end
rho_tr = max(min(rho_tr, 2.0), -2.0);
end

% =========================================================================
function local_print_summary(tag, eval_out, cfg)
model = eval_out.model;
f0 = model.mpc.userdata.dynamic.base_frequency_hz;
resp_load_hz = eval_out.response.resp_load_pu * f0;
qss_hz = resp_load_hz(:, end);
[worst_qss_abs, qss_row] = max(abs(qss_hz));
summary = eval_out.metrics.summary;

fprintf('[%s]\n', tag);
fprintf('  Worst nadir bus        : %d, %.4f Hz\n', summary.worst_nadir_bus, summary.worst_nadir_hz);
fprintf('  Worst |RoCoF| bus      : %d, %.4f Hz/s\n', summary.worst_rocof_bus, abs(summary.worst_rocof_hz_per_s));
fprintf('  Worst |QSS| bus        : %d, %.4f Hz\n', qss_row, worst_qss_abs);
fprintf('  COI nadir              : %.4f Hz @ %.4f s\n', summary.coi_nadir_hz, summary.coi_nadir_time_s);
fprintf('  Security check         : RoCoF<=%.3f, |nadir|<=%.3f, |QSS|<=%.3f\n', ...
    cfg.rocof_limit_hz_per_s, cfg.nadir_limit_hz, cfg.qss_limit_hz);
end

% =========================================================================
function T = local_build_update_table(before_eval, after_eval, ctrl_rows)
g0 = before_eval.model.ddvpp.gen_dynamic_table;
Tb = table();
Tb.gen_row = g0.gen_index(ctrl_rows);
Tb.host_bus = g0.host_bus(ctrl_rows);
Tb.M_before = before_eval.model.M(ctrl_rows);
Tb.M_after = after_eval.model.M(ctrl_rows);
Tb.Delta_M = Tb.M_after - Tb.M_before;
Tb.D_before = before_eval.model.D(ctrl_rows);
Tb.D_after = after_eval.model.D(ctrl_rows);
Tb.Delta_D = Tb.D_after - Tb.D_before;
[~, ord] = sort(abs(Tb.Delta_M) + abs(Tb.Delta_D), 'descend');
T = Tb(ord, :);
end

% =========================================================================
function figs = local_save_figures(before_eval, after_eval, cfg)
save_dir = cfg.save_dir;
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

f0 = before_eval.model.mpc.userdata.dynamic.base_frequency_hz;
t_before = before_eval.response.t;
resp_before = before_eval.response.resp_load_pu * f0;
t_after = after_eval.response.t;
resp_after = after_eval.response.resp_load_pu * f0;

fig1 = figure('Color', 'w', 'Name', 'All buses before');
plot(t_before, resp_before, 'LineWidth', 0.6);
grid on; xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('All load-bus frequency responses before optimization');
path1 = fullfile(save_dir, 'ddvpp_freq_before_all_nodes.png');
exportgraphics(fig1, path1, 'Resolution', 220);
% close(fig1);

fig2 = figure('Color', 'w', 'Name', 'All buses after');
plot(t_after, resp_after, 'LineWidth', 0.6);
grid on; xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('All load-bus frequency responses after optimization');
path2 = fullfile(save_dir, 'ddvpp_freq_after_all_nodes.png');
exportgraphics(fig2, path2, 'Resolution', 220);
% close(fig2);

[~, wb_before] = min(before_eval.metrics.per_bus.nadir_hz);
[~, wb_after] = min(after_eval.metrics.per_bus.nadir_hz);
idx1 = before_eval.metrics.per_bus.bus(wb_before);
idx2 = after_eval.metrics.per_bus.bus(wb_after);

fig3 = figure('Color', 'w', 'Name', 'Worst buses compare');
plot(t_before, resp_before(idx1, :), 'LineWidth', 1.6); hold on;
plot(t_after, resp_after(idx1, :), '--', 'LineWidth', 1.6);
if idx2 ~= idx1
    plot(t_before, resp_before(idx2, :), 'LineWidth', 1.2);
    plot(t_after, resp_after(idx2, :), '--', 'LineWidth', 1.2);
    legend({sprintf('Before worst bus %d', idx1), sprintf('After on bus %d', idx1), ...
            sprintf('Before bus %d', idx2), sprintf('After worst bus %d', idx2)}, ...
            'Location', 'best');
else
    legend({sprintf('Before bus %d', idx1), sprintf('After bus %d', idx1)}, 'Location', 'best');
end
grid on; xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Worst-bus frequency response comparison');
path3 = fullfile(save_dir, 'ddvpp_freq_compare_worst_buses.png');
exportgraphics(fig3, path3, 'Resolution', 220);
% close(fig3);

figs = struct();
figs.before_all_nodes = path1;
figs.after_all_nodes = path2;
figs.compare_worst_buses = path3;
end

% =========================================================================
function [M_clip, D_clip] = local_clip_MD_to_bounds(M, D, gtab)
M_clip = min(max(M(:), gtab.m_min(:)), gtab.m_max(:));
D_clip = min(max(D(:), gtab.d_min(:)), gtab.d_max(:));
end

% =========================================================================
function x = local_box_qp_2d(H, g, lb, ub)
% 精确求解 2 维 box-constrained convex QP:
%   min 0.5 x'Hx + g'x, s.t. lb <= x <= ub
cands = zeros(0, 2);

% unconstrained
xu = -H \ g;
if all(xu >= lb - 1e-12) && all(xu <= ub + 1e-12)
    cands(end+1,:) = xu.'; %#ok<AGROW>
end

% x1 fixed, solve x2
for x1 = [lb(1), ub(1)]
    if H(2,2) > 0
        x2 = -(H(2,1) * x1 + g(2)) / H(2,2);
        x2 = min(max(x2, lb(2)), ub(2));
        cands(end+1,:) = [x1, x2]; %#ok<AGROW>
    end
end

% x2 fixed, solve x1
for x2 = [lb(2), ub(2)]
    if H(1,1) > 0
        x1 = -(H(1,2) * x2 + g(1)) / H(1,1);
        x1 = min(max(x1, lb(1)), ub(1));
        cands(end+1,:) = [x1, x2]; %#ok<AGROW>
    end
end

% corners
cands(end+1,:) = [lb(1), lb(2)]; %#ok<AGROW>
cands(end+1,:) = [lb(1), ub(2)]; %#ok<AGROW>
cands(end+1,:) = [ub(1), lb(2)]; %#ok<AGROW>
cands(end+1,:) = [ub(1), ub(2)]; %#ok<AGROW>

best_val = inf;
best_x = [0; 0];
for ii = 1:size(cands, 1)
    xk = cands(ii,:).';
    val = 0.5 * xk.' * H * xk + g.' * xk;
    if val < best_val
        best_val = val;
        best_x = xk;
    end
end
x = best_x;
end

% =========================================================================
function score = local_mac(v1, v2)
num = abs(v1' * v2)^2;
den = real((v1' * v1) * (v2' * v2));
if den <= 1e-20
    score = 0.0;
else
    score = real(num / den);
end
end

% =========================================================================
function row = local_empty_mode_row()
row = struct('mode_index', 0, 'lambda', 0, 'sigma', 0, 'omega_d', 0, 'vr_freq', []);
end

% =========================================================================
function row = local_empty_history_row()
row = struct( ...
    'outer_iter', 0, ...
    'primal_res', NaN, ...
    'dual_res', NaN, ...
    'rho_final', NaN, ...
    'rho_tr', NaN, ...
    'trust_radius', NaN, ...
    'total_merit_before', NaN, ...
    'predicted_merit_after', NaN, ...
    'total_merit_after', NaN, ...
    'accepted', false);
end

% =========================================================================
function v = local_get_field(s, name, default_val)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default_val;
end
end
