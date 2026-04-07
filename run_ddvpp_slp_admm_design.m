
function result = run_ddvpp_slp_admm_design(user_opt)
%RUN_DDVPP_SLP_ADMM_DESIGN
% DDVPP paper driver that extends the existing frequency-response evaluator
% with the remaining optimization framework:
%   - Section 3.2: eigenvalue sensitivity of virtual inertia / damping
%   - Section 3.3: RoCoF / QSS / nadir security constraints
%   - Section 4.1: outer SLP loop with MAC mode tracking and trust region
%   - Section 4.2: inner distributed ADMM-like coordination
%
% This file is intentionally written as ONE self-contained .m file so that
% you can drop it into the current GitHub repository and run it directly.
%
% Required existing files in the same MATLAB path:
%   evaluate_ddvpp_frequency_response.m
%   modified_case118.m
%
% Typical use
% -------------------------------------------------------------------------
%   addpath(genpath(pwd));
%   result = run_ddvpp_slp_admm_design();
%
%   % Or override a few options:
%   opt = struct();
%   opt.disturbance_id = 1;          % use ddvpp.disturbance_set row 1
%   opt.n_critical_modes = 3;        % track top-3 dangerous oscillatory modes
%   opt.save_dir = 'ddvpp_outputs';
%   result = run_ddvpp_slp_admm_design(opt);
%
% Output
% -------------------------------------------------------------------------
%   result.base_eval      : baseline response before optimization
%   result.final_eval     : final response after optimization
%   result.solution       : optimized M / D and increment table
%   result.history        : outer-loop iteration history
%   result.figures        : figure handles
%   result.options        : resolved options
%
% Notes
% -------------------------------------------------------------------------
% 1) Only M and D are optimized here, which matches the current paper focus
%    on IBR virtual inertia / virtual damping contribution.
% 2) To keep the script easy to run, the inner QP is solved with a custom
%    2-D box-constrained solver, so no Optimization Toolbox is required.
% 3) The script prints key results directly in the MATLAB command window and
%    saves the before/after frequency-response figures to opt.save_dir.

    if nargin < 1 || isempty(user_opt)
        user_opt = struct();
    end

    opt = local_default_options();
    opt = local_merge_struct(opt, user_opt);

    % ---------------------------------------------------------------------
    % Pre-checks against the existing repository interface
    % ---------------------------------------------------------------------
    if exist('evaluate_ddvpp_frequency_response', 'file') ~= 2
        error(['evaluate_ddvpp_frequency_response.m was not found in the MATLAB path. ' ...
               'Please put this new file in the same repository folder and addpath(genpath(pwd)).']);
    end
    if exist(char(opt.case_function), 'file') ~= 2
        error(['The case function "', char(opt.case_function), '" was not found. ' ...
               'Please make sure modified_case118.m is on the MATLAB path.']);
    end

    base_eval = local_run_eval(opt, [], []);
    local_assert_required_fields(base_eval);

    limits = base_eval.model.ddvpp.security_limits;
    gtab = base_eval.model.ddvpp.gen_dynamic_table;
    ctrl_idx = find(gtab.is_controllable);
    if isempty(ctrl_idx)
        error('No controllable IBR nodes were found in ddvpp.gen_dynamic_table.');
    end

    M_cur = base_eval.model.M;
    D_cur = base_eval.model.D;

    trust_radius = opt.trust_region_init;
    rho = opt.rho_init;
    prev_modes = [];

    % 用 0x1 的同字段结构体数组初始化 history，避免后续
    % history(end+1)=hist_row 时触发“不同结构体之间进行下标赋值”。
    history = repmat(local_empty_history_row(), 0, 1);

    fprintf('\n===============================================================\n');
    fprintf('DDVPP SLP-ADMM design started\n');
    fprintf('Case function            : %s\n', char(opt.case_function));
    fprintf('Disturbance mode         : %s\n', char(opt.disturbance_mode));
    if strcmpi(opt.disturbance_mode, 'table_id')
        fprintf('Disturbance table row    : %d\n', opt.disturbance_id);
    end
    fprintf('Controllable IBR count   : %d\n', numel(ctrl_idx));
    fprintf('Tracked critical modes   : %d\n', opt.n_critical_modes);
    fprintf('Initial trust radius     : %.4f\n', trust_radius);
    fprintf('Initial ADMM rho         : %.4f\n', rho);
    fprintf('===============================================================\n\n');

    local_print_screening('Baseline', base_eval, limits);

    % ---------------------------------------------------------------------
    % Outer loop: SLP (paper Section 4.1)
    % ---------------------------------------------------------------------
    accepted_eval = base_eval;
    accepted_modes = [];
    stop_reason = 'maximum outer iterations reached';

    for v = 1:opt.max_outer_iter
        cur_eval = local_run_eval(opt, M_cur, D_cur);

        if isempty(prev_modes)
            cur_modes = local_select_initial_modes(cur_eval, opt.n_critical_modes);
        else
            cur_modes = local_track_modes_by_mac(cur_eval, prev_modes, opt.n_critical_modes);
        end

        sens = local_build_mode_sensitivities(cur_eval, cur_modes, ctrl_idx);
        targets = local_build_security_targets(cur_eval, cur_modes, ctrl_idx, opt);

        fprintf('-------------------- Outer iteration %d --------------------\n', v);
        local_print_modes(cur_modes, targets);

        if targets.already_secure
            fprintf('Current operating point already satisfies RoCoF / nadir / QSS limits.\n');
            accepted_eval = cur_eval;
            accepted_modes = cur_modes;
            stop_reason = 'security limits already satisfied';
            break;
        end

        inner = local_run_inner_admm(cur_eval, ctrl_idx, sens, targets, trust_radius, rho, opt);
        fprintf('  inner ADMM            : it=%d, primal=%.3e, dual=%.3e, rho_final=%.4f\n', ...
            inner.iterations, inner.primal_residual, inner.dual_residual, inner.rho_final);

        dM = zeros(size(M_cur));
        dD = zeros(size(D_cur));
        dM(ctrl_idx) = inner.dx(:,1);
        dD(ctrl_idx) = inner.dx(:,2);

        M_trial = M_cur + dM;
        D_trial = D_cur + dD;

        trial_eval = local_run_eval(opt, M_trial, D_trial);
        trial_modes = local_track_modes_by_mac(trial_eval, cur_modes, opt.n_critical_modes);

        pred_sigma_shift = inner.sum_contrib(1:numel(cur_modes));
        act_sigma_shift = local_actual_sigma_shift(cur_modes, trial_modes);
        [step_accepted, trust_radius_next, rho_tr] = local_trust_region_update( ...
            pred_sigma_shift, act_sigma_shift, trust_radius, opt);

        hist_row = local_empty_history_row();
        hist_row.outer_iter = v;
        hist_row.trust_radius_before = trust_radius;
        hist_row.trust_ratio = rho_tr;
        hist_row.inner_iterations = inner.iterations;
        hist_row.rho_after = inner.rho_final;
        hist_row.pred_shift_norm = norm(pred_sigma_shift, 2);
        hist_row.act_shift_norm = norm(act_sigma_shift, 2);
        hist_row.sigma_req_min = min(targets.req_sigma_shift);
        hist_row.worst_nadir_before = min(cur_eval.metrics.per_bus.nadir_hz);
        load_hz_cur = local_get_resp_load_hz(cur_eval);
        load_hz_trial = local_get_resp_load_hz(trial_eval);
        hist_row.worst_nadir_after = min(load_hz_trial(:));
        hist_row.qss_abs_before = max(abs(load_hz_cur(:, end)));
        hist_row.qss_abs_after = max(abs(load_hz_trial(:, end)));
        hist_row.accepted = step_accepted;

        if step_accepted
            M_cur = M_trial;
            D_cur = D_trial;
            accepted_eval = trial_eval;
            accepted_modes = trial_modes;
            prev_modes = trial_modes;
            trust_radius = trust_radius_next;
            rho = inner.rho_final;

            fprintf('Accepted step: rho_TR = %.4f, trust radius -> %.4f, rho -> %.4f\n', ...
                rho_tr, trust_radius, rho);
            local_print_screening('Accepted iterate', trial_eval, limits);

            if local_is_secure(trial_eval, limits, opt)
                stop_reason = 'all security limits satisfied after accepted step';
                history(end+1) = hist_row; %#ok<AGROW>
                break;
            end
        else
            prev_modes = cur_modes;
            trust_radius = trust_radius_next;

            fprintf('Rejected step: rho_TR = %.4f, trust radius shrunk -> %.4f\n', ...
                rho_tr, trust_radius);

            if trust_radius <= opt.trust_region_min + 1e-12
                stop_reason = 'trust region reached minimum after a rejected step';
                history(end+1) = hist_row; %#ok<AGROW>
                break;
            end
        end

        history(end+1) = hist_row; %#ok<AGROW>

        if v == opt.max_outer_iter
            stop_reason = 'maximum outer iterations reached';
        end
    end

    final_eval = local_run_eval(opt, M_cur, D_cur);
    if isempty(accepted_modes)
        accepted_modes = local_select_initial_modes(final_eval, opt.n_critical_modes); %#ok<NASGU>
    end

    figures = local_plot_before_after(base_eval, final_eval, opt);
    solution = local_build_solution_table(base_eval, final_eval, ctrl_idx);
    local_print_final_report(final_eval, limits, solution, stop_reason, opt);

    result = struct();
    result.base_eval = base_eval;
    result.final_eval = final_eval;
    result.solution = solution;
    result.history = history;
    result.figures = figures;
    result.options = opt;
    result.stop_reason = stop_reason;
end

%% ========================================================================
function opt = local_default_options()
% Centralized option block for the new DDVPP driver.
% The defaults are chosen to match the existing repository interface.

    opt = struct();

    % Existing evaluator interface
    opt.case_function = 'modified_case118';
    opt.disturbance_mode = 'table_id';
    opt.disturbance_id = 1;
    opt.disturbance_load_bus = [];
    opt.disturbance_mw = [];
    opt.disturbance_bus = [];
    opt.disturbance_pu = [];
    opt.duration_s = [];
    opt.n_points = [];
    opt.enforce_local_bounds = true;

    % Paper Section 4.1: outer SLP settings
    opt.max_outer_iter = 12;
    opt.n_critical_modes = 3;
    opt.trust_region_init = 0.10;
    opt.trust_region_min = 0.01;
    opt.trust_region_max = 0.50;
    opt.trust_expand_factor = 1.50;
    opt.trust_shrink_factor = 0.50;
    opt.trust_accept_ratio = 0.10;
    opt.trust_expand_ratio = 0.75;

    % Paper Section 4.2: inner ADMM settings
    opt.max_inner_iter = 150;
    opt.rho_init = 1.00;
    opt.rho_min = 1e-3;
    opt.rho_max = 1e3;
    opt.rho_mu = 10.0;
    opt.rho_tau = 2.0;
    opt.primal_tol = 1e-4;
    opt.dual_tol = 1e-4;

    % Local quadratic cost weights
    opt.weight_M = 1.0;
    opt.weight_D = 1.0;
    opt.min_quad_weight = 1e-6;

    % Numerical guards
    opt.min_margin_hz = 1e-4;
    opt.min_time_s = 0.10;

    % Plot / output
    opt.save_dir = 'ddvpp_outputs';
    opt.topk_show = 10;
end

%% ========================================================================
function eval_out = local_run_eval(opt, M_override, D_override)
% Thin wrapper around the existing evaluator in the repository.

    p = struct();
    p.case_function = opt.case_function;
    p.disturbance_mode = opt.disturbance_mode;
    p.disturbance_id = opt.disturbance_id;
    p.disturbance_load_bus = opt.disturbance_load_bus;
    p.disturbance_mw = opt.disturbance_mw;
    p.disturbance_bus = opt.disturbance_bus;
    p.disturbance_pu = opt.disturbance_pu;
    p.enforce_local_bounds = opt.enforce_local_bounds;

    if ~isempty(opt.duration_s), p.duration_s = opt.duration_s; end
    if ~isempty(opt.n_points),   p.n_points = opt.n_points;   end
    if ~isempty(M_override),     p.override_M = M_override(:); end
    if ~isempty(D_override),     p.override_D = D_override(:); end

    eval_out = evaluate_ddvpp_frequency_response(p);
end

%% ========================================================================
function local_assert_required_fields(eval_out)
% Make sure the current repository output contains the fields used below.

    req_model = {'M','D','K','Tau','Gamma','J','D_tilde','n_gen','M_inv','N','F'};
    req_resp  = {'Vr','Vl','eigvals','bstep','t','coi_hz'};
    req_modal = {'mode_step_gain_load','dangerous_complex_modes'};

    for k = 1:numel(req_model)
        if ~isfield(eval_out.model, req_model{k})
            error('Missing field eval_out.model.%s', req_model{k});
        end
    end
    for k = 1:numel(req_resp)
        if ~isfield(eval_out.response, req_resp{k})
            error('Missing field eval_out.response.%s', req_resp{k});
        end
    end
    for k = 1:numel(req_modal)
        if ~isfield(eval_out.modal, req_modal{k})
            error('Missing field eval_out.modal.%s', req_modal{k});
        end
    end
end

%% ========================================================================
function y = local_get_resp_load_hz(eval_out)
% Compatibility helper: the current repository stores resp_load_pu and not
% necessarily resp_load_hz directly.

    if isfield(eval_out.response, 'resp_load_hz')
        y = eval_out.response.resp_load_hz;
    elseif isfield(eval_out.response, 'resp_load_pu')
        base_freq = eval_out.model.mpc.userdata.dynamic.base_frequency_hz;
        y = eval_out.response.resp_load_pu * base_freq;
    else
        error('Neither resp_load_hz nor resp_load_pu was found in eval_out.response.');
    end
end

%% ========================================================================
function modes = local_select_initial_modes(eval_out, n_keep)
% Paper Section 3.1 / 4.1
% Initial selection of dangerous oscillatory modes.
% We only keep the positive-imaginary half to avoid counting conjugate pairs
% twice.

    eigvals = eval_out.response.eigvals;
    gain = eval_out.modal.mode_step_gain_load;

    idx_complex_pos = find(imag(eigvals) > 1e-8);
    if isempty(idx_complex_pos)
        error('No complex oscillatory modes with positive imaginary part were found.');
    end

    score = zeros(numel(idx_complex_pos),1);
    for k = 1:numel(idx_complex_pos)
        score(k) = max(abs(gain(:, idx_complex_pos(k))));
    end

    [~, ord] = sort(score, 'descend');
    take = idx_complex_pos(ord(1:min(n_keep, numel(ord))));

    modes = repmat(local_empty_mode(), numel(take), 1);
    for ii = 1:numel(take)
        modes(ii) = local_build_mode_struct(eval_out, take(ii));
    end
end

%% ========================================================================
function modes_new = local_track_modes_by_mac(eval_out, modes_old, n_keep)
% Paper Section 4.1 implementation detail:
%   Modal Assurance Criterion (MAC) is used so that the tracked oscillatory
%   mode does not get "lost" when eigenvalue ordering changes.
%
% The user explicitly requested that the new code must not hard-code the
% eigenvalue index. This function implements that requirement.

    if isempty(modes_old)
        modes_new = local_select_initial_modes(eval_out, n_keep);
        return;
    end

    eigvals = eval_out.response.eigvals;
    idx_complex_pos = find(imag(eigvals) > 1e-8);
    if isempty(idx_complex_pos)
        modes_new = local_select_initial_modes(eval_out, n_keep);
        return;
    end

    taken = false(numel(idx_complex_pos),1);
    modes_new = repmat(local_empty_mode(), numel(modes_old), 1);

    for ii = 1:numel(modes_old)
        best_score = -Inf;
        best_pos = 0;

        for jj = 1:numel(idx_complex_pos)
            if taken(jj)
                continue;
            end

            cand_idx = idx_complex_pos(jj);
            cand_vr = eval_out.response.Vr(:, cand_idx);
            mac = local_mac(modes_old(ii).vr, cand_vr);

            lam_gap = abs(eigvals(cand_idx) - modes_old(ii).lambda);
            score = mac + 0.05 / (1.0 + lam_gap);

            if score > best_score
                best_score = score;
                best_pos = jj;
            end
        end

        if best_pos == 0
            modes_new(ii) = local_build_mode_struct(eval_out, idx_complex_pos(1));
            taken(1) = true;
        else
            modes_new(ii) = local_build_mode_struct(eval_out, idx_complex_pos(best_pos));
            taken(best_pos) = true;
        end
    end
end

%% ========================================================================
function mode = local_build_mode_struct(eval_out, mode_idx)
    eigvals = eval_out.response.eigvals;
    gain = eval_out.modal.mode_step_gain_load;

    mode = local_empty_mode();
    mode.index = mode_idx;
    mode.lambda = eigvals(mode_idx);
    mode.sigma = real(eigvals(mode_idx));
    mode.omega_d = imag(eigvals(mode_idx));
    mode.vr = eval_out.response.Vr(:, mode_idx);
    mode.vl = eval_out.response.Vl(:, mode_idx);
    mode.max_load_gain = max(abs(gain(:, mode_idx)));
end

%% ========================================================================
function mode = local_empty_mode()
    mode = struct('index', [], 'lambda', [], 'sigma', [], 'omega_d', [], ...
                  'vr', [], 'vl', [], 'max_load_gain', []);
end

%% ========================================================================
function mac = local_mac(v1, v2)
% Standard MAC based on complex inner product.

    den = (norm(v1)^2) * (norm(v2)^2);
    if den < eps
        mac = 0.0;
    else
        mac = abs(v1' * v2)^2 / den;
    end
end

%% ========================================================================
function sens = local_build_mode_sensitivities(eval_out, modes, ctrl_idx)
% Paper Section 3.2 / equations (7)-(9)
% Compute eigenvalue sensitivities of the tracked modes w.r.t. each
% controllable node's virtual inertia M_i and damping D_i.
%
% Implementation choice:
% Instead of hand-coding only the final scalar formula, we build dA/dM_i and
% dA/dD_i directly from the reduced state matrix. This is easier to audit and
% matches the matrix-perturbation principle:
%   d lambda_k / d p = u_k^T (dA/dp) v_k / (u_k^T v_k)

    model = eval_out.model;
    nm = numel(modes);
    nc = numel(ctrl_idx);

    sigma_m = zeros(nc, nm);
    sigma_d = zeros(nc, nm);

    for kk = 1:nm
        mode_idx = modes(kk).index;
        v = eval_out.response.Vr(:, mode_idx);
        u = eval_out.response.Vl(:, mode_idx);
        denom = u.' * v;
        if abs(denom) < 1e-12
            denom = 1.0;
        end

        for ii = 1:nc
            gi = ctrl_idx(ii);

            dA_dm = local_dA_dM(model, gi);
            dA_dd = local_dA_dD(model, gi);

            dlam_dm = (u.' * (dA_dm * v)) / denom;
            dlam_dd = (u.' * (dA_dd * v)) / denom;

            sigma_m(ii,kk) = real(dlam_dm);
            sigma_d(ii,kk) = real(dlam_dd);
        end
    end

    sens = struct();
    sens.sigma_m = sigma_m;
    sens.sigma_d = sigma_d;
end

%% ========================================================================
function dA = local_dA_dM(model, gi)
% Derivative of the reduced state matrix A with respect to M_i.
% Corresponds to paper Section 3.2 (virtual inertia sensitivity).

    n = model.n_gen;
    E = zeros(n);
    E(gi, gi) = 1.0;

    alpha = model.K(gi) * model.Gamma(gi) / (model.M(gi)^2);

    Z = zeros(n);
    dA21 = (1.0 / model.M(gi)^2) * (E * model.J);
    dA22 = (1.0 / model.M(gi)^2) * (E * model.D_tilde);
    dA23 = -(1.0 / model.M(gi)^2) * E;

    dA31 = -alpha * (E * model.J);
    dA32 = -alpha * (E * model.D_tilde);
    dA33 =  alpha * E;

    dA = [Z,   Z,   Z; ...
          dA21, dA22, dA23; ...
          dA31, dA32, dA33];
end

%% ========================================================================
function dA = local_dA_dD(model, gi)
% Derivative of the reduced state matrix A with respect to D_i.
% Corresponds to paper Section 3.2 (virtual damping sensitivity).

    n = model.n_gen;
    E = zeros(n);
    E(gi, gi) = 1.0;

    Z = zeros(n);
    dA22 = -(1.0 / model.M(gi)) * E;
    dA32 =  (model.K(gi) * model.Gamma(gi) / model.M(gi)) * E;

    dA = [Z, Z, Z; ...
          Z, dA22, Z; ...
          Z, dA32, Z];
end

%% ========================================================================
function targets = local_build_security_targets(eval_out, modes, ctrl_idx, opt)
% Paper Section 3.3 / equations (10)-(16)
%
% Build the linearized safety targets used by the inner ADMM:
%   1) RoCoF lower bound -> local lower bound on M_i
%   2) QSS target        -> global lower bound on total damping increment
%   3) Nadir target      -> global upper bound on modal real-part increments
%
% Notes on implementation:
% - The current repository model keeps K fixed. To remain dimensionally
%   consistent with the current reduced model, the QSS damping budget uses
%   (D + K + mu) as the total static restoring term, while only D is allowed
%   to move in the optimization step.
% - The nadir constraint uses the current worst-bus nadir time and the modal
%   envelope amplitude extracted from mode_step_gain_load.

    limits = eval_out.model.ddvpp.security_limits;
    base_freq = eval_out.model.mpc.userdata.dynamic.base_frequency_hz;
    nm = numel(modes);

    % ----- RoCoF bound (paper eq. 10) -----
    rocof_limit_pu = limits.rocof_limit_hz_per_s / base_freq;
    u = eval_out.model.M(:) .* eval_out.response.bstep(eval_out.model.n_gen + (1:eval_out.model.n_gen));
    u = u(:);  % because bstep(second block) = M^{-1} * u

    rocof_m_floor = abs(u) / max(rocof_limit_pu, 1e-12);

    % ----- QSS bound (paper eq. 11) -----
    qss_limit_pu = limits.qss_limit_hz / base_freq;
    disturbance_pu = abs(eval_out.response.disturbance.disturbance_pu);
    total_required_support = disturbance_pu / max(qss_limit_pu, 1e-12);
    total_current_support = sum(eval_out.model.D) + sum(eval_out.model.K) + sum(eval_out.model.mu);
    req_d_sum_raw = max(0.0, total_required_support - total_current_support);

    gtab = eval_out.model.ddvpp.gen_dynamic_table;
    max_feasible_d = sum(max(gtab.d_max(ctrl_idx) - eval_out.model.D(ctrl_idx), 0.0));
    req_d_sum = min(req_d_sum_raw, max_feasible_d);

    % ----- Nadir / modal-damping bound (paper eqs. 12-16) -----
    load_hz = local_get_resp_load_hz(eval_out);
    [worst_nadir_hz, worst_flat_idx] = min(load_hz(:));
    [worst_bus_row, worst_time_col] = ind2sub(size(load_hz), worst_flat_idx);
    t_nadir = eval_out.response.t(worst_time_col);
    t_nadir = max(t_nadir, opt.min_time_s);

    coi_at_nadir = abs(eval_out.response.coi_hz(worst_time_col));
    margin_hz = max(limits.nadir_limit_hz - coi_at_nadir, opt.min_margin_hz);

    req_sigma_shift = zeros(nm,1);
    sigma_safe = zeros(nm,1);

    for kk = 1:nm
        mode_idx = modes(kk).index;
        local_amp = 2.0 * abs(eval_out.modal.mode_step_gain_load(:, mode_idx));
        local_amp = max(local_amp, 1e-12);

        sigma_safe_bus = (1.0 / t_nadir) * log(margin_hz ./ local_amp);
        sigma_safe(kk) = min(sigma_safe_bus);
        req_sigma_shift(kk) = min(0.0, sigma_safe(kk) - real(modes(kk).lambda));
    end

    targets = struct();
    targets.rocof_m_floor = rocof_m_floor(:);
    targets.req_d_sum = req_d_sum;
    targets.req_sigma_shift = req_sigma_shift(:);
    targets.sigma_safe = sigma_safe(:);
    targets.worst_nadir_bus = eval_out.metrics.per_bus.bus(worst_bus_row);
    targets.worst_nadir_hz = worst_nadir_hz;
    targets.t_nadir = t_nadir;
    targets.margin_hz = margin_hz;
    targets.already_secure = local_is_secure(eval_out, limits, opt);
end

%% ========================================================================
function inner = local_run_inner_admm(eval_out, ctrl_idx, sens, targets, trust_radius, rho0, opt)
% Paper Section 4.2 / ADMM Steps 1-3
%
% Global variable definition used here:
%   z = [ Delta sigma_1; ...; Delta sigma_K; sum(Delta D_i) ]
%
% Safe-set projection:
%   z(1:K) <= req_sigma_shift
%   z(K+1) >= req_d_sum
%
% Each controllable node i owns a 2-D local increment:
%   dx_i = [Delta M_i; Delta D_i]
%
% Its local mapping S_i follows the paper's linearization logic:
%   [Delta sigma_1]
%   [Delta sigma_2]   = S_i * dx_i
%   [   ...       ]
%   [Delta sigma_K]
%   [Delta D_sum  ]
%
% where the last row is [0 1], i.e., only Delta D contributes to the QSS
% damping sum in this simplified implementation.

    gtab = eval_out.model.ddvpp.gen_dynamic_table;
    M = eval_out.model.M;
    D = eval_out.model.D;

    nm = size(sens.sigma_m, 2);
    nc = numel(ctrl_idx);
    nz = nm + 1;

    dx = zeros(nc, 2);
    z = zeros(nz, 1);
    y = zeros(nz, 1);

    rho = rho0;

    % Pre-build local linear maps and local box bounds
    S = cell(nc,1);
    lb = zeros(nc,2);
    ub = zeros(nc,2);
    H0 = cell(nc,1);

    for ii = 1:nc
        gi = ctrl_idx(ii);

        S{ii} = [sens.sigma_m(ii,:).', sens.sigma_d(ii,:).'; ...
                 0.0, 1.0];

        m_lb = max([gtab.m_min(gi) - M(gi), ...
                    targets.rocof_m_floor(gi) - M(gi), ...
                   -trust_radius]);
        m_ub = min([gtab.m_max(gi) - M(gi), trust_radius]);

        d_lb = max([gtab.d_min(gi) - D(gi), -trust_radius]);
        d_ub = min([gtab.d_max(gi) - D(gi), trust_radius]);

        if m_lb > m_ub
            m_lb = m_ub;
        end
        if d_lb > d_ub
            d_lb = d_ub;
        end

        lb(ii,:) = [m_lb, d_lb];
        ub(ii,:) = [m_ub, d_ub];

        w_base = max(gtab.local_cost_quad(gi), opt.min_quad_weight);
        w_m = opt.weight_M * w_base / max(M(gi)^2, 1e-8);
        w_d = opt.weight_D * w_base / max(D(gi)^2, 1e-8);
        H0{ii} = diag([w_m, w_d]);
    end

    sum_contrib = zeros(nz,1);
    r_pri = Inf;
    r_dual = Inf;

    for it = 1:opt.max_inner_iter
        % -------------------------------------------------------------
        % Step 1: x-update (node-level local QP)
        % -------------------------------------------------------------
        for ii = 1:nc
            Si = S{ii};
            xi_old = dx(ii,:).';
            r_other = sum_contrib - Si * xi_old;

            target = z - y - r_other;
            H = H0{ii} + rho * (Si' * Si);
            q = -rho * (Si' * target);

            xi_new = local_solve_box_qp_2d(H, q, lb(ii,:).', ub(ii,:).');
            dx(ii,:) = xi_new.';

            sum_contrib = r_other + Si * xi_new;
        end

        % -------------------------------------------------------------
        % Step 2: z-update (projection onto the safe polyhedral set)
        % -------------------------------------------------------------
        z_old = z;
        z = local_project_to_safe_set(sum_contrib + y, targets);

        % -------------------------------------------------------------
        % Step 3: y-update (dual variable / shadow price)
        % -------------------------------------------------------------
        primal_residual = sum_contrib - z;
        y = y + primal_residual;

        r_pri = norm(primal_residual, 2);
        r_dual = rho * norm(z - z_old, 2);

        if r_pri <= opt.primal_tol && r_dual <= opt.dual_tol
            break;
        end

        % -------------------------------------------------------------
        % Residual balancing for rho (user explicitly requested this)
        % -------------------------------------------------------------
        if r_pri > opt.rho_mu * r_dual
            rho_new = min(rho * opt.rho_tau, opt.rho_max);
            if rho_new ~= rho
                y = y / opt.rho_tau;
                rho = rho_new;
            end
        elseif r_dual > opt.rho_mu * r_pri
            rho_new = max(rho / opt.rho_tau, opt.rho_min);
            if rho_new ~= rho
                y = y * opt.rho_tau;
                rho = rho_new;
            end
        end
    end

    inner = struct();
    inner.dx = dx;
    inner.z = z;
    inner.y = y;
    inner.iterations = it;
    inner.rho_final = rho;
    inner.sum_contrib = sum_contrib;
    inner.primal_residual = r_pri;
    inner.dual_residual = r_dual;
end

%% ========================================================================
function z = local_project_to_safe_set(raw, targets)
% Projection onto the polyhedral safe set used by the ADMM z-update.
% For this implementation the set is separable:
%   sigma shifts : upper bounds
%   QSS damping  : lower bound

    nm = numel(targets.req_sigma_shift);
    z = raw;

    for kk = 1:nm
        z(kk) = min(raw(kk), targets.req_sigma_shift(kk));
    end
    z(nm+1) = max(raw(nm+1), targets.req_d_sum);
end

%% ========================================================================
function x = local_solve_box_qp_2d(H, q, lb, ub)
% Solve
%   min 0.5*x'*H*x + q'*x
%   s.t. lb <= x <= ub
%
% The variable dimension is always 2 in this DDVPP implementation, so a
% tiny custom active-set enumeration is enough and avoids any toolbox call.

    candidates = zeros(2, 9);
    n_cand = 0;

    % 1) Unconstrained candidate
    x0 = -H \ q;
    if all(x0 >= lb - 1e-12) && all(x0 <= ub + 1e-12)
        n_cand = n_cand + 1;
        candidates(:, n_cand) = x0;
    end

    % 2-5) One variable fixed, the other free
    for fixed_var = 1:2
        other = 3 - fixed_var;
        for bound_type = 1:2
            x_tmp = zeros(2,1);
            if bound_type == 1
                x_tmp(fixed_var) = lb(fixed_var);
            else
                x_tmp(fixed_var) = ub(fixed_var);
            end

            x_other = -(H(other, fixed_var) * x_tmp(fixed_var) + q(other)) / H(other, other);
            x_tmp(other) = min(max(x_other, lb(other)), ub(other));

            n_cand = n_cand + 1;
            candidates(:, n_cand) = x_tmp;
        end
    end

    % 6-9) Both variables fixed at corners
    corners = [lb(1), lb(2); ...
               lb(1), ub(2); ...
               ub(1), lb(2); ...
               ub(1), ub(2)];
    for kk = 1:size(corners,1)
        n_cand = n_cand + 1;
        candidates(:, n_cand) = corners(kk,:).';
    end

    best_val = Inf;
    best_x = candidates(:,1);
    for kk = 1:n_cand
        xc = candidates(:, kk);
        val = 0.5 * xc.' * H * xc + q.' * xc;
        if val < best_val
            best_val = val;
            best_x = xc;
        end
    end

    x = best_x;
end

%% ========================================================================
function [step_ok, trust_next, rho_tr] = local_trust_region_update(pred_shift, act_shift, trust_cur, opt)
% Paper Section 4.1 / trust-region mechanism requested by the user.
%
% rho_TR compares actual modal-real-part motion against the linearized
% prediction. Using the inner-product ratio preserves sign information:
% if the actual motion is opposite to the predicted left-shift direction,
% rho_TR becomes negative and the step is rejected.

    pred_norm = norm(pred_shift, 2);

    if pred_norm < 1e-12
        rho_tr = 0.0;
    else
        rho_tr = real(pred_shift(:).' * act_shift(:)) / (pred_norm^2);
    end

    step_ok = (rho_tr >= opt.trust_accept_ratio) || (norm(act_shift, 2) < 1e-10);

    if rho_tr > opt.trust_expand_ratio
        trust_next = min(opt.trust_region_max, trust_cur * opt.trust_expand_factor);
    elseif rho_tr < 0.25
        trust_next = max(opt.trust_region_min, trust_cur * opt.trust_shrink_factor);
    else
        trust_next = trust_cur;
    end
end

%% ========================================================================
function act_shift = local_actual_sigma_shift(modes_before, modes_after)
    nm = min(numel(modes_before), numel(modes_after));
    act_shift = zeros(nm,1);
    for kk = 1:nm
        act_shift(kk) = real(modes_after(kk).lambda) - real(modes_before(kk).lambda);
    end
end

%% ========================================================================
function tf = local_is_secure(eval_out, limits, ~)
% Actual nonlinear security check using the time-domain response.
% This is the "real" acceptance criterion after each outer update.

    load_hz = local_get_resp_load_hz(eval_out);

    worst_rocof = max(abs(eval_out.metrics.per_bus.rocof0_hz_per_s));
    worst_nadir = max(abs(load_hz(:)));
    worst_qss = max(abs(load_hz(:, end)));

    tf = (worst_rocof <= limits.rocof_limit_hz_per_s + 1e-6) && ...
         (worst_nadir <= limits.nadir_limit_hz + 1e-6) && ...
         (worst_qss <= limits.qss_limit_hz + 1e-6);
end

%% ========================================================================
function figures = local_plot_before_after(base_eval, final_eval, opt)
% Save before/after figures:
%   1) all-node frequency curves before optimization
%   2) all-node frequency curves after optimization
%   3) worst buses before-vs-after comparison

    if ~exist(opt.save_dir, 'dir')
        mkdir(opt.save_dir);
    end

    t0 = base_eval.response.t(:);
    y0 = local_get_resp_load_hz(base_eval);
    tf = final_eval.response.t(:);
    yf = local_get_resp_load_hz(final_eval);

    T0 = base_eval.metrics.per_bus;
    [~, order0] = sort(T0.nadir_hz, 'ascend');
    topk = min(opt.topk_show, height(T0));
    show_idx = order0(1:topk);

    figures = struct();

    % -------------------------------------------------------------
    % Figure 1: all load-bus frequency curves before optimization
    % -------------------------------------------------------------
    figures.before = figure('Name', 'All-node frequency response BEFORE optimization', 'Color', 'w');
    plot(t0, y0.', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5); hold on;
    h = gobjects(1, topk + 1);
    h(1) = plot(t0, base_eval.response.coi_hz, 'k-', 'LineWidth', 2.0);
    for kk = 1:topk
        h(1+kk) = plot(t0, y0(show_idx(kk),:), 'LineWidth', 1.2);
    end
    yline(-base_eval.model.ddvpp.security_limits.nadir_limit_hz, 'r--', 'LineWidth', 1.0);
    yline( base_eval.model.ddvpp.security_limits.nadir_limit_hz, 'r--', 'LineWidth', 1.0);
    grid on;
    xlabel('Time (s)');
    ylabel('Load-bus frequency deviation (Hz)');
    title('All load-bus frequency responses before optimization');
    legend(h, local_make_curve_legend(topk, base_eval.metrics.per_bus.bus(show_idx)), ...
        'Location', 'eastoutside');
    local_save_figure(figures.before, fullfile(opt.save_dir, 'ddvpp_freq_before_all_nodes.png'));

    % -------------------------------------------------------------
    % Figure 2: all load-bus frequency curves after optimization
    % -------------------------------------------------------------
    figures.after = figure('Name', 'All-node frequency response AFTER optimization', 'Color', 'w');
    plot(tf, yf.', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5); hold on;
    h = gobjects(1, topk + 1);
    h(1) = plot(tf, final_eval.response.coi_hz, 'k-', 'LineWidth', 2.0);
    for kk = 1:topk
        h(1+kk) = plot(tf, yf(show_idx(kk),:), 'LineWidth', 1.2);
    end
    yline(-final_eval.model.ddvpp.security_limits.nadir_limit_hz, 'r--', 'LineWidth', 1.0);
    yline( final_eval.model.ddvpp.security_limits.nadir_limit_hz, 'r--', 'LineWidth', 1.0);
    grid on;
    xlabel('Time (s)');
    ylabel('Load-bus frequency deviation (Hz)');
    title('All load-bus frequency responses after optimization');
    legend(h, local_make_curve_legend(topk, final_eval.metrics.per_bus.bus(show_idx)), ...
        'Location', 'eastoutside');
    local_save_figure(figures.after, fullfile(opt.save_dir, 'ddvpp_freq_after_all_nodes.png'));

    % -------------------------------------------------------------
    % Figure 3: top-k worst buses before vs after
    % -------------------------------------------------------------
    figures.compare = figure('Name', 'Worst-bus comparison before vs after', 'Color', 'w');
    nrow = ceil(topk / 2);
    for kk = 1:topk
        subplot(nrow, 2, kk); %#ok<LAXES>
        bus_id = base_eval.metrics.per_bus.bus(show_idx(kk));
        plot(t0, y0(show_idx(kk),:), 'LineWidth', 1.4); hold on;
        row_after = find(final_eval.metrics.per_bus.bus == bus_id, 1, 'first');
        if isempty(row_after)
            row_after = show_idx(kk);
        end
        plot(tf, yf(row_after,:), '--', 'LineWidth', 1.4);
        yline(-base_eval.model.ddvpp.security_limits.nadir_limit_hz, 'r--');
        grid on;
        title(sprintf('Bus %d', bus_id));
        xlabel('Time (s)');
        ylabel('\Delta f (Hz)');
        legend({'before', 'after', 'limit'}, 'Location', 'best');
    end
    local_save_figure(figures.compare, fullfile(opt.save_dir, 'ddvpp_freq_compare_worst_buses.png'));
end

%% ========================================================================
function entries = local_make_curve_legend(topk, bus_ids)
    entries = cell(1, 1 + topk);
    entries{1} = 'COI';
    for kk = 1:topk
        entries{1+kk} = sprintf('highlight bus %d', bus_ids(kk));
    end
end

%% ========================================================================
function local_save_figure(fig, filename)
    try
        exportgraphics(fig, filename, 'Resolution', 250);
    catch
        saveas(fig, filename);
    end
end

%% ========================================================================
function solution = local_build_solution_table(base_eval, final_eval, ctrl_idx)
% Build a readable node-wise summary for the command window.

    gtab = final_eval.model.ddvpp.gen_dynamic_table;

    bus = gtab.host_bus(ctrl_idx);
    M0 = base_eval.model.M(ctrl_idx);
    D0 = base_eval.model.D(ctrl_idx);
    M1 = final_eval.model.M(ctrl_idx);
    D1 = final_eval.model.D(ctrl_idx);

    solution = table( ...
        ctrl_idx(:), ...
        bus(:), ...
        M0(:), M1(:), (M1(:)-M0(:)), ...
        D0(:), D1(:), (D1(:)-D0(:)), ...
        'VariableNames', {'gen_row','host_bus','M_before','M_after','Delta_M', ...
                          'D_before','D_after','Delta_D'} ...
        );

    [~, ord] = sort(abs(solution.Delta_M) + abs(solution.Delta_D), 'descend');
    solution = solution(ord,:);
end

%% ========================================================================
function local_print_screening(title_str, eval_out, limits)
    T = eval_out.metrics.per_bus;
    load_hz = local_get_resp_load_hz(eval_out);

    [worst_nadir_hz, worst_flat_idx] = min(load_hz(:));
    [worst_bus_row, ~] = ind2sub(size(load_hz), worst_flat_idx);
    [worst_abs_rocof, idx_rocof] = max(abs(T.rocof0_hz_per_s));
    [worst_abs_qss, idx_qss] = max(abs(load_hz(:, end)));

    fprintf('\n[%s]\n', title_str);
    fprintf('  Worst nadir bus        : %d, %.4f Hz\n', T.bus(worst_bus_row), worst_nadir_hz);
    fprintf('  Worst |RoCoF| bus      : %d, %.4f Hz/s\n', T.bus(idx_rocof), worst_abs_rocof);
    fprintf('  Worst |QSS| bus        : %d, %.4f Hz\n', T.bus(idx_qss), worst_abs_qss);
    fprintf('  COI nadir              : %.4f Hz @ %.4f s\n', ...
        eval_out.response.coi_nadir_hz, eval_out.response.coi_nadir_time_s);
    fprintf('  Security check         : RoCoF<=%.3f, |nadir|<=%.3f, |QSS|<=%.3f\n', ...
        limits.rocof_limit_hz_per_s, limits.nadir_limit_hz, limits.qss_limit_hz);
end

%% ========================================================================
function local_print_modes(modes, targets)
    fprintf('Tracked critical modes and required left shift:\n');
    for kk = 1:numel(modes)
        fprintf('  mode %d : idx=%d, lambda=%.6f %+.6fj, req Delta sigma <= %.6f\n', ...
            kk, modes(kk).index, real(modes(kk).lambda), imag(modes(kk).lambda), ...
            targets.req_sigma_shift(kk));
    end
    fprintf('  worst-nadir bus        : %d\n', targets.worst_nadir_bus);
    fprintf('  current nadir          : %.6f Hz\n', targets.worst_nadir_hz);
    fprintf('  nadir time used        : %.6f s\n', targets.t_nadir);
    fprintf('  COI margin for local oscillation envelope : %.6f Hz\n', targets.margin_hz);
    fprintf('  required sum(Delta D)  : %.6f\n', targets.req_d_sum);
end

%% ========================================================================
function local_print_final_report(final_eval, limits, solution, stop_reason, opt)
    fprintf('\n===============================================================\n');
    fprintf('DDVPP SLP-ADMM design finished\n');
    fprintf('Stop reason             : %s\n', stop_reason);
    local_print_screening('Final', final_eval, limits);

    fprintf('\nTop controllable nodes with the largest parameter updates:\n');
    disp(solution(1:min(opt.topk_show, height(solution)), :));

    fprintf('Saved figures:\n');
    fprintf('  %s\n', fullfile(opt.save_dir, 'ddvpp_freq_before_all_nodes.png'));
    fprintf('  %s\n', fullfile(opt.save_dir, 'ddvpp_freq_after_all_nodes.png'));
    fprintf('  %s\n', fullfile(opt.save_dir, 'ddvpp_freq_compare_worst_buses.png'));
    fprintf('===============================================================\n');
end

%% ========================================================================
function row = local_empty_history_row()
    row = struct( ...
        'outer_iter', [], ...
        'trust_radius_before', [], ...
        'trust_ratio', [], ...
        'inner_iterations', [], ...
        'rho_after', [], ...
        'pred_shift_norm', [], ...
        'act_shift_norm', [], ...
        'sigma_req_min', [], ...
        'worst_nadir_before', [], ...
        'worst_nadir_after', [], ...
        'qss_abs_before', [], ...
        'qss_abs_after', [], ...
        'accepted', []);
end

%% ========================================================================
function out = local_merge_struct(base, override)
    out = base;
    if isempty(override)
        return;
    end
    fn = fieldnames(override);
    for k = 1:numel(fn)
        out.(fn{k}) = override.(fn{k});
    end
end
