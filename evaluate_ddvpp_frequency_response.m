function eval = evaluate_ddvpp_frequency_response(params)
% DDVPP-oriented unified evaluator for the reduced nodal-frequency model.
% Typical use
%   eval = evaluate_ddvpp_frequency_response();
%
%   p = struct();
%   p.case_function = 'modified_case118';
%   base = evaluate_ddvpp_frequency_response(p);
%   p.override_M = base.model.M;
%   p.override_D = base.model.D;
%   idx = base.model.ddvpp.gen_dynamic_table.is_controllable;
%   p.override_M(idx) = 1.10 * p.override_M(idx);
%   eval2 = evaluate_ddvpp_frequency_response(p);
%
% Input fields in params
%   case_function           : function handle or case name string
%   duration_s              : simulation horizon
%   n_points                : number of time samples
%
% Disturbance selection
%   disturbance_mode        : 'table_id', 'load_bus', 'gen_bus_equiv'
%   disturbance_id          : row index in ddvpp.disturbance_set
%   disturbance_load_bus    : load-side bus index
%   disturbance_mw          : disturbance size in MW
%   disturbance_bus         : generator host bus for gen_bus_equiv mode
%   disturbance_pu          : equivalent generator-side disturbance in p.u.
%
% Dynamic-parameter overrides
%   override_M              : replace node-wise M
%   override_D              : replace node-wise D
%   override_K              : replace node-wise K
%   override_Tau            : replace node-wise Tau
%   override_Gamma          : replace node-wise Gamma
%   override_mu             : replace load-side mu
%
% Bound control
%   enforce_local_bounds    : if true, clip overrides to local bounds
%
% Output
%   eval.cfg                : resolved evaluator settings
%   eval.model              : network / dynamic matrices and DDVPP tables
%   eval.response           : time-domain response quantities
%   eval.modal              : eigenstructure and modal gains
%   eval.metrics            : per-bus and summary metrics

    if nargin < 1 || isempty(params)
        params = struct();
    end

    cfg = default_eval_config();
    cfg = merge_struct(cfg, params);

    model = build_eval_model(cfg);
    response = simulate_eval_model(model, cfg);
    modal = build_modal_outputs(model, response, cfg);
    metrics = build_eval_metrics(model, response, modal);

    eval = struct();
    eval.cfg = cfg;
    eval.model = model;
    eval.response = response;
    eval.modal = modal;
    eval.metrics = metrics;
end

%% ========================================================================
function cfg = default_eval_config()
% Default evaluator settings.
% The case file provides most static data, while this struct only stores
% evaluator-level choices such as simulation length and disturbance mode.
    cfg = struct();

    % Use the new DDVPP-ready case file by default.
    cfg.case_function = 'modified_case118';

    % Time-domain simulation settings.
    cfg.duration_s = [];
    cfg.n_points = [];

    % Disturbance selection.
    cfg.disturbance_mode = 'table_id';
    cfg.disturbance_id = 1;
    cfg.disturbance_load_bus = [];
    cfg.disturbance_mw = [];
    cfg.disturbance_bus = [];
    cfg.disturbance_pu = [];

    % Dynamic-parameter overrides.
    cfg.override_M = [];
    cfg.override_D = [];
    cfg.override_K = [];
    cfg.override_Tau = [];
    cfg.override_Gamma = [];
    cfg.override_mu = [];

    % If true, overrides are clipped to local bounds stored in the DDVPP case.
    cfg.enforce_local_bounds = false;
end

%% ========================================================================
function model = build_eval_model(cfg)
% Read case data and construct the reduced nodal-frequency matrices.
% Unlike the older version, this function does not regenerate dynamic priors
% from heuristic templates unless needed. It reads node-wise priors directly
% from mpc.userdata.ddvpp.gen_dynamic_table.

    mpc = load_case_from_cfg(cfg);

    if ~isfield(mpc, 'userdata') || ~isfield(mpc.userdata, 'ddvpp')
        error('The selected case file does not contain mpc.userdata.ddvpp.');
    end
    if ~isfield(mpc.userdata.ddvpp, 'gen_dynamic_table')
        error('The selected case file does not contain ddvpp.gen_dynamic_table.');
    end

    ddvpp = mpc.userdata.ddvpp;
    gtab = ddvpp.gen_dynamic_table;

    bus = mpc.bus;
    gen = mpc.gen;
    branch = mpc.branch;
    host_buses = gen(:,1);

    n_load = size(bus,1);
    n_gen = size(gen,1);
    n_total = n_load + n_gen;

    if height(gtab) ~= n_gen
        error('gen_dynamic_table height (%d) does not match number of generators (%d).', height(gtab), n_gen);
    end

    % Time-domain defaults fall back to case metadata if not explicitly given.
    if isempty(cfg.duration_s)
        cfg.duration_s = mpc.userdata.dynamic.duration_s;
    end
    if isempty(cfg.n_points)
        cfg.n_points = mpc.userdata.dynamic.n_points;
    end

    % Build extended susceptance matrix.
    [B_full, transformer_x] = build_network_from_case(mpc, gtab);

    G = (n_load+1):n_total;
    L = 1:n_load;
    B_GG = B_full(G,G);
    B_GL = B_full(G,L);
    B_LG = B_full(L,G);
    B_LL = B_full(L,L);
    B_LL_inv = inv(B_LL);

    J = B_GG - B_GL * B_LL_inv * B_LG;
    F = -B_LL_inv * B_LG;
    Ldist = B_GL * B_LL_inv;

    % Read node-wise priors directly from DDVPP table.
    M = gtab.M0;
    D = gtab.D0;
    K = gtab.K0;
    Tau = gtab.Tau0;
    Gamma = gtab.Gamma0;

    % Load damping remains on the load side.
    if isfield(mpc.userdata.dynamic, 'mu')
        mu = mpc.userdata.dynamic.mu;
    else
        error('mpc.userdata.dynamic.mu is required.');
    end

    % Apply direct overrides, which is the key interface for later DDVPP design.
    if ~isempty(cfg.override_M),     M = cfg.override_M(:); end
    if ~isempty(cfg.override_D),     D = cfg.override_D(:); end
    if ~isempty(cfg.override_K),     K = cfg.override_K(:); end
    if ~isempty(cfg.override_Tau),   Tau = cfg.override_Tau(:); end
    if ~isempty(cfg.override_Gamma), Gamma = cfg.override_Gamma(:); end
    if ~isempty(cfg.override_mu),    mu = cfg.override_mu(:); end

    validate_size(M, n_gen, 'M');
    validate_size(D, n_gen, 'D');
    validate_size(K, n_gen, 'K');
    validate_size(Tau, n_gen, 'Tau');
    validate_size(Gamma, n_gen, 'Gamma');
    validate_size(mu, n_load, 'mu');

    if cfg.enforce_local_bounds
        [M, D, K, Tau, Gamma] = clip_to_local_bounds(M, D, K, Tau, Gamma, gtab);
    end

    D_tilde = diag(D) - Ldist * diag(mu) * F;
    M_inv = diag(1 ./ M);
    N = diag(K .* Gamma ./ M);

    omega0 = mpc.userdata.dynamic.base_frequency_hz;
    A = [zeros(n_gen), omega0 * eye(n_gen), zeros(n_gen); ...
         -M_inv * J, -M_inv * D_tilde, M_inv; ...
          N * J, N * D_tilde - diag(K ./ Tau), -N - diag(1 ./ Tau)];

    % Store resolved values back into a model struct for downstream use.
    model = struct();
    model.cfg = cfg;
    model.mpc = mpc;
    model.ddvpp = ddvpp;
    model.bus = bus;
    model.gen = gen;
    model.branch = branch;
    model.host_buses = host_buses;
    model.n_load = n_load;
    model.n_gen = n_gen;
    model.n_total = n_total;
    model.transformer_x = transformer_x;

    model.B_full = B_full;
    model.B_GG = B_GG;
    model.B_GL = B_GL;
    model.B_LG = B_LG;
    model.B_LL = B_LL;
    model.B_LL_inv = B_LL_inv;

    model.J = J;
    model.F = F;
    model.Ldist = Ldist;

    model.M = M;
    model.D = D;
    model.K = K;
    model.Tau = Tau;
    model.Gamma = Gamma;
    model.mu = mu;
    model.D_tilde = D_tilde;
    model.M_inv = M_inv;
    model.N = N;
    model.A = A;
end

%% ========================================================================
function response = simulate_eval_model(model, cfg)
% Compute the step response for a selected disturbance.
% The dynamic model is linear and time invariant, so modal superposition is
% used directly after eigendecomposition.

    [u, disturbance_meta] = build_disturbance_input(model, cfg);

    bstep = [zeros(model.n_gen,1); model.M_inv * u; -model.N * u];
    [Vr, Lambda] = eig(model.A);
    eigvals = diag(Lambda);
    Vl = inv(Vr).';

    coeff_r = Vr \ bstep;

    Cg = [zeros(model.n_gen), eye(model.n_gen), zeros(model.n_gen)];
    Cl = model.F * Cg;
    Ct = [eye(model.n_gen), zeros(model.n_gen), zeros(model.n_gen)];

    Xg = Cg * Vr;
    Xl = Cl * Vr;
    Xt = Ct * Vr;

    t = linspace(0.0, model.cfg.duration_s, model.cfg.n_points);
    modal_step = zeros(numel(eigvals), numel(t));
    for k = 1:numel(eigvals)
        lam = eigvals(k);
        if abs(lam) < 1e-10
            modal_step(k,:) = t;
        else
            modal_step(k,:) = (exp(lam * t) - 1.0) / lam;
        end
    end

    resp_gen_pu = real(Xg * (modal_step .* coeff_r));
    resp_load_pu = real(Xl * (modal_step .* coeff_r));
    theta_gen = real(Xt * (modal_step .* coeff_r));

    coi_pu = sum(model.M .* resp_gen_pu, 1) / sum(model.M);

    rocof0_gen_pu_per_s = model.M_inv * u;
    rocof0_load_pu_per_s = model.F * rocof0_gen_pu_per_s;

    [nadir_load_pu, nadir_idx] = min(resp_load_pu, [], 2);
    nadir_time_s = t(nadir_idx).';
    [nadir_gen_pu, nadir_gen_idx] = min(resp_gen_pu, [], 2);
    nadir_time_gen_s = t(nadir_gen_idx).';

    [coi_nadir_pu, coi_nadir_idx] = min(coi_pu);
    coi_nadir_time_s = t(coi_nadir_idx);

    omega0 = model.mpc.userdata.dynamic.base_frequency_hz;

    response = struct();
    response.t = t;
    response.disturbance = disturbance_meta;
    response.bstep = bstep;

    response.Vr = Vr;
    response.Vl = Vl;
    response.eigvals = eigvals;
    response.coeff_r = coeff_r;

    response.theta_gen = theta_gen;
    response.resp_gen_pu = resp_gen_pu;
    response.resp_load_pu = resp_load_pu;

    response.coi_pu = coi_pu;
    response.coi_hz = coi_pu * omega0;

    response.rocof0_gen_pu_per_s = rocof0_gen_pu_per_s;
    response.rocof0_load_pu_per_s = rocof0_load_pu_per_s;
    response.rocof0_gen_hz_per_s = rocof0_gen_pu_per_s * omega0;
    response.rocof0_load_hz_per_s = rocof0_load_pu_per_s * omega0;

    response.nadir_load_pu = nadir_load_pu;
    response.nadir_load_hz = nadir_load_pu * omega0;
    response.nadir_time_s = nadir_time_s;

    response.nadir_gen_pu = nadir_gen_pu;
    response.nadir_gen_hz = nadir_gen_pu * omega0;
    response.nadir_time_gen_s = nadir_time_gen_s;

    response.coi_nadir_pu = coi_nadir_pu;
    response.coi_nadir_hz = coi_nadir_pu * omega0;
    response.coi_nadir_time_s = coi_nadir_time_s;
end

%% ========================================================================
function modal = build_modal_outputs(model, response, ~)
% Collect modal quantities used by later DDVPP stages.
% This function does not yet implement inter-iteration mode tracking, but it
% already provides the basic modal objects needed for the next step:
%   - eigenvalues
%   - left/right eigenvectors
%   - load-side modal residues / gains
%   - a simple ranking of dangerous complex modes

    eigvals = response.eigvals;
    Vr = response.Vr;
    Vl = response.Vl;
    coeff_r = response.coeff_r;

    Cg = [zeros(model.n_gen), eye(model.n_gen), zeros(model.n_gen)];
    Cl = model.F * Cg;

    Xg = Cg * Vr;
    Xl = Cl * Vr;

    n_modes = numel(eigvals);
    n_load = model.n_load;
    n_gen = model.n_gen;

    residue_gen = zeros(n_gen, n_modes);
    residue_load = zeros(n_load, n_modes);
    mode_step_gain_gen = zeros(n_gen, n_modes);
    mode_step_gain_load = zeros(n_load, n_modes);

    for k = 1:n_modes
        residue_gen(:,k) = Xg(:,k) * coeff_r(k);
        residue_load(:,k) = Xl(:,k) * coeff_r(k);

        if abs(eigvals(k)) < 1e-10
            mode_step_gain_gen(:,k) = Xg(:,k) * coeff_r(k);
            mode_step_gain_load(:,k) = Xl(:,k) * coeff_r(k);
        else
            mode_step_gain_gen(:,k) = Xg(:,k) * coeff_r(k) / eigvals(k);
            mode_step_gain_load(:,k) = Xl(:,k) * coeff_r(k) / eigvals(k);
        end
    end

    idx_real = find(abs(imag(eigvals)) < 1e-10);
    idx_complex = find(abs(imag(eigvals)) >= 1e-10);

    if isempty(idx_real)
        global_mode_amplitude = zeros(model.n_load,1);
    else
        global_mode_amplitude = vecnorm(mode_step_gain_load(:,idx_real), 2, 2);
    end
    if isempty(idx_complex)
        local_mode_amplitude = zeros(model.n_load,1);
    else
        local_mode_amplitude = vecnorm(mode_step_gain_load(:,idx_complex), 2, 2);
    end

    complex_mode_score = zeros(numel(idx_complex),1);
    for kk = 1:numel(idx_complex)
        k = idx_complex(kk);
        complex_mode_score(kk) = max(abs(mode_step_gain_load(:,k)));
    end
    [~, ord] = sort(complex_mode_score, 'descend');

    dangerous_complex_modes = table();
    if ~isempty(idx_complex)
        sel = idx_complex(ord);
        dangerous_complex_modes = table( ...
            sel(:), ...
            eigvals(sel), ...
            real(eigvals(sel)), ...
            imag(eigvals(sel)), ...
            complex_mode_score(ord), ...
            'VariableNames', {'mode_index','lambda','sigma','omega_d','max_load_mode_gain'});
    end

    modal = struct();
    modal.eigvals = eigvals;
    modal.Vr = Vr;
    modal.Vl = Vl;
    modal.idx_real = idx_real;
    modal.idx_complex = idx_complex;
    modal.residue_gen = residue_gen;
    modal.residue_load = residue_load;
    modal.mode_step_gain_gen = mode_step_gain_gen;
    modal.mode_step_gain_load = mode_step_gain_load;
    modal.global_mode_amplitude = global_mode_amplitude;
    modal.local_mode_amplitude = local_mode_amplitude;
    modal.dangerous_complex_modes = dangerous_complex_modes;
end

%% ========================================================================
function metrics = build_eval_metrics(model, response, modal)
% Summarize the most useful quantities for quick screening:
%   - nadir
%   - nadir time
%   - RoCoF(0+)
%   - local modal indicator
%
% Also attach the generator-side DDVPP table so later scripts can match
% bus-level responses with controllable nodes.

    bus = (1:model.n_load).';
    nadir_hz = response.nadir_load_hz;
    nadir_time_s = response.nadir_time_s;
    rocof0_hz_per_s = response.rocof0_load_hz_per_s;
    local_indicator = real(modal.local_mode_amplitude);
    global_indicator = real(modal.global_mode_amplitude);

    per_bus = table( ...
        bus, ...
        nadir_hz, ...
        nadir_time_s, ...
        rocof0_hz_per_s, ...
        abs(nadir_hz), ...
        abs(rocof0_hz_per_s), ...
        local_indicator, ...
        global_indicator, ...
        'VariableNames', {'bus','nadir_hz','nadir_time_s','rocof0_hz_per_s','abs_nadir_hz','abs_rocof0_hz_per_s','local_indicator','global_indicator'} ...
        );

    per_bus.nadir_rank = dense_rank(per_bus.nadir_hz, 'ascend');
    per_bus.rocof_rank = dense_rank(per_bus.abs_rocof0_hz_per_s, 'descend');
    per_bus.local_rank = dense_rank(per_bus.local_indicator, 'descend');

    [~, worst_nadir_row] = min(per_bus.nadir_hz);
    [~, worst_rocof_row] = max(per_bus.abs_rocof0_hz_per_s);
    [~, worst_local_row] = max(per_bus.local_indicator);

    summary = struct();
    summary.coi_nadir_hz = response.coi_nadir_hz;
    summary.coi_nadir_time_s = response.coi_nadir_time_s;
    summary.worst_nadir_bus = per_bus.bus(worst_nadir_row);
    summary.worst_nadir_hz = per_bus.nadir_hz(worst_nadir_row);
    summary.worst_rocof_bus = per_bus.bus(worst_rocof_row);
    summary.worst_rocof_hz_per_s = per_bus.rocof0_hz_per_s(worst_rocof_row);
    summary.worst_local_bus = per_bus.bus(worst_local_row);
    summary.worst_local_indicator = per_bus.local_indicator(worst_local_row);

    metrics = struct();
    metrics.per_bus = per_bus;
    metrics.summary = summary;
    metrics.gen_dynamic_table = model.ddvpp.gen_dynamic_table;
end

%% ========================================================================
function [u, meta] = build_disturbance_input(model, cfg)
% Convert a selected disturbance description into the equivalent generation-
% side input vector u used by the reduced model.
%
% Three modes are supported:
%   1) table_id       : read one candidate disturbance from ddvpp.disturbance_set
%   2) load_bus       : user specifies load-side bus and MW size
%   3) gen_bus_equiv  : user directly specifies an equivalent generation-side
%                       disturbance at one generator host bus

    meta = struct();
    dset = model.ddvpp.disturbance_set;

    switch lower(cfg.disturbance_mode)
        case 'table_id'
            if isempty(dset)
                error('ddvpp.disturbance_set is empty, so table_id mode cannot be used.');
            end
            rid = cfg.disturbance_id;
            if rid < 1 || rid > height(dset)
                error('disturbance_id is out of range.');
            end

            if ~dset.enabled(rid)
                warning('Selected disturbance row is marked disabled.');
            end

            side = char(dset.side{rid});
            if strcmpi(side, 'load')
                disturbance_pu = dset.deltaP_pu(rid);
                bus_id = dset.bus(rid);
                u = model.Ldist(:, bus_id) * disturbance_pu;

                meta.type = 'load_bus';
                meta.source = 'ddvpp.disturbance_set';
                meta.row = rid;
                meta.bus = bus_id;
                meta.disturbance_pu = disturbance_pu;
                meta.disturbance_mw = dset.deltaP_mw(rid);
            else
                error('Only load-side entries are currently implemented in table_id mode.');
            end

        case 'load_bus'
            if isempty(cfg.disturbance_load_bus) || isempty(cfg.disturbance_mw)
                error('disturbance_load_bus and disturbance_mw are required in load_bus mode.');
            end
            disturbance_pu = cfg.disturbance_mw / model.mpc.baseMVA;
            u = model.Ldist(:, cfg.disturbance_load_bus) * disturbance_pu;

            meta.type = 'load_bus';
            meta.source = 'direct_cfg';
            meta.bus = cfg.disturbance_load_bus;
            meta.disturbance_pu = disturbance_pu;
            meta.disturbance_mw = cfg.disturbance_mw;

        case 'gen_bus_equiv'
            if isempty(cfg.disturbance_bus) || isempty(cfg.disturbance_pu)
                error('disturbance_bus and disturbance_pu are required in gen_bus_equiv mode.');
            end
            u = zeros(model.n_gen,1);
            gen_idx = find(model.host_buses == cfg.disturbance_bus, 1, 'first');
            if isempty(gen_idx)
                error('disturbance_bus does not match any generator host bus.');
            end
            u(gen_idx) = cfg.disturbance_pu;

            meta.type = 'gen_bus_equiv';
            meta.source = 'direct_cfg';
            meta.bus = cfg.disturbance_bus;
            meta.disturbance_pu = cfg.disturbance_pu;
            meta.disturbance_mw = cfg.disturbance_pu * model.mpc.baseMVA;

        otherwise
            error('Unsupported disturbance_mode.');
    end
end

%% ========================================================================
function mpc = load_case_from_cfg(cfg)
    if isa(cfg.case_function, 'function_handle')
        mpc = cfg.case_function();
    elseif ischar(cfg.case_function) || isstring(cfg.case_function)
        mpc = feval(char(cfg.case_function));
    else
        error('cfg.case_function must be a function handle or case name string.');
    end
end

%% ========================================================================
function [B, transformer_x] = build_network_from_case(mpc, gtab)
% Build the extended susceptance matrix.
% The original IEEE-118 load-side network comes from mpc.branch.
% Each generator is then connected to its host bus through the effective
% coupling reactance stored in gen_dynamic_table.xdp_effective_initial.

    bus = mpc.bus;
    gen = mpc.gen;
    branch = mpc.branch;

    n_load = size(bus,1);
    n_gen = size(gen,1);
    n_total = n_load + n_gen;
    B = zeros(n_total, n_total);

    line_x_scale = mpc.userdata.dynamic.line_x_scale;
    north_edges = mpc.userdata.dynamic.north_edges;
    north_corridor_scale = mpc.userdata.dynamic.north_corridor_scale;

    for rr = 1:size(branch,1)
        f = branch(rr,1);
        t = branch(rr,2);
        x = branch(rr,4) * line_x_scale;

        tap = branch(rr,9);
        if tap == 0
            tap = 1.0;
        end

        b = 1.0 / (x * tap);
        if is_weak_edge(f, t, north_edges)
            b = b * north_corridor_scale;
        end

        B(f,f) = B(f,f) + b;
        B(t,t) = B(t,t) + b;
        B(f,t) = B(f,t) - b;
        B(t,f) = B(t,f) - b;
    end

    transformer_x = gtab.xdp_effective_initial;
    host_buses = gen(:,1);

    for k = 1:n_gen
        gi = n_load + k;
        li = host_buses(k);
        x = transformer_x(k);
        b = 1.0 / x;

        B(gi,gi) = B(gi,gi) + b;
        B(li,li) = B(li,li) + b;
        B(gi,li) = B(gi,li) - b;
        B(li,gi) = B(li,gi) - b;
    end
end

%% ========================================================================
function [M, D, K, Tau, Gamma] = clip_to_local_bounds(M, D, K, Tau, Gamma, gtab)
% Clip node-wise parameters to the bounds stored in the DDVPP case file.
% This is optional but useful for quick tests before a proper optimizer is
% implemented.

    M = min(max(M, gtab.m_min), gtab.m_max);
    D = min(max(D, gtab.d_min), gtab.d_max);
    K = min(max(K, gtab.k_min), gtab.k_max);
    Tau = min(max(Tau, gtab.tau_min), gtab.tau_max);
    Gamma = min(max(Gamma, gtab.gamma_min), gtab.gamma_max);
end

%% ========================================================================
function validate_size(x, n, name)
    if numel(x) ~= n
        error('%s must have length %d, but got %d.', name, n, numel(x));
    end
end

%% ========================================================================
function tf = is_weak_edge(f, t, edges)
    pair = sort([f t]);
    tf = any(edges(:,1) == pair(1) & edges(:,2) == pair(2));
end

%% ========================================================================
function r = dense_rank(x, direction)
    if nargin < 2
        direction = 'ascend';
    end

    if strcmpi(direction, 'descend')
        vals = unique(x(:), 'sorted');
        vals = flipud(vals);
    else
        vals = unique(x(:), 'sorted');
    end

    r = zeros(size(x));
    for ii = 1:numel(vals)
        r(x == vals(ii)) = ii;
    end
end

%% ========================================================================
function out = merge_struct(base, override)
    out = base;
    if isempty(override)
        return;
    end
    fn = fieldnames(override);
    for k = 1:numel(fn)
        out.(fn{k}) = override.(fn{k});
    end
end
