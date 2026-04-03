function results = ddvpp_slp_admm_case118(mpc, opts)
% Core nested SLP-ADMM driver.

    data = ddvpp_build_model_data(mpc);
    x = ddvpp_pack_x(data);
    delta = mpc.admm.delta_init;
    rho = mpc.admm.rho_init;

    [worst_bus, scan_info] = ddvpp_find_worst_disturbance(data, x, opts);
    if opts.verbose
        fprintf('Selected worst disturbance bus = %d\n', worst_bus);
    end

    hist = struct('outer', [], 'cost', [], 'delta', [], 'rho', [], ...
                  'pred_push', [], 'actual_push', [], 'rhoTR', [], ...
                  'signed_mode_idx', [], 'signed_pred', [], 'signed_actual', [], ...
                  'worst_bus', worst_bus, 'admm_iters', []);

    tracked_prev = [];
    accept_count = 0;
    req = [];

    for v = 1:opts.max_outer
        model = ddvpp_linear_model(data, x);
        dist = ddvpp_make_disturbance(data, worst_bus, opts.disturbance_size);
        sim0 = ddvpp_simulate_frequency(data, model, dist, opts.time_horizon, opts.dt, data.f_base);
        crit = ddvpp_identify_critical_modes(model, dist, sim0, opts.n_modes, tracked_prev, data);
        tracked_prev = crit.tracked;

        req = ddvpp_build_requirements(data, x, sim0, crit, worst_bus, opts.disturbance_size);
        if opts.verbose
            fprintf('Outer %2d | max local nadir = %.4f Hz | required mode push max = %.4e | qss deficit = %.4e | max ibr roc deficit = %.4e\n', ...
                v, max(sim0.nadir_hz), max([req.mode_req; 0]), req.qss_req, max([req.roc_deficit_ibr; 0]));
        end

        admm = ddvpp_run_admm(data, x, crit, req, delta, rho, mpc, opts);
        if admm.infeasible
            x_trial = x;
            predicted_push = zeros(size(req.mode_req));
            actual_push = zeros(size(req.mode_req));
            binding_idx = 0;
            pred_signed = 0;
            act_signed = 0;
            rhoTR = -Inf;
            accept = false;
            if opts.verbose
                fprintf('  ADMM declared infeasible: %s\n', admm.message);
            end
        else
            x_trial = x + admm.dx_global;

            model1 = ddvpp_linear_model(data, x_trial);
            sim1 = ddvpp_simulate_frequency(data, model1, dist, opts.time_horizon, opts.dt, data.f_base);
            crit1 = ddvpp_identify_critical_modes(model1, dist, sim1, opts.n_modes, tracked_prev, data);

            predicted_push = sum(crit.S_sigma .* admm.dx_full, 2);
            actual_push = real(crit1.lambda) - real(crit.lambda);

            [binding_idx, pred_signed, act_signed, rhoTR] = signed_trust_ratio(crit, req, predicted_push, actual_push);
            accept = isfinite(rhoTR) && (pred_signed > 0) && (act_signed > 0) && (rhoTR >= 0.10);
        end

        if accept
            x = x_trial;
            accept_count = accept_count + 1;
            tracked_prev = crit1.tracked;
            if rhoTR > 0.75
                delta = min(opts.delta_max, delta * opts.rho_expand);
            elseif rhoTR < 0.25
                delta = max(opts.delta_min, delta * opts.rho_shrink);
            end
        else
            delta = max(opts.delta_min, delta * opts.rho_shrink);
        end

        rho = admm.rho_final;

        hist.outer(end+1,1) = v;
        hist.cost(end+1,1) = admm.cost;
        hist.delta(end+1,1) = delta;
        hist.rho(end+1,1) = rho;
        hist.pred_push(end+1,1) = norm(predicted_push);
        hist.actual_push(end+1,1) = norm(actual_push);
        hist.rhoTR(end+1,1) = rhoTR;
        hist.signed_mode_idx(end+1,1) = binding_idx;
        hist.signed_pred(end+1,1) = pred_signed;
        hist.signed_actual(end+1,1) = act_signed;
        hist.admm_iters(end+1,1) = admm.iters;

        if opts.verbose
            fprintf(['  ADMM iters = %d | binding mode = %d | predicted left shift = %.3e | ', ...
                     'actual left shift = %.3e | rhoTR = %.3f | accept = %d | delta = %.3f | rho = %.3f\n'], ...
                     admm.iters, binding_idx, pred_signed, act_signed, rhoTR, accept, delta, rho);
        end

        if req.all_satisfied && admm.converged && ~admm.infeasible
            if opts.verbose
                fprintf('All controllable constraints satisfied with converged ADMM. Stopping outer loop.\n');
            end
            break;
        end
    end

    modelF = ddvpp_linear_model(data, x);
    distF = ddvpp_make_disturbance(data, worst_bus, opts.disturbance_size);
    simF = ddvpp_simulate_frequency(data, modelF, distF, opts.time_horizon, opts.dt, data.f_base);
    critF = ddvpp_identify_critical_modes(modelF, distF, simF, opts.n_modes, tracked_prev, data);
    metrics = ddvpp_evaluate_security(data, x, simF, worst_bus, opts.disturbance_size);
    diagnostics = ddvpp_diagnostics(data, x, modelF, distF, simF, critF, req);

    results = struct();
    results.data = data;
    results.mpc = mpc;
    results.opts = opts;
    results.x_final = x;
    results.delta_final = delta;
    results.rho_final = rho;
    results.worst_bus = worst_bus;
    results.scan_info = scan_info;
    results.history = hist;
    results.model = modelF;
    results.sim = simF;
    results.critical = critF;
    results.metrics = metrics;
    results.diagnostics = diagnostics;
    results.accept_count = accept_count;

    fprintf('S_sigma max abs = %.3e\n', max(abs(critF.S_sigma(:))));
    nz = abs(critF.S_sigma) > 0;
    if any(nz(:))
        fprintf('S_sigma min abs nonzero = %.3e\n', min(abs(critF.S_sigma(nz))));
    end
    fprintf('mode_req = '); fprintf('%.3e ', req.mode_req); fprintf('\n');
    [~, ib] = max(simF.nadir_hz);
    fprintf('worst nadir bus = %d, t_nadir = %.4f s\n', data.bus_ids(ib), simF.nadir_t(ib));
    fprintf('max sensitivity FD relative error = %.3e\n', diagnostics.max_sensitivity_relerr);
    fprintf('max RoCoF initial-condition mismatch = %.3e Hz/s\n', diagnostics.rocof0_max_err);
end

function [binding_idx, pred_left_shift, act_left_shift, rhoTR] = signed_trust_ratio(crit, req, predicted_push, actual_push)
% A beneficial move is a leftward eigenvalue shift, i.e. negative d sigma.

    [req_max, binding_idx] = max(req.mode_req);
    if isempty(binding_idx) || req_max <= 0
        [~, binding_idx] = max(real(crit.lambda));
    end

    pred_left_shift = -predicted_push(binding_idx);
    act_left_shift = -actual_push(binding_idx);

    if pred_left_shift <= 1e-10
        rhoTR = -Inf;
    else
        rhoTR = act_left_shift / pred_left_shift;
    end
end
