function crit = ddvpp_identify_critical_modes(model, dist, sim, n_modes, tracked_prev, data)
% Identify critical oscillatory modes and compute rigorous modal envelopes.
% For one representative eigenvalue from each conjugate pair with imag>0,
% the node-wise envelope amplitude is
%   C_{i,k} = | (u_k^H B d) * v_{omega,i}^{(k)} / lambda_k |
% and the real oscillatory contribution envelope is 2*C_{i,k}.

    [Vr, L] = eig(model.A);
    lambda = diag(L);
    [Vl, ~] = eig(model.A.');
    Vl = conj(Vl);

    idx_osc = find(imag(lambda) > 1e-6 & real(lambda) < -1e-10);
    if isempty(idx_osc)
        error('No stable oscillatory modes with positive imaginary part were found.');
    end

    n = data.ngen;
    U = model.B * dist;
    residue_score = zeros(numel(idx_osc),1);
    for kk = 1:numel(idx_osc)
        j = idx_osc(kk);
        v = Vr(:,j);
        u = Vl(:,j) / (Vl(:,j)' * Vr(:,j));
        alpha = u' * U;
        Cg = abs((alpha / lambda(j)) * v(n+1:2*n));
        Cl = abs(data.F * ((alpha / lambda(j)) * v(n+1:2*n)));
        C_all = zeros(data.nbus,1);
        C_all(data.gen_pos) = Cg;
        C_all(data.load_pos) = Cl;
        residue_score(kk) = max(2 * C_all) / max(-real(lambda(j)), 1e-6);
    end

    [~, order0] = sort(residue_score, 'descend');
    candidates = idx_osc(order0);

    tracked_idx = zeros(min(n_modes, numel(candidates)),1);
    if isempty(tracked_prev)
        tracked_idx = candidates(1:numel(tracked_idx));
    else
        used = false(numel(candidates),1);
        for k = 1:numel(tracked_idx)
            prev_v = tracked_prev.Vr(:,k);
            mac = -inf(numel(candidates),1);
            for j = 1:numel(candidates)
                if used(j)
                    continue;
                end
                vv = Vr(:, candidates(j));
                mac(j) = abs(prev_v' * vv)^2 / (((prev_v' * prev_v) * (vv' * vv)) + 1e-12);
            end
            [~, idx_best] = max(real(mac));
            tracked_idx(k) = candidates(idx_best);
            used(idx_best) = true;
        end
    end

    tracked_lambda = lambda(tracked_idx);
    [~, order] = sort(real(tracked_lambda), 'descend');
    tracked_idx = tracked_idx(order);

    K = numel(tracked_idx);
    S_sigma = zeros(K, 2 * data.nibr);
    Cik_abs = zeros(data.nbus, K);
    alpha_vec = zeros(K,1);
    norm_check = zeros(K,1);
    for k = 1:K
        j = tracked_idx(k);
        u = Vl(:,j) / (Vl(:,j)' * Vr(:,j));
        v = Vr(:,j);
        alpha = u' * U;
        alpha_vec(k) = alpha;
        norm_check(k) = Vl(:,j)' * Vr(:,j);
        S_sigma(k,:) = ddvpp_mode_sensitivity(data, model, u, v);

        Cg = abs((alpha / lambda(j)) * v(n+1:2*n));
        Cl = abs(data.F * ((alpha / lambda(j)) * v(n+1:2*n)));
        C_all = zeros(data.nbus,1);
        C_all(data.gen_pos) = Cg;
        C_all(data.load_pos) = Cl;
        Cik_abs(:,k) = C_all;
    end

    crit = struct();
    crit.idx = tracked_idx;
    crit.lambda = lambda(tracked_idx);
    crit.Vr = Vr(:, tracked_idx);
    crit.Vl = Vl(:, tracked_idx);
    crit.tracked = struct('Vr', Vr(:, tracked_idx), 'Vl', Vl(:, tracked_idx));
    crit.S_sigma = S_sigma;
    crit.Cik_abs = Cik_abs;
    crit.residue_score = residue_score;
    crit.sim_nadir = sim.nadir_hz;
    crit.alpha = alpha_vec;
    crit.normalization = norm_check;
    crit.envelope_hz = 2 * Cik_abs * data.f_base;
end
