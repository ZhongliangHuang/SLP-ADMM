function req = ddvpp_build_requirements(data, x, sim, crit, worst_bus, disturbance_size)
% Build linearized security requirements.
% The oscillatory safety margin uses the rigorous envelope 2*|C_{i,k}|*f_base.

    [m_all, d_all] = ddvpp_unpack_x(data, x);
    f_lim = data.security.Nadir_lim;
    qss_lim = data.security.QSS_lim;
    roc_lim = data.security.RoCoF_lim;

    [max_nadir, i_bus] = max(sim.nadir_hz);
    t_nadir = sim.nadir_t(i_bus);
    [~, t_idx] = min(abs(sim.t - t_nadir));
    coi_nadir = sim.coi_hz_drop(t_idx);
    margin = max(f_lim - coi_nadir, 1e-9);

    K = numel(crit.lambda);
    mode_req = zeros(K,1);
    sigma_safe = nan(K,1);
    active_bus = nan(K,1);
    t_eff = max(t_nadir, data.security.nadir_time_floor);
    for k = 1:K
        env_all = max(crit.envelope_hz(:,k), 1e-12);
        ratio = margin ./ env_all;
        active = env_all >= (1 + data.security.envelope_tol) * margin;
        if any(active)
            sigma_safe_all = log(ratio(active)) ./ t_eff;
            [sigma_safe(k), idx_local] = min(sigma_safe_all);
            buses = find(active);
            active_bus(k) = buses(idx_local);
            mode_req(k) = max(0, real(crit.lambda(k)) - sigma_safe(k) + data.security.mode_margin_eps);
        else
            sigma_safe(k) = real(crit.lambda(k));
            mode_req(k) = 0;
        end
    end

    total_d = sum(d_all) + sum(data.mu_load);
    required_total_d = abs(disturbance_size) * data.f_base / max(qss_lim, 1e-9);
    qss_req = max(0, required_total_d - total_d);

    eqP = abs(ddvpp_make_disturbance(data, worst_bus, disturbance_size));
    roc_need = eqP * data.f_base / max(roc_lim, 1e-9);
    roc_deficit = max(0, roc_need - m_all);
    roc_deficit_ibr = roc_deficit(data.ibr_idx);
    roc_deficit_nonibr = roc_deficit(~data.ibr_mask);
    uncontrollable_roc = any(roc_deficit_nonibr > 1e-8);

    req = struct();
    req.mode_req = mode_req;
    req.qss_req = qss_req;
    req.roc_deficit = roc_deficit;
    req.roc_deficit_ibr = roc_deficit_ibr;
    req.required_total_d = required_total_d;
    req.max_nadir = max_nadir;
    req.margin = margin;
    req.t_eff = t_eff;
    req.sigma_safe = sigma_safe;
    req.active_bus = active_bus;
    req.uncontrollable_roc = uncontrollable_roc;
    req.all_satisfied = (max(mode_req) <= 1e-6) && (qss_req <= 1e-6) && all(roc_deficit_ibr <= 1e-8) && ~uncontrollable_roc;
end
