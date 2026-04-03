function metrics = ddvpp_evaluate_security(data, x, sim, worst_bus, disturbance_size)
    [m_all, d_all] = ddvpp_unpack_x(data, x);
    total_d = sum(d_all) + sum(data.mu_load);
    qss_hz = abs(disturbance_size) * data.f_base / max(total_d, 1e-12);

    eqP = abs(ddvpp_make_disturbance(data, worst_bus, disturbance_size));
    rocof_gen = data.f_base * eqP ./ max(m_all, 1e-12);
    rocof_all = zeros(data.nbus,1);
    rocof_all(data.gen_pos) = rocof_gen;
    rocof_all(data.load_pos) = abs(data.F * rocof_gen);

    metrics = struct();
    metrics.max_nadir_hz = max(sim.nadir_hz);
    metrics.max_coi_nadir_hz = max(sim.coi_hz);
    metrics.qss_hz = qss_hz;
    metrics.max_rocof_proxy = max(abs(rocof_all));
    metrics.max_rocof_sim = max(sim.rocof_max_hz_s);
    metrics.safe_nadir = metrics.max_nadir_hz <= data.security.Nadir_lim;
    metrics.safe_qss = qss_hz <= data.security.QSS_lim;
    metrics.safe_rocof = metrics.max_rocof_proxy <= data.security.RoCoF_lim;
end
