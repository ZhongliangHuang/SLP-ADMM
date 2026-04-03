function sim = ddvpp_simulate_frequency(data, model, dist, T_end, dt, f_base)
% Simulate the reduced model and reconstruct the full-bus frequency.
% A Van Loan block exponential is used so that the singular zero eigenvalue
% does not require pinv(A).

    A = model.A;
    B = model.B;
    n3 = size(A,1);
    n = n3 / 3;
    t = 0:dt:T_end;
    nt = numel(t);
    U = B * dist;
    X = zeros(n3, nt);

    Aug = [A, U; zeros(1, n3), 0];
    for k = 1:nt
        Ek = expm(Aug * t(k));
        X(:,k) = Ek(1:n3, end);
    end

    omegaG_pu = X(n+1:2*n, :);
    omegaL_pu = data.F * omegaG_pu;

    omega_all_pu = zeros(data.nbus, nt);
    omega_all_pu(data.gen_pos, :) = omegaG_pu;
    omega_all_pu(data.load_pos, :) = omegaL_pu;

    freq_all_hz_signed = f_base * omega_all_pu;
    freqG_hz_signed = f_base * omegaG_pu;
    freqL_hz_signed = f_base * omegaL_pu;

    freq_drop_hz = max(0, -freq_all_hz_signed);
    freq_drop_hz_G = max(0, -freqG_hz_signed);
    freq_drop_hz_L = max(0, -freqL_hz_signed);

    [nadir_hz, nadir_idx] = max(freq_drop_hz, [], 2);
    nadir_t = t(nadir_idx(:));
    nadir_signed = zeros(data.nbus,1);
    for i = 1:data.nbus
        nadir_signed(i) = freq_all_hz_signed(i, nadir_idx(i));
    end

    coi_dev_pu = (model.m(:).' * omegaG_pu) / max(sum(model.m), 1e-12);
    coi_hz_signed = f_base * coi_dev_pu(:);
    coi_hz_drop = max(0, -coi_hz_signed);

    rocof_signed = zeros(data.nbus, nt);
    if nt >= 2
        rocof_signed(:,1) = (freq_all_hz_signed(:,2) - freq_all_hz_signed(:,1)) / dt;
        for k = 2:nt-1
            rocof_signed(:,k) = (freq_all_hz_signed(:,k+1) - freq_all_hz_signed(:,k-1)) / (2*dt);
        end
        rocof_signed(:,nt) = (freq_all_hz_signed(:,nt) - freq_all_hz_signed(:,nt-1)) / dt;
    end
    rocof_abs = abs(rocof_signed);
    rocof_max_hz_s = max(rocof_abs, [], 2);

    rocof0_gen = f_base * (dist(:) ./ max(model.m(:), 1e-12));
    rocof0_all = zeros(data.nbus,1);
    rocof0_all(data.gen_pos) = rocof0_gen;
    rocof0_all(data.load_pos) = f_base * (data.F * rocof0_gen);

    sim = struct();
    sim.t = t;
    sim.X = X;
    sim.omegaG = omegaG_pu;
    sim.omegaL = omegaL_pu;
    sim.omega_all = omega_all_pu;
    sim.freq_hz = freq_drop_hz;
    sim.freq_hz_signed = freq_all_hz_signed;
    sim.freq_hz_G = freq_drop_hz_G;
    sim.freq_hz_L = freq_drop_hz_L;
    sim.freq_hz_G_signed = freqG_hz_signed;
    sim.freq_hz_L_signed = freqL_hz_signed;
    sim.nadir_hz = nadir_hz;
    sim.nadir_signed = nadir_signed;
    sim.nadir_t = nadir_t;
    sim.coi_hz = coi_hz_drop;
    sim.coi_hz_drop = coi_hz_drop;
    sim.coi_hz_signed = coi_hz_signed;
    sim.rocof_hz_s_signed = rocof_signed;
    sim.rocof_hz_s = rocof_abs;
    sim.rocof_max_hz_s = rocof_max_hz_s;
    sim.rocof0_hz_s = abs(rocof0_all);
    sim.rocof0_hz_s_signed = rocof0_all;
end
