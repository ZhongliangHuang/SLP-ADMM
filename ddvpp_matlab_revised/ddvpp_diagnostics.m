function diag = ddvpp_diagnostics(data, x, model, dist, sim, crit, req)
% Diagnostics that make implementation errors easy to spot from copied logs.

    ncheck = min(3, numel(crit.idx));
    sens_relerr = nan(ncheck, 2 * data.nibr);
    eps_fd = 1e-6;
    for k = 1:ncheck
        lam0 = crit.lambda(k);
        for ii = 1:2*data.nibr
            xp = x;
            xm = x;
            xp(ii) = xp(ii) + eps_fd;
            xm(ii) = max(1e-6, xm(ii) - eps_fd);
            mp = ddvpp_linear_model(data, xp);
            mm = ddvpp_linear_model(data, xm);
            lp = local_match_eig(mp.A, lam0);
            lm = local_match_eig(mm.A, lam0);
            fd = real((lp - lm) / (2 * eps_fd));
            den = max(1e-8, abs(fd));
            sens_relerr(k, ii) = abs(fd - crit.S_sigma(k, ii)) / den;
        end
    end

    qss_terminal_hz = max(abs(sim.freq_hz_signed(:, end)));
    rocof0_err = max(abs(sim.rocof_hz_s_signed(:,1) - sim.rocof0_hz_s_signed));
    norm_err = max(abs(crit.normalization - 1));
    zero_eigs = sum(abs(eig(model.A)) < 1e-8);

    diag = struct();
    diag.qss_terminal_hz = qss_terminal_hz;
    diag.qss_req = req.qss_req;
    diag.rocof0_max_err = rocof0_err;
    diag.modal_norm_max_err = norm_err;
    diag.zero_eig_count = zero_eigs;
    diag.max_sensitivity_relerr = max(sens_relerr(:));
    diag.sensitivity_relerr = sens_relerr;
    diag.top_nadir_buses = data.bus_ids(local_topk(sim.nadir_hz, 5));
    diag.top_rocof_buses = data.bus_ids(local_topk(sim.rocof_max_hz_s, 5));
end

function lam = local_match_eig(A, lam_ref)
    vals = eig(A);
    [~, idx] = min(abs(vals - lam_ref));
    lam = vals(idx);
end

function idx = local_topk(x, k)
    [~, order] = sort(x, 'descend');
    idx = order(1:min(k, numel(order)));
end
