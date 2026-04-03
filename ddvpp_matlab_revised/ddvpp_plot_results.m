function ddvpp_plot_results(results, outdir)
    data = results.data;
    sim = results.sim;
    hist = results.history;
    crit = results.critical;

    % 1) Frequency trajectories at worst buses
    [~, idx] = sort(sim.nadir_hz, 'descend');
    topk = idx(1:min(8, numel(idx)));
    f1 = figure('Color', 'w'); hold on; grid on;
    for i = 1:numel(topk)
        plot(sim.t, sim.freq_hz(topk(i), :), 'LineWidth', 1.2);
    end
    yline(data.security.Nadir_lim, 'r--', 'LineWidth', 1.4);
    xlabel('Time (s)'); ylabel('|\Delta f_i| (Hz)');
    title('Critical nodal frequency deviations');
    legend(cellstr(num2str(data.bus_ids(topk))), 'Location', 'best');
    saveas(f1, fullfile(outdir, 'fig1_frequency_trajectories.png'));

    % 2) Spatial nodal nadir bar chart
    f2 = figure('Color', 'w');
    bar(data.bus_ids, sim.nadir_hz); hold on; grid on;
    yline(data.security.Nadir_lim, 'r--', 'LineWidth', 1.4);
    xlabel('Bus'); ylabel('Local nadir (Hz)');
    title('Spatial distribution of nodal nadir');
    saveas(f2, fullfile(outdir, 'fig2_spatial_nadir.png'));

    % 3) Eigenvalue scatter
    f3 = figure('Color', 'w'); hold on; grid on;
    lam = eig(results.model.A);
    scatter(real(lam), imag(lam), 20, 'filled');
    scatter(real(crit.lambda), imag(crit.lambda), 60, 'd', 'filled');
    xlabel('Real part'); ylabel('Imaginary part');
    title('State matrix eigenvalues and tracked critical modes');
    saveas(f3, fullfile(outdir, 'fig3_eigenvalues.png'));

    % 4) Outer loop diagnostics
    f4 = figure('Color', 'w'); hold on; grid on;
    yyaxis left; plot(hist.outer, hist.rhoTR, 'o-', 'LineWidth', 1.5);
    ylabel('\rho_{TR}');
    yyaxis right; plot(hist.outer, hist.delta, 's-', 'LineWidth', 1.5);
    ylabel('\delta');
    xlabel('Outer iteration');
    title('Trust-region adaptation diagnostics');
    saveas(f4, fullfile(outdir, 'fig4_trust_region.png'));

    % 5) ADMM iterations and penalty evolution
    f5 = figure('Color', 'w'); hold on; grid on;
    yyaxis left; plot(hist.outer, hist.admm_iters, 'o-', 'LineWidth', 1.5);
    ylabel('ADMM iterations');
    yyaxis right; plot(hist.outer, hist.rho, 's-', 'LineWidth', 1.5);
    ylabel('\rho');
    xlabel('Outer iteration');
    title('Residual balancing effect on ADMM');
    saveas(f5, fullfile(outdir, 'fig5_admm_balance.png'));

    % 6) Final parameter allocation
    [m_all, d_all] = ddvpp_unpack_x(data, results.x_final);
    f6 = figure('Color', 'w'); hold on; grid on;
    stairs(data.ddvpp_bounds(:,1), m_all(data.ibr_idx), 'LineWidth', 1.5);
    stairs(data.ddvpp_bounds(:,1), d_all(data.ibr_idx), 'LineWidth', 1.5);
    xlabel('IBR bus'); ylabel('Allocated value');
    title('Final DDVPP allocation on IBR nodes');
    legend({'Virtual inertia m', 'Virtual damping d'}, 'Location', 'best');
    saveas(f6, fullfile(outdir, 'fig6_final_allocation.png'));

    % 7) RoCoFmax spatial distribution
    f7 = figure('Color', 'w');
    bar(data.bus_ids, sim.rocof_max_hz_s); hold on; grid on;
    yline(data.security.RoCoF_lim, 'r--', 'LineWidth', 1.4);
    xlabel('Bus'); ylabel('RoCoF_{max} (Hz/s)');
    title('Spatial distribution of RoCoF_{max}');
    saveas(f7, fullfile(outdir, 'fig7_rocof_spatial.png'));

    % 8) Text summary
    fid = fopen(fullfile(outdir, 'summary.txt'), 'w');
    fprintf(fid, 'Worst disturbance bus: %d\n', results.worst_bus);
    fprintf(fid, 'Max nodal nadir (Hz): %.6f\n', results.metrics.max_nadir_hz);
    fprintf(fid, 'Max COI nadir (Hz): %.6f\n', results.metrics.max_coi_nadir_hz);
    fprintf(fid, 'QSS deviation (Hz): %.6f\n', results.metrics.qss_hz);
    fprintf(fid, 'Max RoCoF proxy (Hz/s): %.6f\n', results.metrics.max_rocof_proxy);
    fprintf(fid, 'Max RoCoF from simulation (Hz/s): %.6f\n', results.metrics.max_rocof_sim);
    fprintf(fid, 'Nadir safe: %d\n', results.metrics.safe_nadir);
    fprintf(fid, 'QSS safe: %d\n', results.metrics.safe_qss);
    fprintf(fid, 'RoCoF safe: %d\n', results.metrics.safe_rocof);
    fclose(fid);
end
