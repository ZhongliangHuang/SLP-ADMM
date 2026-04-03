function results = run_ddvpp_case118_demo()
% RUN_DDVPP_CASE118_DEMO
% Main entry for the DDVPP case study on IEEE 118-bus system.
% It implements:
%   1) outer SLP loop with adaptive trust-region
%   2) inner ADMM with residual balancing
%   3) modal tracking via MAC
%   4) result export and figure generation
%
% Required toolboxes/packages:
%   - MATPOWER (for case118 and makeBdc)
% Optional:
%   - Optimization Toolbox (quadprog). A lightweight fallback is provided.

    clc;
    close all;

    mpc = setup_ddvpp_case118_ddvpp();

    opts = struct();
    opts.max_outer = 12;
    opts.max_admm = 200;
    opts.n_modes = 3;             % number of critical oscillatory modes to secure
    opts.time_horizon = 12.0;     % seconds
    opts.dt = 0.02;
    opts.disturbance_size = -0.18; % pu power step
    opts.n_worst_loads = 8;       % scan top load buses to find worst disturbance location
    opts.rho_expand = 2.0;
    opts.rho_shrink = 0.5;
    opts.delta_max = 0.50;
    opts.delta_min = 0.01;
    opts.verbose = true;

    outdir = fullfile(pwd, 'ddvpp_outputs');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    results = ddvpp_slp_admm_case118(mpc, opts);
    ddvpp_plot_results(results, outdir);
    save(fullfile(outdir, 'ddvpp_results.mat'), 'results');

    fprintf('\nSaved outputs to: %s\n', outdir);
end
