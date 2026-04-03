function mpc = setup_ddvpp_case118_ddvpp()
% SETUP_DDVPP_CASE118_DDVPP
% Revised deterministic IEEE 118-bus DDVPP initialization.
%
% This version only modifies DDVPP-related settings and leaves MATPOWER
% case118 untouched. Compared with the previous setup draft, this file:
%   1) spreads IBRs over a wider electrical area,
%   2) avoids large groups of exactly identical dynamic parameters,
%   3) lifts hotspot inertia capacity while keeping QSS moderate,
%   4) preserves compatibility with the existing codebase.

    mpc = case118;

    %% basic indexing
    gen_buses = mpc.gen(:, 1);
    num_gen = length(gen_buses);
    all_buses = mpc.bus(:, 1);
    load_buses = setdiff(all_buses, gen_buses, 'stable');

    %% DDVPP IBR deployment
    % Use a mixed deployment: several hotspot buses near the vulnerable area,
    % plus backbone support buses with large generation capacity.
    ibr_buses = [59 61 65 66 80 89 90 92 99 100 103 107 110 112 116]';
    ibr_buses = intersect(ibr_buses(:), gen_buses(:), 'stable');

    hot_buses  = [80 89 90 100 116]';
    mid_buses  = [92 99 103 107 110 112]';
    back_buses = [59 61 65 66]';

    %% dynamic model parameters
    % dyn columns: [bus_id, type(1=SG,2=IBR), m0, d0, k, tau, gamma]
    mpc.dyn = zeros(num_gen, 7);
    for i = 1:num_gen
        bus_id = gen_buses(i);
        mpc.dyn(i, 1) = bus_id;

        % small deterministic spread to avoid repeated identical modes
        eta = 0.01 * mod(bus_id, 7);
        zeta = 0.005 * mod(bus_id, 5);

        if ismember(bus_id, hot_buses)
            % hotspot suppression layer
            mpc.dyn(i, 2) = 2;
            mpc.dyn(i, 3) = 0.34 + eta;          % m0
            mpc.dyn(i, 4) = 1.00 + 0.20*zeta;    % d0 
            mpc.dyn(i, 5) = 14.5 + 0.3*mod(bus_id,3); % k
            mpc.dyn(i, 6) = 0.20 + 0.01*zeta;    % tau
            mpc.dyn(i, 7) = 0.00;                % gamma
        elseif ismember(bus_id, mid_buses)
            % modal damping layer
            mpc.dyn(i, 2) = 2;
            mpc.dyn(i, 3) = 0.28 + eta;
            mpc.dyn(i, 4) = 0.88 + 0.20*zeta;
            mpc.dyn(i, 5) = 13.2 + 0.25*mod(bus_id,4);
            mpc.dyn(i, 6) = 0.24 + 0.01*zeta;
            mpc.dyn(i, 7) = 0.00;
        elseif ismember(bus_id, back_buses)
            % backbone support layer
            mpc.dyn(i, 2) = 2;
            mpc.dyn(i, 3) = 0.30 + eta;
            mpc.dyn(i, 4) = 0.95 + 0.20*zeta;
            mpc.dyn(i, 5) = 13.8 + 0.2*mod(bus_id,5);
            mpc.dyn(i, 6) = 0.22 + 0.01*zeta;
            mpc.dyn(i, 7) = 0.00;
        else
            % synchronous-generator baseline
            mpc.dyn(i, 2) = 1;
            mpc.dyn(i, 3) = 6.30 + 0.03*mod(bus_id,4);
            mpc.dyn(i, 4) = 1.85 + 0.04*mod(bus_id,3);
            mpc.dyn(i, 5) = 18.0 + 0.2*mod(bus_id,5);
            mpc.dyn(i, 6) = 0.44 + 0.01*zeta;
            mpc.dyn(i, 7) = 0.12 + 0.01*mod(bus_id,3);
        end
    end

    %% DDVPP optimization bounds and costs
    % columns: [bus_id, m_max, d_max, c_m, c_d]
    % The current solver interprets m_max and d_max as absolute upper bounds.
    num_ibr = length(ibr_buses);
    mpc.ddvpp_bounds = zeros(num_ibr, 5);
    for i = 1:num_ibr
        bus_id = ibr_buses(i);
        mpc.ddvpp_bounds(i, 1) = bus_id;

        if ismember(bus_id, hot_buses)
            mpc.ddvpp_bounds(i, 2) = 3.78;  % absolute m upper bound
            mpc.ddvpp_bounds(i, 3) = 1.65;  % absolute d upper bound
            mpc.ddvpp_bounds(i, 4) = 0.95;
            mpc.ddvpp_bounds(i, 5) = 0.48;
        elseif ismember(bus_id, mid_buses)
            mpc.ddvpp_bounds(i, 2) = 3.62;
            mpc.ddvpp_bounds(i, 3) = 1.40;
            mpc.ddvpp_bounds(i, 4) = 1.08;
            mpc.ddvpp_bounds(i, 5) = 0.52;
        else
            mpc.ddvpp_bounds(i, 2) = 3.68;
            mpc.ddvpp_bounds(i, 3) = 1.45;
            mpc.ddvpp_bounds(i, 4) = 1.00;
            mpc.ddvpp_bounds(i, 5) = 0.50;
        end
    end

    %% explicit load-damping profile on eliminated load buses
    Pd = max(mpc.bus(:, 3), 0) / mpc.baseMVA;
    [~, load_pos] = ismember(load_buses, all_buses);
    PdL = Pd(load_pos);
    mu_total = 14.0;
    if sum(PdL) > 0
        mpc.mu_load = mu_total * PdL / sum(PdL);
    else
        mpc.mu_load = zeros(numel(load_buses), 1);
    end

    %% security settings
    mpc.security = struct( ...
        'RoCoF_lim', 1.5, ...
        'Nadir_lim', 0.8, ...
        'QSS_lim',   0.2, ...
        'f_base',    50.0, ...
        'mode_margin_eps', 1e-3, ...
        'nadir_time_floor', 0.20, ...
        'mode_time_fraction', 0.50, ...
        'envelope_tol', 1e-3, ...
        'hard_abort_on_roc_infeasible', true);

    %% ADMM settings
    mpc.admm = struct( ...
        'rho_init', 1.0, ...
        'delta_init', 0.04, ...
        'tol_primal', 1e-4, ...
        'tol_dual', 1e-4, ...
        'mu_balance', 10.0, ...
        'tau_inc', 2.0, ...
        'tau_dec', 2.0);
end
