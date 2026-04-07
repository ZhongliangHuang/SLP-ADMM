function mpc = modified_case118()
%
%   Design principles in this version
%   ---------------------------------
%   1) The electrical network stored in mpc.bus / mpc.gen / mpc.branch
%      remains valid MATPOWER data.
%   2) The nodal-frequency model still uses the reduced linear state-space
%      formulation built outside this file.
%   3) Cluster multipliers are kept only as INITIAL PRIORS used to build
%      node-wise initial values. Later optimization should act directly on
%      node-wise M and D, rather than through cluster multipliers.
%   4) Devices are distinguished as SG or IBR. For optimization purposes,
%      all IBRs are treated as GFM-IBR and use the same parameter type.
%
%   Returned fields
%   ---------------
%   Standard MATPOWER fields:
%      mpc.bus, mpc.gen, mpc.branch
%   DDVPP metadata:
%      mpc.userdata.dynamic           : compatible reconstruction metadata
%      mpc.userdata.generator_map     : reproducible host-to-template mapping
%      mpc.userdata.modified_branch   : branch scaling traceability table
%      mpc.userdata.ddvpp             : DDVPP-ready data bundle
%         .gen_dynamic_table          : per-node dynamic prior and bounds
%         .disturbance_set            : candidate disturbances
%         .security_limits            : frequency-security settings
%         .optimization_notes         : modeling notes for later scripts

mpc = case118();

%% Basic dynamic / reconstruction metadata kept for compatibility
% ----------------------------------------------------------------
mpc.userdata = struct();
mpc.userdata.dynamic = struct();
mpc.userdata.dynamic.name = 'ddvpp_ieee118_prior_case';
mpc.userdata.dynamic.base_frequency_hz = 50.0;
mpc.userdata.dynamic.line_x_scale = 0.4;
mpc.userdata.dynamic.disturbance_load_bus = 50;
mpc.userdata.dynamic.disturbance_mw = 1200.0;
mpc.userdata.dynamic.duration_s = 20.0;
mpc.userdata.dynamic.n_points = 2001;
mpc.userdata.dynamic.governor_gamma = 0.20;
mpc.userdata.dynamic.governor_tau_s = 0.75;
mpc.userdata.dynamic.governor_gain_scale = 0.90;
mpc.userdata.dynamic.droop_R = 0.05;
mpc.userdata.dynamic.inertia_scale = 1.05;
mpc.userdata.dynamic.damping_pu_fraction = 0.02;
mpc.userdata.dynamic.load_damping_fraction = 0.01;
mpc.userdata.dynamic.min_load_for_mu_mw = 10.0;
mpc.userdata.dynamic.north_corridor_scale = 0.55;
mpc.userdata.dynamic.north_edges = [ ...
    8 9;
    8 5;
    8 30;
    9 10;
    4 11;
    5 11;
    11 12;
    11 13;
    30 17;
    26 30;
    30 38];

%% Initial-prior clusters
% ------------------------------------------------------------------------
% They are used here to build node-wise INITIAL values for M, D, K, Tau.
clusters = struct([]);
clusters(1).hosts = [8 10 12];
clusters(1).transformer_x_mult = 1.25;
clusters(1).inertia_mult = 1.60;
clusters(1).damping_mult = 0.90;
clusters(1).governor_gain_mult = 0.42;
clusters(1).tau_mult = 2.00;

clusters(2).hosts = 46;
clusters(2).transformer_x_mult = 0.95;
clusters(2).inertia_mult = 0.85;
clusters(2).damping_mult = 1.00;
clusters(2).governor_gain_mult = 1.40;
clusters(2).tau_mult = 0.75;

clusters(3).hosts = 49;
clusters(3).transformer_x_mult = 0.72;
clusters(3).inertia_mult = 0.24;
clusters(3).damping_mult = 1.10;
clusters(3).governor_gain_mult = 3.00;
clusters(3).tau_mult = 0.40;

clusters(4).hosts = [54 55 56];
clusters(4).transformer_x_mult = 0.85;
clusters(4).inertia_mult = 2.80;
clusters(4).damping_mult = 1.20;
clusters(4).governor_gain_mult = 1.90;
clusters(4).tau_mult = 0.70;

clusters(5).hosts = [31 32 34 36 40 42];
clusters(5).transformer_x_mult = 1.00;
clusters(5).inertia_mult = 1.00;
clusters(5).damping_mult = 1.00;
clusters(5).governor_gain_mult = 1.15;
clusters(5).tau_mult = 0.90;

mpc.userdata.dynamic.clusters = clusters;

%% Generator fleet templates without deadband / delay
% -----------------------------------------------------------
generator_fleet = struct( ...
    'name',        {'TB6','TC6','TB10','TC10','TG3'}, ...
    'count',       {27, 9, 5, 2, 11}, ...
    'capacity_mw', {667.0, 667.0, 1050.0, 1050.0, 300.0}, ...
    'inertia_s',   {6.0, 8.0, 10.0, 10.6, 9.0}, ...
    'xdp',         {0.041, 0.039, 0.024, 0.021, 0.065} ...
);
mpc.userdata.dynamic.generator_fleet = generator_fleet;

%% Reproducible generator-to-template mapping
% -----------------------------------------------------------
host_buses = mpc.gen(:,1);
[~, order] = sort(mpc.gen(:,9), 'descend');

rows = [];
row_type = cell(0,1);
start_idx = 1;
for kk = 1:numel(generator_fleet)
    templ = generator_fleet(kk);
    subset = order(start_idx : start_idx + templ.count - 1);
    for ii = 1:numel(subset)
        idx = subset(ii);
        rows = [rows; ...
            idx, ...                                  % gen row index
            host_buses(idx), ...                      % host bus
            templ.capacity_mw, ...                    % assigned unit capacity MW
            templ.capacity_mw / mpc.baseMVA, ...      % assigned unit capacity p.u.
            templ.inertia_s, ...                      % H in seconds
            templ.xdp];                               % x'd used as coupling prior
        row_type{end+1,1} = templ.name; %#ok<AGROW>
    end
    start_idx = start_idx + templ.count;
end

generator_map = array2table(sortrows(rows,1), ...
    'VariableNames', {'gen_index','host_bus','capacity_mw','capacity_pu','inertia_s','xdp'});
generator_map.type = row_type;
generator_map = movevars(generator_map, 'type', 'After', 'host_bus');
mpc.userdata.generator_map = generator_map;

%% Device classification for DDVPP
% ------------------------------------------------------------------------
% Practical prior used here:
%   buses highlighted by the reconstruction hotspot / converter-support prior
%   are treated as IBR hosts; the remaining generator hosts are treated as SG.
%   This is only an INITIAL classification prior and can be revised later if
%   the paper framework or data source provides a stronger bus-by-bus roster.

ibr_host_buses = unique([8 10 12 46 49 54 55 56 31 32 34 36 40 42]).';
is_ibr = ismember(generator_map.host_bus, ibr_host_buses);
device_class = repmat({'SG'}, height(generator_map), 1);
device_class(is_ibr) = {'IBR'};
is_controllable = is_ibr;   % optimize IBRs only in the current DDVPP setup

%% Build node-wise initial dynamic priors
% ------------------------------------------------------------------------
capacity_pu = generator_map.capacity_pu;
base_H = generator_map.inertia_s;

M0 = mpc.userdata.dynamic.inertia_scale * 2.0 .* base_H .* capacity_pu;
D0 = mpc.userdata.dynamic.damping_pu_fraction .* capacity_pu;
K0 = mpc.userdata.dynamic.governor_gain_scale .* capacity_pu ./ mpc.userdata.dynamic.droop_R;
Tau0 = mpc.userdata.dynamic.governor_tau_s * ones(size(capacity_pu));
Gamma0 = mpc.userdata.dynamic.governor_gamma * ones(size(capacity_pu));

transformer_x0 = generator_map.xdp;

for ii = 1:numel(host_buses)
    host = host_buses(ii);
    rule = local_get_rule(host, clusters);
    if isempty(rule)
        continue;
    end
    M0(ii) = M0(ii) * rule.inertia_mult;
    D0(ii) = D0(ii) * rule.damping_mult;
    K0(ii) = K0(ii) * rule.governor_gain_mult;
    Tau0(ii) = Tau0(ii) * rule.tau_mult;
    transformer_x0(ii) = transformer_x0(ii) * rule.transformer_x_mult;
end

mu0 = mpc.userdata.dynamic.load_damping_fraction .* ...
    (max(mpc.bus(:,3), mpc.userdata.dynamic.min_load_for_mu_mw) ./ mpc.baseMVA);

%% Local feasible bounds for DDVPP optimization
% ------------------------------------------------------------------------
% Design choice:
%   - SGs are treated as fixed-support units in the current DDVPP layer.
%   - IBRs are controllable GFM-IBR nodes.
%   - Bounds are static engineering priors intended to support the local
%     feasible set required by the SLP-ADMM framework.
%
% IBR bounds are intentionally moderate:
%   M can move within [0.4, 2.5] of initial prior
%   D can move within [0.5, 3.0] of initial prior
%   K can move within [0.5, 2.0] of initial prior
%   Tau can move within [0.4, 2.0] of initial prior
%   Gamma remains within [0, 1]
%
% SG bounds are fixed at their initial value so they do not enter the current
% DDVPP optimization vector. They can be relaxed later if the paper scope
% expands to joint SG + IBR optimization.

n_gen = height(generator_map);
m_min = M0;
m_max = M0;
d_min = D0;
d_max = D0;
k_min = K0;
k_max = K0;
tau_min = Tau0;
tau_max = Tau0;
gamma_min = Gamma0;
gamma_max = Gamma0;

for ii = 1:n_gen
    if is_controllable(ii)
        m_min(ii) = 0.40 * M0(ii);
        m_max(ii) = 2.50 * M0(ii);

        d_min(ii) = max(0.50 * D0(ii), 1e-4);
        d_max(ii) = 3.00 * D0(ii);

        k_min(ii) = 0.50 * K0(ii);
        k_max(ii) = 2.00 * K0(ii);

        tau_min(ii) = 0.40 * Tau0(ii);
        tau_max(ii) = 2.00 * Tau0(ii);

        gamma_min(ii) = 0.00;
        gamma_max(ii) = 1.00;
    end
end

% Simple convex local cost prior:
% SGs fixed -> zero cost weights because not optimized here
% IBRs controllable -> positive cost weights
local_cost_quad = zeros(n_gen,1);
local_cost_lin = zeros(n_gen,1);
for ii = 1:n_gen
    if is_controllable(ii)
        % Larger flexible headroom at stronger IBR nodes can be modeled by
        % slightly smaller quadratic penalty. Here use a simple capacity-based
        % prior to keep the file self-contained.
        local_cost_quad(ii) = 1.0 / max(generator_map.capacity_mw(ii), 1.0);
        local_cost_lin(ii) = 0.01 * generator_map.capacity_pu(ii);
    end
end

%% DDVPP data bundle
% ------------------------------------------------------------------------
ddvpp = struct();

ddvpp.gen_dynamic_table = table( ...
    generator_map.gen_index, ...
    generator_map.host_bus, ...
    generator_map.type, ...
    device_class, ...
    is_controllable, ...
    generator_map.capacity_mw, ...
    generator_map.capacity_pu, ...
    generator_map.inertia_s, ...
    generator_map.xdp, ...
    transformer_x0, ...
    M0, D0, K0, Tau0, Gamma0, ...
    m_min, m_max, ...
    d_min, d_max, ...
    k_min, k_max, ...
    tau_min, tau_max, ...
    gamma_min, gamma_max, ...
    local_cost_quad, local_cost_lin, ...
    'VariableNames', { ...
        'gen_index','host_bus','template_type','device_class','is_controllable', ...
        'rated_mw','rated_pu','base_H_s','xdp_nominal','xdp_effective_initial', ...
        'M0','D0','K0','Tau0','Gamma0', ...
        'm_min','m_max','d_min','d_max','k_min','k_max','tau_min','tau_max','gamma_min','gamma_max', ...
        'local_cost_quad','local_cost_lin'} ...
    );

% Disturbance set
disturbance_id = (1:4).';
bus = [50; 49; 37; 8];
side = {'load'; 'load'; 'load'; 'load'};
deltaP_mw = [1200; 1000; 800; 800];
deltaP_pu = deltaP_mw / mpc.baseMVA;
weight = [1.00; 0.90; 0.70; 0.60];
enabled = true(size(disturbance_id));

ddvpp.disturbance_set = table(disturbance_id, bus, side, deltaP_mw, deltaP_pu, weight, enabled);

% Security limits
ddvpp.security_limits = struct();
ddvpp.security_limits.rocof_limit_hz_per_s = 1.0;
ddvpp.security_limits.nadir_limit_hz = 0.5;
ddvpp.security_limits.qss_limit_hz = 0.2;
ddvpp.security_limits.max_iterations_slp = 20;
ddvpp.security_limits.max_iterations_admm = 200;
ddvpp.security_limits.trust_region_init = 0.10;
ddvpp.security_limits.trust_region_min = 0.01;
ddvpp.security_limits.trust_region_max = 0.50;
ddvpp.security_limits.rho_init = 1.0;
ddvpp.security_limits.primal_tol = 1e-4;
ddvpp.security_limits.dual_tol = 1e-4;

ddvpp.optimization_notes = { ...
    'Clusters are initial priors only. Optimization should act on node-wise M and D directly.'; ...
    'All IBRs are treated as GFM-IBR in the current DDVPP layer.'; ...
    'SG units are fixed-support units in the current optimization setup.'; ...
    'deadband and delay are excluded from the reduced linear optimization model.'; ...
    'Local bounds are engineering priors and can be tightened after validation.' ...
    };

mpc.userdata.ddvpp = ddvpp;

%% Keep dynamic priors accessible for backward compatibility
% ------------------------------------------------------------------------
mpc.userdata.dynamic.M = M0;
mpc.userdata.dynamic.D = D0;
mpc.userdata.dynamic.K = K0;
mpc.userdata.dynamic.Tau = Tau0;
mpc.userdata.dynamic.Gamma = Gamma0;
mpc.userdata.dynamic.mu = mu0;

%% Branch modification summary for traceability
% ------------------------------------------------------------------------
branch_index = (1:size(mpc.branch,1)).';
from_bus = mpc.branch(:,1);
to_bus = mpc.branch(:,2);
x_original = mpc.branch(:,4);
x_modified_nominal = x_original * mpc.userdata.dynamic.line_x_scale;
north_corridor_edge = false(size(branch_index));
b_effective_used = zeros(size(branch_index));
x_effective_if_using_B = x_modified_nominal;

for rr = 1:size(mpc.branch,1)
    pair = sort([from_bus(rr), to_bus(rr)]);
    is_weak = any(all(mpc.userdata.dynamic.north_edges == repmat(pair, size(mpc.userdata.dynamic.north_edges,1),1), 2));
    north_corridor_edge(rr) = is_weak;
    if is_weak
        x_effective_if_using_B(rr) = x_modified_nominal(rr) / mpc.userdata.dynamic.north_corridor_scale;
    end
    tap = mpc.branch(rr,9);
    if tap == 0
        tap = 1.0;
    end
    b_effective_used(rr) = 1.0 / (x_effective_if_using_B(rr) * tap);
end

mpc.userdata.modified_branch = table(branch_index, from_bus, to_bus, x_original, ...
    x_modified_nominal, north_corridor_edge, x_effective_if_using_B, b_effective_used);

end

function rule = local_get_rule(host, clusters)
    rule = [];
    for kk = 1:numel(clusters)
        if any(clusters(kk).hosts == host)
            rule = clusters(kk);
            return;
        end
    end
end
