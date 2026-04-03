function data = ddvpp_build_model_data(mpc)
% Build the reduced nodal-frequency model data following Wang et al.
% Generation-side buses are retained in the dynamic model, and load-side
% buses are eliminated through Kron reduction.

    if ~exist('makeBdc', 'file')
        error('MATPOWER function makeBdc not found. Please add MATPOWER to MATLAB path.');
    end

    baseMVA = mpc.baseMVA;
    bus = mpc.bus;
    branch = mpc.branch;
    gen = mpc.gen;
    gen_buses = gen(:,1);
    nbus = size(bus,1);
    ngen = length(gen_buses);

    all_buses = bus(:,1);
    load_buses = setdiff(all_buses, gen_buses, 'stable');
    nload = numel(load_buses);

    [Bbus, ~, ~, ~] = makeBdc(baseMVA, bus, branch);
    Bbus = full(Bbus);

    [~, gen_pos] = ismember(gen_buses, all_buses);
    [~, load_pos] = ismember(load_buses, all_buses);

    BGG = Bbus(gen_pos, gen_pos);
    BGL = Bbus(gen_pos, load_pos);
    BLG = Bbus(load_pos, gen_pos);
    BLL = Bbus(load_pos, load_pos);

    if rcond(BLL) < 1e-12
        warning('BLL is poorly conditioned. Using pseudo-inverse fallback.');
        BLL_inv = pinv(BLL);
    else
        BLL_inv = inv(BLL);
    end

    J = BGG - BGL * BLL_inv * BLG;
    J = 0.5 * (J + J.');

    L = BGL * BLL_inv;
    F = -BLL_inv * BLG;

    Pd = max(bus(:,3), 0) / baseMVA;
    PdL = Pd(load_pos);
    if isfield(mpc, 'mu_load')
        mu_load = mpc.mu_load(:);
        if numel(mu_load) ~= nload
            error('mpc.mu_load must have one entry per load-side bus.');
        end
    else
        if sum(PdL) > 0
            mu_total = 0.8 + 0.4 * sum(PdL);
            mu_load = mu_total * PdL / sum(PdL);
        else
            mu_load = zeros(nload,1);
        end
    end

    Dload_equiv = -L * diag(mu_load) * F;

    data = struct();
    data.baseMVA = baseMVA;
    data.bus = bus;
    data.branch = branch;
    data.gen = gen;
    data.bus_ids = all_buses(:);
    data.gen_buses = gen_buses(:);
    data.load_buses = load_buses(:);
    data.ngen = ngen;
    data.nload = nload;
    data.nbus = nbus;
    data.gen_pos = gen_pos(:);
    data.load_pos = load_pos(:);
    data.Bbus = Bbus;
    data.BGG = BGG;
    data.BGL = BGL;
    data.BLG = BLG;
    data.BLL = BLL;
    data.J = J;
    data.L = L;
    data.F = F;
    data.mu_load = mu_load(:);
    data.Dload_equiv = Dload_equiv;
    data.f_base = mpc.security.f_base;
    data.omega0 = 2 * pi * mpc.security.f_base;
    data.security = mpc.security;
    data.ddvpp_bounds = mpc.ddvpp_bounds;
    data.dyn = mpc.dyn;
    data.ibr_mask = (mpc.dyn(:,2) == 2);
    data.sg_mask = (mpc.dyn(:,2) == 1);

    ibr_buses = mpc.ddvpp_bounds(:,1);
    [lia, loc] = ismember(ibr_buses, gen_buses);
    if ~all(lia)
        error('Some DDVPP buses are not generator buses.');
    end
    data.ibr_idx = loc(:);
    data.nibr = numel(loc);
    data.ibr_buses = ibr_buses(:);

    data.cost_cm = mpc.ddvpp_bounds(:,4);
    data.cost_cd = mpc.ddvpp_bounds(:,5);
    data.m_max = mpc.ddvpp_bounds(:,2);
    data.d_max = mpc.ddvpp_bounds(:,3);
end
