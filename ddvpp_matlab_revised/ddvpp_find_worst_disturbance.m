function [worst_bus, info] = ddvpp_find_worst_disturbance(data, x, opts)
% Scan a few high-load buses and pick the disturbance location causing the
% worst full-network nodal nadir.

    Pd = max(data.bus(:,3), 0);
    [~, idx] = sort(Pd, 'descend');
    candidates = data.bus(idx(1:min(opts.n_worst_loads, numel(idx))), 1);

    model = ddvpp_linear_model(data, x);

    best_score = -inf;
    worst_bus = candidates(1);
    info = struct('candidates', candidates, 'scores', zeros(numel(candidates),1));

    for i = 1:numel(candidates)
        b = candidates(i);
        dist = ddvpp_make_disturbance(data, b, opts.disturbance_size);
        sim = ddvpp_simulate_frequency(data, model, dist, min(opts.time_horizon, 8), opts.dt, data.f_base);
        score = max(sim.nadir_hz) + 0.25 * std(sim.nadir_hz);
        info.scores(i) = score;
        if score > best_score
            best_score = score;
            worst_bus = b;
        end
    end
end
