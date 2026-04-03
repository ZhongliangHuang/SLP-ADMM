function dist = ddvpp_make_disturbance(data, bus_id, step_pu)
% Map a disturbance at any physical bus to the generation-side disturbance
% vector used by the reduced model.

    n = data.ngen;
    dist = zeros(n,1);

    [tfG, idxG] = ismember(bus_id, data.gen_buses);
    if tfG
        dist(idxG) = step_pu;
        return;
    end

    [tfL, idxL] = ismember(bus_id, data.load_buses);
    if tfL
        dist = step_pu * data.L(:, idxL);
        return;
    end

    error('Bus %d is not found in the system bus list.', bus_id);
end
