function x = ddvpp_pack_x_current(model, data)
    m = model.m;
    d = model.d;
    x = [m(data.ibr_idx); d(data.ibr_idx)];
end
