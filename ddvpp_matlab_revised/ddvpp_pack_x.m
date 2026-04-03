function x = ddvpp_pack_x(data)
    m = data.dyn(:,3);
    d = data.dyn(:,4);
    x = [m(data.ibr_idx); d(data.ibr_idx)];
end
