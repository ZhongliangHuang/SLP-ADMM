function [m_all, d_all] = ddvpp_unpack_x(data, x)
    m_all = data.dyn(:,3);
    d_all = data.dyn(:,4);
    n = data.nibr;
    m_all(data.ibr_idx) = x(1:n);
    d_all(data.ibr_idx) = x(n+1:2*n);
end
