function s = ddvpp_mode_sensitivity(data, model, u, v)
% Analytic real-part sensitivity of a simple eigenvalue with respect to each
% IBR virtual inertia m_i and damping d_i.
% Uses d lambda / d p = u^H (dA/dp) v with u^H v = 1.

    nibr = data.nibr;
    s = zeros(1, 2*nibr);

    Dtilde = model.Dtilde;
    J = data.J;

    uomega = u(data.ngen+1:2*data.ngen);
    ug = u(2*data.ngen+1:3*data.ngen);

    vtheta = v(1:data.ngen);
    vomega = v(data.ngen+1:2*data.ngen);
    vg = v(2*data.ngen+1:3*data.ngen);

    for ii = 1:nibr
        gi = data.ibr_idx(ii);
        mi = model.m(gi);
        kii = model.k(gi);
        gammai = model.gamma(gi);
        Nii = model.N(gi, gi);

        common = J(gi,:) * vtheta + Dtilde(gi,:) * vomega - vg(gi);

        dlam_dm = ((conj(uomega(gi)) - kii * gammai * conj(ug(gi))) / (mi^2)) * common;
        dlam_dd = (-conj(uomega(gi)) / mi + Nii * conj(ug(gi))) * vomega(gi);

        s(ii) = real(dlam_dm);
        s(nibr + ii) = real(dlam_dd);
    end
end
