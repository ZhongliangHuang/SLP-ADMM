function model = ddvpp_linear_model(data, x)
% Build the Wang-style nodal-frequency state matrix.
% xdot = A x + B_in * DeltaP_G, where DeltaP_G is the disturbance vector
% referred to the generation side.

    [m, d] = ddvpp_unpack_x(data, x);
    k = data.dyn(:,5);
    tau = data.dyn(:,6);
    gamma = data.dyn(:,7);
    n = data.ngen;

    if any(m <= 0)
        error('All inertia parameters must remain positive.');
    end
    if any(tau <= 0)
        error('All time constants must remain positive.');
    end

    M = diag(m);
    Minv = diag(1 ./ m);
    D = diag(d);
    K = diag(k);
    T = diag(tau);
    Gamma = diag(gamma);
    N = K * Gamma * Minv;
    Dtilde = D - data.L * diag(data.mu_load) * data.F;

    A = zeros(3*n, 3*n);
    A(1:n, n+1:2*n) = data.omega0 * eye(n);
    A(n+1:2*n, 1:n) = -Minv * data.J;
    A(n+1:2*n, n+1:2*n) = -Minv * Dtilde;
    A(n+1:2*n, 2*n+1:3*n) = Minv;
    A(2*n+1:3*n, 1:n) = N * data.J;
    A(2*n+1:3*n, n+1:2*n) = N * Dtilde - (T \ K);
    A(2*n+1:3*n, 2*n+1:3*n) = -N - inv(T);

    B_in = zeros(3*n, n);
    B_in(n+1:2*n, :) = Minv;
    B_in(2*n+1:3*n, :) = -N;

    model = struct();
    model.A = A;
    model.B = B_in;
    model.m = m(:);
    model.d = d(:);
    model.k = k(:);
    model.tau = tau(:);
    model.gamma = gamma(:);
    model.M = M;
    model.Minv = Minv;
    model.D = D;
    model.Dtilde = Dtilde;
    model.K = K;
    model.T = T;
    model.Gamma = Gamma;
    model.N = N;
end
