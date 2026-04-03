function out = ddvpp_run_admm(data, x, crit, req, delta, rho, mpc, opts)
% Inner ADMM loop for the separable quadratic allocation subproblem.
% Constraints:
%   -S_sigma * dx >= mode_req
%   sum_i Delta d_i >= qss_req
%   Delta m_i >= roc_deficit_ibr(i)
% with local trust-region and box bounds.

    nibr = data.nibr;
    p = 2 * nibr;
    K = numel(req.mode_req);

    if req.uncontrollable_roc && data.security.hard_abort_on_roc_infeasible
        out = make_infeasible_output(p, opts.max_admm, rho, 'RoCoF lower bound is violated on non-IBR buses and cannot be corrected by DDVPP variables.');
        return;
    end

    z = zeros(K+1,1);
    y = zeros(K+1,1);
    z_lb = [req.mode_req; req.qss_req];

    S = [-crit.S_sigma; [zeros(1,nibr), ones(1,nibr)]];

    dx_all = zeros(p, nibr);
    cost = 0;
    converged = false;
    infeasible = false;
    infeasible_msg = '';
    primal_hist = zeros(opts.max_admm,1);
    dual_hist = zeros(opts.max_admm,1);
    rho_hist = zeros(opts.max_admm,1);

    [m0, d0] = ddvpp_unpack_x(data, x);
    m0_ibr = m0(data.ibr_idx);
    d0_ibr = d0(data.ibr_idx);

    for t = 1:opts.max_admm
        z_old = z;

        for i = 1:nibr
            Si = S(:, [i, nibr+i]);
            ci = [data.cost_cm(i), data.cost_cd(i)];
            lb = [max(-delta, 1e-6 - m0_ibr(i)); max(-delta, 1e-6 - d0_ibr(i))];
            ub = [min(delta, data.m_max(i) - m0_ibr(i)); min(delta, data.d_max(i) - d0_ibr(i))];
            lb(1) = max(lb(1), req.roc_deficit_ibr(i));

            if ub(1) + 1e-12 < lb(1)
                infeasible = true;
                infeasible_msg = sprintf('RoCoF lower bound infeasible at IBR bus %d.', data.ibr_buses(i));
                break;
            end
            if ub(2) + 1e-12 < lb(2)
                infeasible = true;
                infeasible_msg = sprintf('Damping trust region infeasible at IBR bus %d.', data.ibr_buses(i));
                break;
            end

            g_others = S * sum(dx_all,2) - Si * dx_all([i, nibr+i], i);
            q = g_others - z;
            dxi = local_qp_solver(ci, Si, q, y, rho, lb, ub);
            dx_all([i, nibr+i], i) = dxi;
        end

        if infeasible
            break;
        end

        g = zeros(K+1,1);
        for i = 1:nibr
            g = g + S(:, [i, nibr+i]) * dx_all([i, nibr+i], i);
        end

        z = max(g + y / rho, z_lb);
        y = y + rho * (g - z);

        r = g - z;
        s = rho * (z - z_old);
        nr = norm(r);
        ns = norm(s);
        primal_hist(t) = nr;
        dual_hist(t) = ns;
        rho_hist(t) = rho;

        if nr > mpc.admm.mu_balance * ns
            rho = rho * mpc.admm.tau_inc;
        elseif ns > mpc.admm.mu_balance * nr
            rho = rho / mpc.admm.tau_dec;
        end

        if nr <= mpc.admm.tol_primal && ns <= mpc.admm.tol_dual
            converged = true;
            break;
        end
    end

    if infeasible
        out = make_infeasible_output(p, max(1, t), rho, infeasible_msg, primal_hist, dual_hist, rho_hist);
        return;
    end

    dx_global = zeros(p,1);
    for i = 1:nibr
        dx_global([i, nibr+i]) = dx_global([i, nibr+i]) + dx_all([i, nibr+i], i);
        cost = cost + data.cost_cm(i) * dx_all(i,i)^2 + data.cost_cd(i) * dx_all(nibr+i,i)^2;
    end

    out = struct();
    out.dx_full = dx_global(:)';
    out.dx_global = dx_global;
    out.cost = cost;
    out.iters = t;
    out.converged = converged;
    out.infeasible = false;
    out.message = '';
    out.rho_final = rho;
    out.primal_residual = primal_hist(1:t);
    out.dual_residual = dual_hist(1:t);
    out.rho_hist = rho_hist(1:t);
end

function x = local_qp_solver(ci, Si, q, y, rho, lb, ub)
% Minimize c1*x1^2 + c2*x2^2 + y'*(Si*x + q) + rho/2 * ||Si*x + q||^2
% subject to box bounds.

    H = 2 * diag(ci) + rho * (Si' * Si);
    f = Si' * y + rho * (Si' * q);

    if exist('quadprog', 'file') == 2
        qp_opts = optimoptions('quadprog', 'Display', 'off');
        x = quadprog(H, f, [], [], [], [], lb, ub, [], qp_opts);
        if isempty(x)
            x = min(max(-H \ f, lb), ub);
        end
    else
        x = min(max(zeros(2,1), lb), ub);
        alpha = 1 / max(real(eig(H)));
        for k = 1:200
            g = H * x + f;
            x = x - alpha * g;
            x = min(max(x, lb), ub);
        end
    end
end

function out = make_infeasible_output(p, t, rho, msg, primal_hist, dual_hist, rho_hist)
    if nargin < 5
        primal_hist = zeros(t,1);
        dual_hist = zeros(t,1);
        rho_hist = rho * ones(t,1);
    else
        primal_hist = primal_hist(1:t);
        dual_hist = dual_hist(1:t);
        rho_hist = rho_hist(1:t);
    end
    out = struct();
    out.dx_full = zeros(1,p);
    out.dx_global = zeros(p,1);
    out.cost = inf;
    out.iters = t;
    out.converged = false;
    out.infeasible = true;
    out.message = msg;
    out.rho_final = rho;
    out.primal_residual = primal_hist;
    out.dual_residual = dual_hist;
    out.rho_hist = rho_hist;
end
