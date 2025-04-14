function cut_dist
    n = 400;
    p1 = 0.2;
    p2 = 0.5;
    A = triu(rand(n) < p1, 1);
    A = A + A';
    B = triu(rand(n) < p2, 1);
    B = B + B';
    s = cut_distance(A, B);
    fprintf('The cut distance between A and B: %f\n', s);
end

function s = cut_distance(G, H)
    n = size(G, 1);
    A = (G - H) / (n * n);
    s = cut_norm(A);
end

function sdp_result = cut_norm(A)
    n1 = size(A, 1);
    A_col_sum = sum(A, 1);
    A_row_sum = sum(A, 2);
    A_tot = sum(A_col_sum);
    A = [A, -1.0 * A_row_sum];
    A = [A; [-1.0 * A_col_sum, A_tot]];
    p = max(min(round(sqrt(2 * n1) / 2), 100), 1);
    n2 = 2 * n1 + 2;
    x0 = randn(p, n2);
    nrmx0 = sum(x0.^2, 1);
    x0 = x0 ./ sqrt(nrmx0);
    [x, ~] = optimize(x0, @cut_norm_quad, A);
    U = x(:, 1:(n2/2));
    V = x(:, (n2/2+1):end);
    sdp_result = abs(sum(sum(A .* (U' * V)))) / 4.0;
end

function [f, g] = cut_norm_quad(V, A)
    n = size(A, 1);
    Us = V(:, 1:n);
    Vs = V(:, (n+1):end);
    g = 2 * [Vs * A', Us * A];
    f = (sum(sum(g(:, 1:n) .* Us)) + sum(sum(g(:, n+1:end) .* Vs))) / 2;
end

function [x, g] = optimize(x, fun, args, xtol, ftol, gtol, rho, eta, gamma, tau, nt, mxitr)
    if nargin < 4, xtol = 1e-8; end
    if nargin < 5, ftol = 1e-10; end
    if nargin < 6, gtol = 1e-8; end
    if nargin < 7, rho = 1e-4; end
    if nargin < 8, eta = 0.1; end
    if nargin < 9, gamma = 0.85; end
    if nargin < 10, tau = 1e-3; end
    if nargin < 11, nt = 5; end
    if nargin < 12, mxitr = 600; end
    crit = ones(mxitr, 3);
    [n, ~] = size(x);
    nrmx = sum(x.^2, 1);
    if norm(nrmx) > 1e-8
        x = x ./ sqrt(nrmx);
    end
    [f, g] = fun(x, args);
    xtg = sum(x .* g, 1);
    gg = sum(g.^2, 1);
    xx = sum(x.^2, 1);
    xxgg = xx .* gg;
    dtX = x .* xtg - g;
    nrmG = norm(dtX, 'fro');
    Q = 1;
    Cval = f;
    tau_orig = tau;
    tau = tau_orig;
    for itr = 1:mxitr
        xp = x;
        fp = f;
        gp = g;
        dtXP = dtX;
        nls = 1;
        deriv = rho * nrmG^2;
        while true
            tau2 = tau / 2;
            beta = 1 + tau2^2 * (-xtg.^2 + xxgg);
            a1 = ((1 + tau2 * xtg).^2 - tau2^2 * xxgg) ./ beta;
            a2 = -tau * xx ./ beta;
            x = xp .* a1 + gp .* a2;
            [f, g] = fun(x, args);
            if f <= Cval - tau * deriv || nls >= 5
                break;
            end
            tau = eta * tau;
            nls = nls + 1;
        end
        xtg = sum(x .* g, 1);
        gg = sum(g.^2, 1);
        xx = sum(x.^2, 1);
        xxgg = xx .* gg;
        dtX = x .* xtg - g;
        nrmG = norm(dtX, 'fro');
        s = x - xp;
        XDiff = norm(s, 'fro') / sqrt(n);
        FDiff = abs(fp - f) / (abs(fp) + 1);
        crit(itr, :) = [nrmG, XDiff, FDiff];
        start_index = max(1, itr - nt + 1);
        mcrit = mean(crit(start_index:itr, :), 1);
        if (XDiff < xtol && FDiff < ftol) || nrmG < gtol || all(mcrit(2:3) < 10 * [xtol, ftol])
            break;
        end
        y = dtX - dtXP;
        sy = abs(sum(sum(s .* y)));
        tau = tau_orig;
        if sy > 0
            if mod(itr, 2) == 0
                tau = sum(sum(s.^2)) / sy;
            else
                tau = sy / sum(sum(y.^2));
            end
            tau = max(min(tau, 1e20), 1e-20);
        end
        Qp = Q;
        Q = gamma * Qp + 1;
        Cval = (gamma * Qp * Cval + f) / Q;
    end
end
