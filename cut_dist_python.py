import numpy as np
import networkx as nx

def cut_distance(G, H):
    n = len(G)
    A = (G - H) / (n * n)
    return cut_norm(A)

def cut_norm(A):
    n1 = len(A)
    A_col_sum = np.sum(A, axis=0)
    A_row_sum = np.sum(A, axis=1)
    A_tot = np.sum(A_col_sum)
    A = np.c_[A, -1.0 * A_row_sum]
    A = np.r_[A, [np.concatenate((-1.0 * A_col_sum, [A_tot]))]]
    p = int(max(min(round(np.sqrt(2 * n1) / 2), 100), 1))
    n2 = 2 * n1 + 2
    x0 = np.random.randn(p, n2)
    nrmx0 = np.sum(x0 * x0, axis=0)
    x0 = np.divide(x0, np.sqrt(nrmx0))
    x, g = optimize(x0, cut_norm_quad, A)
    U = x[:, :n2 // 2]
    V = x[:, n2 // 2:]
    sdp_result = np.abs(np.sum(A * np.matmul(U.T, V))) / 4.0
    return sdp_result

def cut_norm_quad(V, A):
    n = len(A)
    Vs = V[:, n:]
    Us = V[:, :n]
    g = 2 * np.c_[np.matmul(Vs, A.T), np.matmul(Us, A)]
    f = (np.sum(g[:, :n] * Us) + np.sum(g[:, n:] * Vs)) / 2
    return f, g

def optimize(x, fun, args, xtol=1e-8, ftol=1e-10, gtol=1e-8, rho=1e-4,
             eta=0.1, gamma=0.85, tau=1e-3, nt=5, mxitr=600):
    crit = np.ones((mxitr, 3))
    (n, p) = x.shape
    nrmx = np.sum(x * x, axis=0)
    if np.linalg.norm(nrmx) > 1e-8:
        x = np.divide(x, np.sqrt(nrmx))
    f, g = fun(x, args)
    xtg = np.sum(x * g, axis=0)
    gg = np.sum(g * g, axis=0)
    xx = np.sum(x * x, axis=0)
    xxgg = xx * gg
    dtX = x * xtg - g
    nrmG = np.linalg.norm(dtX, 'fro')
    Q = 1
    Cval = f
    tau_orig = tau
    tau = tau_orig
    for itr in range(mxitr):
        xp = x
        fp = f
        gp = g
        dtXP = dtX
        nls = 1
        deriv = rho * nrmG**2
        while True:
            tau2 = tau / 2
            beta = 1 + tau2**2 * (-xtg**2 + xxgg)
            a1 = ((1 + tau2 * xtg)**2 - tau2**2 * xxgg) / beta
            a2 = -tau * xx / beta
            x = xp * a1 + gp * a2
            f, g = fun(x, args)
            if f <= Cval - tau * deriv or nls >= 5:
                break
            tau = eta * tau
            nls = nls + 1
        xtg = np.sum(x * g, axis=0)
        gg = np.sum(g * g, axis=0)
        xx = np.sum(x * x, axis=0)
        xxgg = xx * gg
        dtX = x * xtg - g
        nrmG = np.linalg.norm(dtX, 'fro')
        s = x - xp
        XDiff = np.linalg.norm(s, 'fro') / np.sqrt(n)
        FDiff = abs(fp - f) / (abs(fp) + 1)
        crit[itr, :] = [nrmG, XDiff, FDiff]
        mcrit = np.mean(crit[itr - min(nt, itr):itr + 1, :], axis=0)
        if ((XDiff < xtol and FDiff < ftol) or nrmG < gtol
                or np.all(mcrit[1:] < 10 * np.array([xtol, ftol]))):
            break
        y = dtX - dtXP
        sy = np.sum(s * y)
        sy = abs(sy)
        tau = tau_orig
        if sy > 0:
            if np.mod(itr, 2) == 0:
                tau = np.sum(s * s) / sy
            else:
                tau = sy / np.sum(y * y)
            tau = max(min(tau, 1e20), 1e-20)
        Qp = Q
        Q = gamma * Qp + 1
        Cval = (gamma * Qp * Cval + f) / Q
    return x, g

n = 30
p_1 = 0.2
p_2 = 0.5

A = nx.erdos_renyi_graph(n, p_1);
A = np.array(nx.adjacency_matrix(A).todense());
B = nx.erdos_renyi_graph(n, p_2);
B = np.array(nx.adjacency_matrix(B).todense());

s = cut_distance(A, B);
print("The cut distance between A and B: ", s)

