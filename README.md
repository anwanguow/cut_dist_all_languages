The Implementation of Cut Distance in All Languages
==============

This repository contains the demonstration of the cut distance in all programming languages, including C, C++, Java, Python, Fortran, Matlab, and Mathematica (both imperative and functional programming).

The implementation of cut distance is based on computing the cut norm via a semidefinite programming (SDP) relaxation on the Stiefel manifold. As the number of nodes increases, the computed results tend to converge more closely to their true values, due to reduced randomness and better structural resolution—reflecting the asymptotic nature of cut distance in graph limit theory. The optimization algorithm relies on a random seed, which in the code is dynamically assigned based on the current time. To ensure consistent results and deterministic output across runs, the random seed should be fixed to a constant value by users. Note that the original method is proposed by Zaiwen Wen and Wotao Yin, see https://link.springer.com/article/10.1007/s10107-012-0584-1.

Algorithm Implementation
-----------------
The implementation of the cut distance (for node-labelled graphs) can be found in [Algorithm.pdf](Algorithm.pdf) of this repo.

Please note that its content is consistent with [Algorithm.md](Algorithm.md), but since GitHub does not support rendering LaTeX formulas in markdown files, the markdown version may not display correctly.

Related Article:
-----------------
The related article "Graph theory based approach to identify phase transitions in condensed matter" is published in Physical Review B (PRB) and can be accessed at https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.054116.

The repo of this article is https://github.com/anwanguow/graph_phase_transition.

In this article, cut distance is used for constructing two kinds of pseudo-order parameters, i.e., structural order parameter $\mathcal{D}_s$ and dynamical order parameter $\mathcal{D}_d$.

This article only outlines the core idea of the algorithm and does not delve into every detail. In the implementation, we use the trick of extending $A$ and $X$ to square matrices and expressing the objective function as taking the strict trace of a squared matrix, which can be found in [Algorithm.pdf](Algorithm.pdf) in this repository.

Contact:
-----------------
An Wang: amturing@outlook.com


