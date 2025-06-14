The Implementation of Cut Distance in All Languages
==============

This repository contains the demonstration of the cut distance in all programming languages, including C, C++, Java, Python, Fortran, Matlab, and Mathematica (both imperative and functional programming).

The related article "Graph theory based approach to identify phase transitions in condensed matter" is published in Physical Review B (PRB) and can be accessed at https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.054116.

The repo of this article is https://github.com/anwanguow/graph_phase_transition.

The implementation of cut distance is based on computing the cut norm via a semidefinite programming (SDP) relaxation on the Stiefel manifold. As the number of nodes increases, the computed results tend to converge more closely to their true values, due to reduced randomness and better structural resolution—reflecting the asymptotic nature of cut distance in graph limit theory. The optimization algorithm relies on a random seed, which in the code is dynamically assigned based on the current time. To ensure consistent results and deterministic output across runs, the random seed should be fixed to a constant value by users. Note that the original method is proposed by Zaiwen Wen and Wotao Yin, see https://link.springer.com/article/10.1007/s10107-012-0584-1.

Reference
-----------------

Please consider adding the following citation if you find it is useful to your research.

```bibtex
@article{PhysRevB.111.054116,
  title = {Graph theory based approach to identify phase transitions in condensed matter},
  author = {Wang, An and Sosso, Gabriele C.},
  journal = {Phys. Rev. B},
  volume = {111},
  issue = {5},
  pages = {054116},
  numpages = {10},
  year = {2025},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.111.054116},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.111.054116}
}
```

Contact:
-----------------
An Wang: amturing@outlook.com


