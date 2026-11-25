# 04 — Arbitrarily Shifted Symmetries ($\mathcal T_2$ and CF-ISM) in 1D, 2D and 3D

Reproduces **Section 7.1–7.4**, **Table 4**, and **Figures 10–13** of the manuscript.

| Script | Description |
| :--- | :--- |
| `test_shifted.m` | Verification of half-grid Type II (DCT-II/DST-II) and $\mathcal T_2$ representations in 1D/2D/3D |
| `exp_shifted_1d.m`, `fig_cfism_horizon.m` | **Table 4** (all 16 cases of Wang 2023 Ex. 5.2(a)–(d) × 4 shifts including irrational $\theta=0.37123456789$) and **Figure 11** (CF-ISM horizon law $T_\ast(d_0)$) |
| `exp_wang52_evol.m`, `fig_wang52_evol.m` | **Figure 10**: 4-method 1D shifted evolution ($T=1000$, Wang Ex. 5.2(a) $\theta=\pi/4$) visually refuting standard ETDRK4, Wang ESP ($S=1, T_e=1$), and LQTY parity filter vs. our $\mathcal T_2$ ESP & CF-ISM |
| `exp_cfism_2d.m`, `fig_cfism_2d.m` | **Figure 12**: 2D irrationally shifted mixed-harmonic X-junction (CF-ISM) |
| `exp_shifted_2d3d.m`, `fig_shifted_2d3d.m` | **Figure 13**: 3D irrationally shifted octant saddle ($\mathcal T_2$, 100% pure double precision, zero shift detection) |
