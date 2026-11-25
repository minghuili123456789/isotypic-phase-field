# 01 — 1D Benchmark (Wang 2023 Ex. 5.1) & Seven High-Order Integrators

Reproduces **Section 5.2–5.3**, **Table 1**, **Figure 2**, and **Figure 3** of the manuscript.

| Script | Description |
| :--- | :--- |
| `explore1d.m`, `exp_law.m` | Breakdown-time law $T_\ast = \lambda_\perp^{-1}\log(1/\mathbf{u})+O(1)$ across IMEX1, ETDRK2, ETDRK4, CN, and ESP ($S=0$ and $S=1$, including ESP $T_e=20$ refutation) |
| `fig_1d_wang.m` | Generates `fig_1d_wang.pdf` (**Figure 2**) |
| `exp_highorder.m` | Convergence-order, non-symmetric equivalence, and symmetry tests for 7 high-order integrators: BDF2-SAV, ARS(4,4,3), SBDF3, SBDF4, IFRK4, SDC4, Strang2 |
| `exp_highorder_evol.m`, `add_snap47.m` | Long-time trajectories ($0 \le t \le 1000$) and mid-collapse snapshot at $t=47$ |
| `fig_highorder_1d.m` | Generates `fig_highorder_1d.pdf` (**Figure 3**) |

Run any script directly in MATLAB; `../../core` is added automatically.
