# Paper Numerical Experiments (`experiments/`)

All scripts below reproduce the figures and tables of the manuscript. Before running any experiment script in MATLAB, add the core library to your path:

```matlab
addpath(genpath('../core'));
```

(Each script writes intermediate `.mat` data to `../results/` and figures to `../figs/`; create those directories if needed, or run from the repository root.)

## Subfolder Guide

| Subfolder | Paper Section / Figures / Tables | Scripts |
| :--- | :--- | :--- |
| [`01_1d_wang51_and_high_order/`](01_1d_wang51_and_high_order/) | Sec. 5.2–5.3, **Table 1**, **Figs. 2–3**: Wang (2023) Ex. 5.1 breakdown law ($S=0,1$, ESP $T_e=20$ refutation) & time-evolving departure of 7 high-order integrators (BDF2-SAV, ARS(4,4,3), SBDF3, SBDF4, IFRK4, SDC4, Strang2) up to $T=1000$ | `explore1d.m`, `exp_law.m`, `exp_highorder.m`, `exp_highorder_evol.m`, `add_snap47.m`, `fig_1d_wang.m`, `fig_highorder_1d.m` |
| [`02_1d_faithfulness_and_cn/`](02_1d_faithfulness_and_cn/) | Sec. 5.4–5.5, **Figs. 4–5**: Sub-roundoff faithfulness ($\delta\in[10^{-200},10^{-2}]$) vs. 120-digit VPA reference; companion adaptive-cG Crank–Nicolson benchmark ($N=4096$, $\eps=1$) & multi-kink drift | `exp_delta.m`, `ref_vpa.m`, `fig_faithful.m`, `exp_draft57.m`, `fig_draft57.m` |
| [`03_2d_3d_centred_and_wang41/`](03_2d_3d_centred_and_wang41/) | Sec. 6.1–6.7, **Tables 2–3**, **Figs. 6–9**: 1D/2D/3D equivalence, FFT-$N$ equivariance table, 2D X-junction, **head-to-head 2D Wang (2023) Ex. 4.1 under identical ESP parameters ($A_{\mathrm R}$ vs. $A_{\mathrm R}\times\mathcal T_2$)**, 2D $C_2$ point-symmetric pattern, 3D octant saddle, Flory–Huggins log potential, wall-clock cost | `test_basic.m`, `test_fftsym.m`, `exp2d3d.m`, `exp2d3d_b.m`, `exp_wang41_2d.m`, `exp_wang41_fullgroup.m`, `exp_log.m`, `exp_cost.m`, `fig_2d.m`, `fig_wang41_2d.m`, `fig_C2d_slim.m`, `fig_3d.m` |
| [`04_shifted_t2_and_cfism/`](04_shifted_t2_and_cfism/) | Sec. 7.1–7.4, **Table 4**, **Figs. 10–13**: Half-grid Type II & $\mathcal T_2$ verification, Wang (2023) Ex. 5.2 refutation & 1D shifted evolution ($T=1000$), CF-ISM horizon law $T_\ast(d_0)$, irrationally shifted 2D mixed-harmonic X-junction (CF-ISM) & 3D octant saddle ($\mathcal T_2$) | `test_shifted.m`, `exp_shifted_1d.m`, `exp_wang52_evol.m`, `exp_cfism_2d.m`, `exp_shifted_2d3d.m`, `fig_wang52_evol.m`, `fig_cfism_horizon.m`, `fig_cfism_2d.m`, `fig_shifted_2d3d.m` |
| [`05_pde_extensions/`](05_pde_extensions/) | Sec. 7.5: Drop-in extensibility to first-order SAV, 2D conservative Allen–Cahn (CAC), and 2D fourth-order Cahn–Hilliard | `exp_extensibility.m` |
