# 03 — 2D & 3D Centred Benchmarks, Wang (2023) 2D Ex. 4.1, Flory–Huggins & Timing

Reproduces **Section 6**, **Tables 2–3**, and **Figures 6–9** of the manuscript.

| Script | Description |
| :--- | :--- |
| `test_basic.m` | Verification of 1D/2D/3D algebraic equivalence ($3\times10^{-16}\text{--}1.3\times10^{-15}$) and bitwise invariance |
| `test_fftsym.m` | **Table 2**: FFT-length equivariance fragility of standard MATLAB FFTW ($N=32,\dots,256$) |
| `exp2d3d.m`, `fig_2d.m` | 2D X-junction ($\varepsilon=0.05$, $N=256$, ETDRK2, **Figure 6**) |
| `exp_wang41_2d.m`, `exp_wang41_fullgroup.m`, `fig_wang41_2d.m` | **Figure 7**: Head-to-head 2D comparison on Wang (2023) Ex. 4.1 ($\kappa=0.1$, $u_0=\sin(2x)\sin(2y)$, $N_p=256$, identical ESP parameters $S=1, T_e=1, N_s=10, N_{\mathrm{Iter}}=20$) extended from $t=20$ to $t=120$ comparing Wang standard ESP, $A_{\mathrm R}$ isotypic ESP, and full 16-channel $A_{\mathrm R}\times\mathcal T_2$ isotypic ESP |
| `exp2d3d_b.m`, `fig_C2d_slim.m`, `fig_3d.m` | 2D $C_2$ point-symmetric pattern ($N=192$, **Figure 8**) and 3D octant saddle ($N=96$, **Figure 9**) |
| `exp_log.m` | Non-polynomial Flory–Huggins logarithmic potential in 1D & 2D via Walsh–Hadamard orbit evaluation |
| `exp_cost.m` | **Table 3**: Wall-clock timing per ETDRK4 step & bitwise identity of full vs. reduced modes |
