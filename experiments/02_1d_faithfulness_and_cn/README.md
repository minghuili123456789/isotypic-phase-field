# 02 — 1D Faithfulness (120-Digit VPA Reference) & Companion Crank–Nicolson Benchmark

Reproduces **Section 5.4–5.5**, **Figure 4**, and **Figure 5** of the manuscript.

| Script | Description |
| :--- | :--- |
| `exp_delta.m` | Faithfulness test for genuine asymmetries $\delta\in[10^{-200},10^{-2}]$ comparing isotypic double, standard double/single, and the LQTY (2022) parity filter |
| `ref_vpa.m` | Independent 120-digit variable-precision arithmetic (VPA) reference (including 120-digit DFT matrix) for $\delta=10^{-10}, 10^{-30}, 10^{-60}$ |
| `fig_faithful.m` | Generates unified 1×3 `fig_faithful.pdf` (**Figure 4**) |
| `exp_draft57.m`, `fig_draft57.m` | Companion adaptive-cG Crank–Nicolson benchmark ($N=4096$, $\varepsilon=1$, $u_0=\tfrac12(\sin x+\sin 4x)$, $T=100$) comparing standard CN, zero-mode CN, and isotypic CN (**Figure 5**) |
