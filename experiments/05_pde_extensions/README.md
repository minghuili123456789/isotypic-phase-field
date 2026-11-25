# 05 — Drop-In PDE & Scheme Extensions (SAV, Conservative Allen–Cahn, Cahn–Hilliard)

Reproduces **Section 7.5** of the manuscript.

| Script | Description |
| :--- | :--- |
| `exp_extensibility.m` | Verifies drop-in compatibility with first-order SAV (1D Allen–Cahn, $2.2\times10^{-16}$ equivalence, standard breaks at $T_\ast=44.00$ vs. isotypic $D\equiv0$), 2D Conservative Allen–Cahn ($N=256$ X-junction: standard breaks at $T_\ast=0.220$ despite exact mass conservation $|\mathrm{mean}|\le 1.5\times10^{-17}$ vs. isotypic $D\equiv0$), and 2D fourth-order Cahn–Hilliard ($N=192$, $5.6\times10^{-16}$ equivalence, isotypic $D\equiv0$) |
