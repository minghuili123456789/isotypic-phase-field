# Core Isotypic Spectral Library (`core/`)

Reusable MATLAB implementation of the non-interventional isotypic spectral backends described in Section 4 of the manuscript:

| File | Description |
| :--- | :--- |
| `sp_setup.m`, `sp_ops.m`, `sp_r2r.m` | Reflection isotypic ($A_{\mathrm R}$, Type I DCT-I/DST-I) and standard Fourier pseudo-spectral operators in 1D, 2D, and 3D |
| `sp_init.m`, `sp_channels.m`, `sp_assemble.m` | Exact bitwise channel initialization, full-to-channel splitting, and channel-to-full assembly |
| `iso_orbit.m` | Odd-exact Walsh–Hadamard orbit evaluation for general non-polynomial odd potentials (Proposition 4.4) |
| `iso_transform_mat.m` | Matrix-form DCT-I/DST-I helper for verification |
| `sp_t2_setup.m`, `sp_t2_ops.m`, `sp_t2_init.m`, `sp_t2_assemble.m` | Half-period translation isotypic ($\mathcal T_2$, wavenumber-parity twiddled FFTs) in 1D, 2D, and 3D |
| `cf_gauge_init.m` | Canonical-Frame Isotypic Spectral Method (CF-ISM, Algorithm 2): 1-FFT circular auto-convolution + $d_0$-digit alias-free Gauss–Newton gauge polish at $t=0$ |
| `sp_integrator.m`, `sp_run.m` | Shared time-stepping drivers (`imex1`, `etdrk2`, `etdrk4`, `cn`, `esp`, `zmode`) |
| `fig_common.m` | Unified LaTeX figure styling defaults |

## Intellectual Property & Citation
Copyright © 2026 **Minghui Li** and **Zixian Liu** (MIT License; see [`../LICENSE`](../LICENSE)). When using or embedding any file from `core/` in research or software, please retain the copyright header and cite our paper (see [`../CITATION.cff`](../CITATION.cff) or the BibTeX entry in [`../README.md`](../README.md)).
