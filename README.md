# Isotypic Spectral Schemes for Exact Long-Time Symmetry Preservation in Phase-Field Simulations

Official code repository for the manuscript:

> **Minghui Li\*, Zixian Liu**, *"Exact long-time symmetry preservation in phase-field simulations: isotypic spectral schemes that are equivariant in floating-point arithmetic"*, School of Mathematical Sciences, Jiangxi Science and Technology Normal University, Nanchang, Jiangxi, China (`minghuili_xmu@163.com`).

## Overview

Standard Fourier pseudo-spectral schemes (IMEX, ETDRK2/4, Crank–Nicolson, ESP, SAV, SBDF, IMEX-RK, SDC, operator splitting) preserve odd and reflection symmetries in exact arithmetic, yet lose them in floating-point arithmetic due to non-equivariant FFT rounding, collapsing to spurious steady states at a universal breakdown time $T_\ast = \lambda_\perp^{-1}\log(1/\mathbf{u}) + O(1)$.

This repository provides a **non-interventional, drop-in symmetry-channel (isotypic) spectral backend** that:
- is **algebraically identical** to the standard Fourier pseudo-spectral method on arbitrary non-symmetric data ($10^{-16}\text{--}10^{-15}$ agreement);
- preserves every reflection, point-inversion ($C_2$), half-period anti-periodic ($\mathcal{T}_2$), and product ($A_{\mathrm R}\times\mathcal{T}_2$) symmetry **bitwise ($D \equiv 0$)** with zero filters, projections, or zero-mode resets;
- propagates genuine sub-roundoff asymmetries faithfully with relative accuracy down to $\delta = 10^{-60}$ (validated against 120-digit VPA);
- handles **unknown arbitrary (even irrational) spatial shifts** in 100% pure IEEE double-precision time stepping via $\mathcal{T}_2$ and **CF-ISM**.

## Repository Structure

| Folder | Description |
| :--- | :--- |
| [`core/`](core/) | Core MATLAB isotypic spectral library ($A_{\mathrm R}$, $\mathcal{T}_2$, Walsh–Hadamard orbit evaluation, CF-ISM, shared integrators) |
| [`experiments/`](experiments/) | All numerical experiments and figure scripts from the manuscript, organized by section (`01`–`05`) |
| [`supplementary-bc-evidence/`](supplementary-bc-evidence/) | **Supplementary Boundary-Condition Evidence**: runnable MATLAB code and comparison figures proving standard DCT/DST solvers fail on **1D/2D homogeneous Neumann** and **1D homogeneous Dirichlet** Allen–Cahn problems while our half-domain isotypic channel backend preserves $D \equiv 0$ bitwise |
| [`dropin-library-backends/`](dropin-library-backends/) | **Ready-to-Use Drop-In Adapters**: minimal replacement modules, worked demos, and output plots for **MATLAB (`ode15s` / spectral)**, **Chebfun (`spin`/`spin2`/`spin3`)**, **Python (NumPy/SciPy)**, **Dedalus v3**, **Shenfun**, and **custom research solvers** |

## Quick Start (MATLAB)

```matlab
addpath('core');
% 1D reflection isotypic backend (drop-in replacement for FFT)
P = sp_setup(struct('d',1,'N',128,'eps2',1,'gam',1,'S',0,'rep','iso'));
O = sp_ops(P);
I = sp_integrator(P, O, 'etdrk4', 0.02);
U = sp_init(P, @(x) sin(x));
vh = O.toSpec(U);
for n = 1:50000, vh = I.step(vh); end
V = O.toPhys(vh);
fprintf('Max symmetry defect at T = 1000: %g\n', max(abs(V(:,1)))); % exact 0!
```

See [`experiments/README.md`](experiments/README.md) to reproduce every figure and table in the paper.

## License
MIT License — see [`LICENSE`](LICENSE).
