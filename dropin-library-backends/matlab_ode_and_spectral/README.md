# MATLAB Drop-In Backend (`ode15s`, `ode45`, and Spectral Solvers)

Drop-in module [`iso_spectral_backend.m`](iso_spectral_backend.m) providing `toSpec`, `toPhys`, `cube`, `orbit`, `split`, `assemble`, and `K2` for both reflection (`'AR'`) and half-period parity (`'T2'`) backends.

## Live Demo & Wall-Clock Efficiency: MATLAB Built-In `ode15s`

Run [`demo_matlab_ode15s.m`](demo_matlab_ode15s.m) (`N = 128`, `RelTol = 1e-10, AbsTol = 1e-12`, JIT-warmed):

| Test Case | Accuracy / Symmetry | Standard `fft` RHS | Drop-In `iso_spectral_backend` RHS | Speed Ratio |
| :--- | :--- | :--- | :--- | :---: |
| **Classic Non-Symmetric 1D** ($T = 1.5$, all channels active) | $\|u_{\mathrm{std}} - u_{\mathrm{iso}}\|_\infty = 2.0\times 10^{-10}$ (matches `RelTol`) | `0.048 s` (`412` steps, `818` fevals) | `0.043 s` (`340` steps, `658` fevals) | **`0.88x`** (on par) |
| **Symmetric Odd 1D** ($u_0 = \sin x$, $T = 150$) | Standard collapses at **$T_\ast = 27$**; drop-in keeps **$D \equiv 0$ bitwise** | `0.493 s` (`2,254` steps, `22,463` fevals) | `0.067 s` (`646` steps, `1,383` fevals) | **`0.14x`** (~7x faster) |

**Why is the drop-in backend equally fast on general non-symmetric problems and even faster on symmetric ones?**
- On **general non-symmetric problems** (all channels nonzero), the DOF count ($2(M+1) = N+2$ vs. $N$) and per-step transform cost match standard `fft`, giving nearly identical wall-clock runtime and identical trajectories.
- On **symmetric problems**, the non-interventional exact-zero channel short-circuit avoids transforming identically zero channels, and eliminating the spurious symmetry-breaking transition prevents `ode15s` from wasting thousands of adaptive steps and Jacobian evaluations tracking a false collapse.

![MATLAB ode15s Comparison](fig_demo_ode15s.png)
