# Drop-In Isotypic Spectral Backends for Popular PDE Libraries

This folder provides minimal, ready-to-use **drop-in backend modules, patches, and worked examples** showing how to embed our floating-point-equivariant symmetry-channel representation into existing scientific-computing libraries **without modifying any time-stepping algorithm** and **with virtually zero runtime penalty on general non-symmetric problems** (and up to **2x–7x speedup** on symmetric ones).

## The Universal 3-Line Drop-In Rule

To upgrade *any* Fourier (or DCT/DST) pseudo-spectral PDE solver so that it preserves symmetries bitwise without filters, zero-mode resets, or accuracy loss on non-symmetric problems:

1. **Replace full-grid `fft`/`ifft`** with channel-wise transforms (`toSpec` / `toPhys` on $A_{\mathrm R}$ or $\mathcal{T}_2$ channels).
2. **Replace pointwise `u.^3` (or odd `F(u)`)** with the exact $\mathbb{Z}_2^d$-graded channel product `cube(U)` (or odd-exact Walsh–Hadamard orbit evaluation `orbit(U, F)`).
3. **Keep everything else untouched**: all diagonal multipliers $m(|\mathbf{k}|^2)$, ODE/PDE time integrators (`ode15s`, `solve_ivp`, ETDRK, IMEX, SAV, SBDF, SDC, ESP, CN/cG), and adaptive step controllers remain 100% unchanged!

---

## Available Library Adapters & Verified Runtime Efficiency

| Library / Framework | Subfolder | Key Files | Verified Accuracy & Efficiency (Symmetric & Non-Symmetric) |
| :--- | :--- | :--- | :--- |
| **MATLAB built-in `ode15s` & spectral scripts** | [`matlab_ode_and_spectral/`](matlab_ode_and_spectral/) | `iso_spectral_backend.m`, `demo_matlab_ode15s.m` | **Symmetric ($T=150$):** standard `fft` breaks at $T_\ast=27$ (`0.493 s`, 2,254 steps); drop-in keeps **$D \equiv 0$** in **`0.067 s` (`0.14x`, ~7x faster!)**. **Classic non-symmetric ($T=1.5$, all channels active):** `0.043 s` vs. `0.048 s` (**`0.88x`**), $\|u_{\mathrm{std}}-u_{\mathrm{iso}}\|_\infty = 2.0\times 10^{-10}$ (matching `RelTol=1e-10`). |
| **Chebfun (`spin` / `spin2` / `spin3`)** | [`chebfun_spin_patch/`](chebfun_spin_patch/) | `README.md`, `spin_iso_etdrk4_demo.m` | **1D symmetric ($T=200$):** `1.05x–1.07x` (`D===0`). **2D symmetric X-junction ($N=256$):** **`0.52x` (`1.63 s` vs. `3.12 s`, ~2x faster!, `D===0`)**. **Classic 2D non-symmetric coarsening ($N=256$, all 4 channels active):** `T2` takes **`1.08x`** (`6.48` vs. `5.98 ms/step`, `3.5e-14` agreement); MATLAB-emulated `AR` takes `1.46x` (`~1.0x` with native FFTW DCT/DST). |
| **Python (NumPy + SciPy `solve_ivp` / `scipy.fft`)** | [`python_numpy_scipy/`](python_numpy_scipy/) | `isotypic_spectral.py`, `demo_scipy_etdrk4.py` | Self-contained 1D/2D/3D Python class (`IsotypicBackend`) using native `scipy.fft.dct`/`dst` (Type I) for $A_{\mathrm R}$ and twiddled FFTs for $\mathcal T_2$. |
| **Dedalus v3 & Shenfun (Python)** | [`python_dedalus_shenfun/`](python_dedalus_shenfun/) | `dedalus_v3_isotypic_allen_cahn.py`, `shenfun_isotypic_allen_cahn.py` | Exact channel-equation formulation using Dedalus's native `SBDF4` and Shenfun's ETDRK2. |
| **Custom Research Solvers (MATLAB / Python)** | [`custom_solver_template/`](custom_solver_template/) | `README.md` | 5-minute retrofit checklist for IMEX, SAV, CN/cG, ESP, and multi-field PDEs. |

---

## Live Demonstration Figures (Accuracy + Wall-Clock Cost)

### 1. MATLAB Built-In `ode15s` (`matlab_ode_and_spectral/demo_matlab_ode15s.m`)
Changing only 3 lines in the `ode15s` RHS prevents the $T_\ast = 27$ collapse (~7x faster on the symmetric run) while matching standard `fft` to `RelTol = 1e-10` at nearly identical runtime on a classic non-symmetric problem:

![MATLAB ode15s Drop-In Comparison](matlab_ode_and_spectral/fig_demo_ode15s.png)

### 2. Chebfun `spin` / `spin2` ETDRK4 Style (`chebfun_spin_patch/spin_iso_etdrk4_demo.m`)
Exact same Cox–Matthews ETDRK4 time stepper on 1D symmetric ($T=200$), 2D symmetric X-junction ($N=256$, ~2x faster!), and a classic 2D non-symmetric spinodal coarsening benchmark ($N=256$, all 4 channels active, `1.08x` cost and `3.5e-14` agreement):

![Chebfun spin/spin2 Drop-In Comparison](chebfun_spin_patch/fig_demo_chebfun_spin.png)

---

## Intellectual Property & Citation

Copyright © 2026 **Minghui Li** and **Zixian Liu** (School of Mathematical Sciences, Jiangxi Science and Technology Normal University), released under the [MIT License](../LICENSE).

All drop-in backend modules (`iso_spectral_backend.m`, `isotypic_spectral.py`, Dedalus v3 / Shenfun / Chebfun recipes) are freely usable in academic and open-source solvers provided that the original copyright notice is retained. **If you embed or adapt any of these backends in your code, toolbox, or paper, please cite our accompanying manuscript**:

```bibtex
@article{LiLiu2026Isotypic,
  author  = {Minghui Li and Zixian Liu},
  title   = {Exact long-time symmetry preservation in phase-field simulations:
             isotypic spectral schemes that are equivariant in floating-point arithmetic},
  year    = {2026},
  note    = {Preprint; code available at
             \url{https://github.com/minghuili123456789/isotypic-phase-field}}
}
```
