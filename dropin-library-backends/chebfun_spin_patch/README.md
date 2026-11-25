# Drop-In Patch for Chebfun (`spin`, `spin2`, `spin3`)

Chebfun's periodic PDE drivers (`spin`, `spin2`, `spin3`, built around `@spinop` / `@spinop2` / `@spinop3` and Cox–Matthews ETDRK4 in `spinscheme/@etdrk4/oneStep.m`) compute every time step using three primitives:
1. Diagonal Fourier multipliers `E`, `Q`, `f1`, `f2`, `f3` on the linear symbol `L(k) = -eps^2 |k|^2`,
2. Forward/inverse `fft`/`ifft` (or `fft2`/`ifft2`, `fftn`/`ifftn`),
3. Pointwise evaluation of the nonlinear operator `N(u) = u - u.^3` in physical space.

## Minimal 3-Line Drop-In Replacement
To make Chebfun's `spin` / `spin2` / `spin3` **bitwise symmetry-preserving** (without filters, zero-mode resets, or accuracy loss on non-symmetric initial data):

1. Add `matlab_ode_and_spectral/iso_spectral_backend.m` to your MATLAB path.
2. Initialize the channel backend `B = iso_spectral_backend(d, N, 'AR')` (or `'T2'` for unknown arbitrary shifts).
3. Inside the ETDRK4 nonlinear evaluation callback (`oneStep.m` / `computeNonlin`):
   ```matlab
   % BEFORE (standard Chebfun FFT evaluation — leaks O(1e-16) across channels):
   %   u  = real(ifftn(uh));
   %   Nu = fftn(u - u.^3);

   % AFTER (drop-in isotypic channel backend — bitwise D === 0, identical algebra):
   U  = B.toPhys(Uh);
   Nu = B.toSpec(U - B.cube(U));   % or B.orbit(U, F_odd) for non-polynomial odd F
   ```
   All ETDRK4 contour-integral weights (`E`, `E2`, `Q`, `f1`, `f2`, `f3`) are evaluated verbatim on `B.K2` (exact same eigenvalues $|k|^2$)!

## Standalone Verification (No Chebfun Installation Required)
Run [`spin_iso_etdrk4_demo.m`](spin_iso_etdrk4_demo.m) to reproduce Chebfun `spin`/`spin2`'s exact Cox–Matthews ETDRK4 loop on 1D ($u_0 = \sin x$, $T = 200$) and 2D X-junction ($N = 256$, $T = 5$), producing [`fig_demo_chebfun_spin.png`](fig_demo_chebfun_spin.png).
