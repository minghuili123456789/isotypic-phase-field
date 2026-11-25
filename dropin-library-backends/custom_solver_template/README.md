# 5-Minute Retrofit Guide for Custom Research Solvers

If your group maintains a custom periodic, Neumann, or Dirichlet spectral solver (IMEX, ETDRK, CN/cG, SAV, SBDF, SDC, or ESP), follow this 3-step checklist:

## Step 1: Choose the Symmetry-Channel Backend
- **Centred reflection or $C_2$ point symmetry** $\to$ `mode = 'AR'` (DCT-I/DST-I on $[0,\pi]^d$).
- **Unknown arbitrary / irrational shift with odd harmonics (anti-periodic)** $\to$ `mode = 'T2'` (pure double precision, zero shift detection).
- **Combined reflection + $\pi$-periodicity** (e.g. $\sin(2x)\sin(2y)$) $\to$ $A_{\mathrm R}\times\mathcal{T}_2$ (see `experiments/03_2d_3d_centred_and_wang41/exp_wang41_fullgroup.m`).
- **Unknown irrational shift with mixed even+odd harmonics** $\to$ CF-ISM (`core/cf_gauge_init.m` once at $t=0$, then `mode = 'AR'` in pure double for all $t>0$).
- **Homogeneous Neumann / Dirichlet on $[0,\pi]^d$** $\to$ half-domain $[0,\pi/2]^d$ midpoint channels (see `supplementary-bc-evidence/`).

## Step 2: Swap Only the Spatial Primitives (Before vs. After)

### MATLAB
```matlab
% --- BEFORE (standard FFT solver) ---
% uh = fftn(u);
% Nu = fftn(u - u.^3);
% uh_new = mult(K2) .* uh + coeff(K2) .* Nu;

% --- AFTER (drop-in isotypic backend) ---
B  = iso_spectral_backend(d, N, 'AR');   % or 'T2'
Uh = B.toSpec(B.split(u0));
% Inside every time step:
U  = B.toPhys(Uh);
Nu = B.toSpec(U - B.cube(U));            % or B.orbit(U, F_odd) for non-polynomial F
Uh_new = mult(B.K2) .* Uh + coeff(B.K2) .* Nu;
```

### Python (NumPy / SciPy)
```python
from isotypic_spectral import IsotypicBackend

B  = IsotypicBackend(d=2, N=256, mode="AR")  # or mode="T2"
Uh = B.to_spec(B.split(u0))
# Inside every time step:
U  = B.to_phys(Uh)
Nu = B.to_spec(U - B.cube(U))                # or B.orbit(U, F_odd)
Uh_new = mult(B.K2) * Uh + coeff(B.K2) * Nu
```

## Step 3: Scalar Inner Products (for SAV / Energy / Mass Conservation)
Because channels are mutually orthogonal, any discrete $L^2$ inner product $\langle a, b \rangle$ is evaluated channel-wise with trapezoidal weights `[1, 2, ..., 2, 1]` on $[0,\pi]^d$ (for `AR`) or uniform weight `2` on $[0,\pi)^d$ (for `T2`):
- In MATLAB: `ip = h * sum(w .* sum(A .* B, 2));` (1D `AR`, exact to machine precision).
- Every scalar multiplier (e.g., SAV auxiliary variable $r^{n+1}$) multiplies all channels diagonally, preserving $0 \cdot r^{n+1} = 0$ bitwise!

## Verification Checklist
1. Run 20 steps on a **random non-symmetric initial field**: verify `max|u_std - B.assemble(U_iso)| ~ 1e-15`.
2. Run long-time integration on your **symmetric initial field**: verify all transverse channels stay **exact `0.0`** (`maxD == 0`).
