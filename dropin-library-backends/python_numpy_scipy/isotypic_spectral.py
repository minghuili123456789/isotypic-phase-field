"""
Drop-in floating-point-equivariant isotypic spectral backend for Python (NumPy + SciPy).

Supports:
  - Reflection isotypic group A_R (1D, 2D, 3D via scipy.fft.dct/dst Type I)
  - Half-period translation group T_2 (1D, 2D, 3D via twiddled FFTs; preserves
    anti-periodic / odd-harmonic symmetries for ANY unknown irrational shift!)
  - Exact Z_2^d-graded cubic nonlinearity AND odd-exact Walsh-Hadamard orbit
    evaluation for general odd nonlinearities (e.g., Flory-Huggins log potential).

Algebraic identity: on arbitrary non-symmetric data, every operator matches
standard NumPy/SciPy Fourier pseudo-spectral operators to ~1e-15 while keeping
every symmetric subspace bitwise invariant (D === 0) without filters or resets.

Copyright (c) 2026 Minghui Li, Zixian Liu (MIT License).
Please cite: M. Li, Z. Liu, "Exact long-time symmetry preservation in
phase-field simulations: isotypic spectral schemes that are equivariant in
floating-point arithmetic" (2026), https://github.com/minghuili123456789/isotypic-phase-field
"""
import numpy as np
from scipy.fft import dct, dst, fft, ifft


class IsotypicBackend:
    def __init__(self, d: int, N: int, mode: str = "AR"):
        """
        Parameters
        ----------
        d : int
            Spatial dimension (1, 2, or 3).
        N : int
            Periodic grid size per axis (even for 'AR', divisible by 4 for 'T2').
        mode : {'AR', 'T2'}
            'AR': reflection isotypic (DCT-I / DST-I on [0, pi]^d)
            'T2': half-period translation isotypic (twiddled FFT on [0, pi)^d)
        """
        self.d = d
        self.N = N
        self.M = N // 2
        self.nc = 1 << d
        self.mode = mode.upper()

        if self.mode == "AR":
            k1 = np.arange(self.M + 1, dtype=float)
            grids = np.meshgrid(*([k1] * d), indexing="ij")
            k2 = sum(g**2 for g in grids)
            self.K2 = np.repeat(k2[..., None], self.nc, axis=-1)
        elif self.mode == "T2":
            assert N % 4 == 0, "N must be divisible by 4 for T2 mode."
            j = np.arange(self.M, dtype=float)
            self.tw = np.exp(-1j * np.pi * j / self.M)
            base = np.concatenate([np.arange(self.M // 2), np.arange(-self.M // 2, 0)])
            ke = 2.0 * base
            ko = 2.0 * base + 1.0
            K1 = (ke, ko)
            self.K2 = np.zeros((self.M,) * d + (self.nc,), dtype=float)
            for c in range(self.nc):
                axes_k = [K1[(c >> ax) & 1] for ax in range(d)]
                grids = np.meshgrid(*axes_k, indexing="ij")
                self.K2[..., c] = sum(g**2 for g in grids)
        else:
            raise ValueError("mode must be 'AR' or 'T2'")

    # ------------------------------------------------------------------
    # Transforms (replace standard np.fft.fftn / ifftn)
    # ------------------------------------------------------------------
    def to_spec(self, U: np.ndarray) -> np.ndarray:
        if self.mode == "AR":
            W = np.zeros_like(U, dtype=float)
            for c in range(self.nc):
                A = U[..., c].copy()
                for ax in range(self.d):
                    if (c >> ax) & 1:
                        # DST-I on interior points 1..M-1 (exact zeros at 0 and M)
                        sl = [slice(None)] * self.d
                        sl[ax] = slice(1, self.M)
                        B = np.zeros_like(A)
                        B[tuple(sl)] = dst(A[tuple(sl)], type=1, axis=ax)
                        A = B
                    else:
                        A = dct(A, type=1, axis=ax)
                W[..., c] = A
            return W
        else:
            W = np.zeros(U.shape, dtype=np.complex128)
            for c in range(self.nc):
                A = U[..., c].astype(np.complex128, copy=True)
                for ax in range(self.d):
                    if (c >> ax) & 1:
                        shp = [1] * self.d
                        shp[ax] = self.M
                        A = fft(A * self.tw.reshape(shp), axis=ax)
                    else:
                        A = fft(A, axis=ax)
                W[..., c] = A
            return W

    def to_phys(self, Uh: np.ndarray) -> np.ndarray:
        if self.mode == "AR":
            W = np.zeros_like(Uh, dtype=float)
            scale = 1.0 / (2.0 * self.M)
            for c in range(self.nc):
                A = Uh[..., c].copy()
                for ax in range(self.d):
                    if (c >> ax) & 1:
                        sl = [slice(None)] * self.d
                        sl[ax] = slice(1, self.M)
                        B = np.zeros_like(A)
                        B[tuple(sl)] = dst(A[tuple(sl)], type=1, axis=ax) * scale
                        A = B
                    else:
                        A = dct(A, type=1, axis=ax) * scale
                W[..., c] = A
            return W
        else:
            W = np.zeros(Uh.shape, dtype=float)
            for c in range(self.nc):
                A = Uh[..., c].copy()
                for ax in range(self.d):
                    if (c >> ax) & 1:
                        shp = [1] * self.d
                        shp[ax] = self.M
                        A = ifft(A, axis=ax) * np.conj(self.tw).reshape(shp)
                    else:
                        A = ifft(A, axis=ax)
                W[..., c] = A.real
            return W

    # ------------------------------------------------------------------
    # Nonlinearities (replace pointwise u**3 or F(u))
    # ------------------------------------------------------------------
    def cube(self, U: np.ndarray) -> np.ndarray:
        """Exact Z_2^d-graded cubic nonlinearity (u^3)_sigma."""
        A = U.reshape(-1, self.nc)
        A2 = np.zeros_like(A)
        A3 = np.zeros_like(A)
        for r in range(self.nc):
            for a in range(self.nc):
                A2[:, r] += A[:, a] * A[:, a ^ r]
        for s in range(self.nc):
            for r in range(self.nc):
                A3[:, s] += A2[:, r] * A[:, s ^ r]
        return A3.reshape(U.shape)

    def orbit(self, U: np.ndarray, F_odd) -> np.ndarray:
        """Odd-exact Walsh-Hadamard orbit evaluation for general odd F(u)."""
        V = U.reshape(-1, self.nc).copy()
        for i in range(self.d):
            step = 1 << i
            W = V.copy()
            for c in range(self.nc):
                if not (c & step):
                    p, q = V[:, c], V[:, c | step]
                    W[:, c] = p + q
                    W[:, c | step] = p - q
            V = W
        V = F_odd(V)
        for i in range(self.d):
            step = 1 << i
            W = V.copy()
            for c in range(self.nc):
                if not (c & step):
                    p, q = V[:, c], V[:, c | step]
                    W[:, c] = 0.5 * (p + q)
                    W[:, c | step] = 0.5 * (p - q)
            V = W
        return V.reshape(U.shape)

    # ------------------------------------------------------------------
    # Split / assemble helpers between full N^d periodic grid and channels
    # ------------------------------------------------------------------
    def split(self, u: np.ndarray) -> np.ndarray:
        if self.mode == "AR":
            refl = np.concatenate([[0], np.arange(self.N - 1, 0, -1)])
            vals = [u.astype(float)]
            # Iterate axes from d-1 down to 0 so axis `ax` maps to bit `ax` in channel index c
            for ax in range(self.d - 1, -1, -1):
                nv = []
                for v in vals:
                    vr = np.take(v, refl, axis=ax)
                    sl = [slice(None)] * self.d
                    sl[ax] = slice(0, self.M + 1)
                    ve = 0.5 * (v + vr)[tuple(sl)]
                    vo = 0.5 * (v - vr)[tuple(sl)]
                    sl0 = [slice(None)] * self.d; sl0[ax] = 0
                    slM = [slice(None)] * self.d; slM[ax] = self.M
                    vo[tuple(sl0)] = 0.0
                    vo[tuple(slM)] = 0.0
                    nv.extend([ve, vo])
                vals = nv
            return np.stack(vals, axis=-1)
        else:
            vals = [u.astype(float)]
            for ax in range(self.d - 1, -1, -1):
                nv = []
                for v in vals:
                    sl1 = [slice(None)] * self.d; sl1[ax] = slice(0, self.M)
                    sl2 = [slice(None)] * self.d; sl2[ax] = slice(self.M, self.N)
                    nv.append(0.5 * (v[tuple(sl1)] + v[tuple(sl2)]))
                    nv.append(0.5 * (v[tuple(sl1)] - v[tuple(sl2)]))
                vals = nv
            return np.stack(vals, axis=-1)

    def assemble(self, U: np.ndarray) -> np.ndarray:
        u = np.zeros((self.N,) * self.d, dtype=float)
        if self.mode == "AR":
            idx = np.concatenate([np.arange(self.M + 1), np.arange(self.M - 1, 0, -1)])
            sgn = np.concatenate([np.ones(self.M + 1), -np.ones(self.M - 1)])
        else:
            idx = np.concatenate([np.arange(self.M), np.arange(self.M)])
            sgn = np.concatenate([np.ones(self.M), -np.ones(self.M)])
        for c in range(self.nc):
            A = U[..., c]
            for ax in range(self.d):
                A = np.take(A, idx, axis=ax)
                if (c >> ax) & 1:
                    shp = [1] * self.d; shp[ax] = self.N
                    A = A * sgn.reshape(shp)
            u += A
        return u
