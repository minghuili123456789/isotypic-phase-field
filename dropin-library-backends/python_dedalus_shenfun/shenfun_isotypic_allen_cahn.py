"""
Drop-in Isotypic Symmetry-Channel Formulation for Shenfun (https://github.com/spectralDNS/shenfun)

Replaces Shenfun's single Fourier FunctionSpace with two orthogonal parity spaces
on [0, pi] (using Shenfun's built-in Neumann/Dirichlet bases or SciPy DCT-I/DST-I),
leaving Shenfun's ETDRK / IMEX / BDF integrators completely unchanged:
    u = u0 (even channel) + u1 (odd channel)
    (u^3)_0 = u0**3 + 3*u0*u1**2
    (u^3)_1 = 3*u0**2*u1 + u1**3

Note: For periodic Fourier spectral methods in Shenfun (d = 1, 2, 3), the fastest
drop-in replacement is importing `IsotypicBackend` from `../python_numpy_scipy/isotypic_spectral.py`
directly inside your Shenfun nonlinear RHS callback `LinearRHS` / `NonlinearRHS`!
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python_numpy_scipy"))
from isotypic_spectral import IsotypicBackend


def shenfun_style_etdrk2_2d(N=128, eps=0.05, tau=1e-3, T=1.0, mode="AR"):
    """
    Minimal 2D Allen-Cahn drop-in showing how a Shenfun-style spectral solver
    swaps only `forward`/`backward` and `u**3` to `IsotypicBackend`.
    """
    gam = 1.0 / eps**2
    S = 2.0 * gam
    B = IsotypicBackend(d=2, N=N, mode=mode)
    Lk = B.K2 + S
    z = -tau * Lk
    E = np.exp(z)
    z_safe = np.where(np.abs(z) > 1e-6, z, 1.0)
    Q1 = tau * np.where(np.abs(z) > 1e-6, (E - 1.0) / z_safe, 1.0 + 0.5 * z)
    Q2 = tau * np.where(np.abs(z) > 1e-6, (E - 1.0 - z) / z_safe**2, 0.5 + z / 6.0)

    def nonlin(Uh):
        U = B.to_phys(Uh)
        return B.to_spec((gam + S) * U - gam * B.cube(U))

    # 2D odd-odd X-junction initial datum (channel 3 = (1, 1))
    xh = np.arange(B.M + 1) * 2 * np.pi / N if mode == "AR" else np.arange(B.M) * 2 * np.pi / N
    X, Y = np.meshgrid(xh, xh, indexing="ij")
    U0 = np.zeros(B.K2.shape)
    U0[..., 3] = np.tanh(np.sin(X) * np.sin(Y) / (np.sqrt(2.0) * eps))
    if mode == "AR":
        U0[[0, -1], :, 3] = 0.0
        U0[:, [0, -1], 3] = 0.0

    Uh = B.to_spec(U0)
    nst = int(round(T / tau))
    for _ in range(nst):
        Nu = nonlin(Uh)
        a = E * Uh + Q1 * Nu
        Uh = a + Q2 * (nonlin(a) - Nu)

    V = B.to_phys(Uh)
    maxD = float(np.max(np.abs(V[..., :3])))
    print(f"Shenfun-style 2D X-junction ({mode}, N={N}, T={T}): max transverse defect D = {maxD:g}")


if __name__ == "__main__":
    shenfun_style_etdrk2_2d()
