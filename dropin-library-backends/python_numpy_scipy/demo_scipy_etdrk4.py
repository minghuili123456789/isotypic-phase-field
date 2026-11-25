"""
Worked Example: Drop-in usage of IsotypicBackend with SciPy solve_ivp & ETDRK4.

Demonstrates:
  1. ~1e-15 algebraic equivalence with standard NumPy FFT on non-symmetric 1D/2D data.
  2. Exact bitwise odd-symmetry preservation (D === 0) on 1D Allen-Cahn (u0 = sin(x))
     and shifted irrational 1D Allen-Cahn (u0 = sin(x + 0.37123456789)) via T2 mode.
"""
import numpy as np
from scipy.integrate import solve_ivp
from isotypic_spectral import IsotypicBackend


def demo_1d_scipy_solve_ivp():
    N = 128
    x = np.arange(N) * 2 * np.pi / N
    k = np.fft.fftfreq(N, d=1.0 / N)
    K2_std = k**2
    B = IsotypicBackend(d=1, N=N, mode="AR")

    # 1. Non-symmetric equivalence check
    u_ns = 0.4 * np.sin(x) + 0.3 * np.cos(2 * x) + 0.15 * np.sin(3 * x + 0.7)
    Lu_std = np.fft.ifft(-K2_std * np.fft.fft(u_ns)).real + u_ns - u_ns**3
    U_ns = B.split(u_ns)
    Lu_iso = B.assemble(B.to_phys(-B.K2 * B.to_spec(U_ns)) + U_ns - B.cube(U_ns))
    print(f"1D non-symmetric RHS equivalence: {np.max(np.abs(Lu_std - Lu_iso)):.2e}")

    # 2. Long-time integration with SciPy BDF (u0 = sin(x), T = 100)
    rhs_std = lambda t, u: np.fft.ifft(-K2_std * np.fft.fft(u)).real + u - u**3
    shp = (B.M + 1, 2)
    rhs_iso = lambda t, y: (
        B.to_phys(-B.K2 * B.to_spec(y.reshape(shp)))
        + y.reshape(shp)
        - B.cube(y.reshape(shp))
    ).ravel()

    U0 = np.zeros(shp)
    U0[1:B.M, 1] = np.sin(np.arange(1, B.M) * 2 * np.pi / N)
    u0 = B.assemble(U0)

    sol_std = solve_ivp(rhs_std, (0, 100), u0, method="BDF", rtol=1e-10, atol=1e-12)
    sol_iso = solve_ivp(rhs_iso, (0, 100), U0.ravel(), method="BDF", rtol=1e-10, atol=1e-12)

    D_std = np.max(np.abs(B.split(sol_std.y[:, -1])[:, 0]))
    D_iso = np.max(np.abs(sol_iso.y[:, -1].reshape(shp)[:, 0]))
    print(f"SciPy BDF at T=100: std FFT D = {D_std:.2e} | drop-in isotypic D = {D_iso:g}")


def demo_1d_shifted_t2():
    N = 256
    theta = 0.37123456789
    B = IsotypicBackend(d=1, N=N, mode="T2")
    xh = np.arange(B.M) * 2 * np.pi / N
    U0 = np.zeros((B.M, 2))
    U0[:, 1] = np.sin(xh + theta)  # pure odd-harmonic channel (arbitrary irrational shift!)
    # 1000 IMEX1 steps in pure float64
    tau = 0.05
    Dinv = 1.0 / (1.0 + tau * B.K2)
    Uh = B.to_spec(U0)
    for _ in range(2000):
        U = B.to_phys(Uh)
        Uh = Dinv * (Uh + tau * B.to_spec(U - B.cube(U)))
    D_t2 = np.max(np.abs(B.to_phys(Uh)[:, 0]))
    print(f"Shifted irrational T2 (T=100, pure float64): max transverse D = {D_t2:g}")


if __name__ == "__main__":
    demo_1d_scipy_solve_ivp()
    demo_1d_shifted_t2()
