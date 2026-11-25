# Drop-In Recipes for Dedalus v3 & Shenfun

- [`dedalus_v3_isotypic_allen_cahn.py`](dedalus_v3_isotypic_allen_cahn.py): Declares the even ($u_0$, Neumann–Neumann on $[0,\pi]$) and odd ($u_1$, Dirichlet–Dirichlet on $[0,\pi]$) symmetry channels with the exact $\mathbb Z_2$-graded cubic equations and Dedalus v3's native `d3.SBDF4` integrator untouched.
- [`shenfun_isotypic_allen_cahn.py`](shenfun_isotypic_allen_cahn.py): Drops `IsotypicBackend` (`../python_numpy_scipy/isotypic_spectral.py`) directly into a 2D Shenfun-style ETDRK2 solver for the odd–odd X-junction ($D \equiv 0$ bitwise).
