# Python Drop-In Backend (NumPy + SciPy)

[`isotypic_spectral.py`](isotypic_spectral.py) implements a self-contained, zero-overhead `IsotypicBackend(d, N, mode)` class using only `numpy` and `scipy.fft` (`dct`/`dst` Type I for `'AR'`; twiddled `fft`/`ifft` for `'T2'`).

## Usage
```python
from isotypic_spectral import IsotypicBackend
B = IsotypicBackend(d=1, N=128, mode="AR")  # or mode="T2" for arbitrary shifts
U = B.split(u0)
Uh = B.to_spec(U)
Nu = B.to_spec(U - B.cube(U))               # or B.orbit(U, F_odd)
```

Run `python demo_scipy_etdrk4.py` to test both `scipy.integrate.solve_ivp(..., method="BDF")` and irrational-shift $\mathcal T_2$ preservation in pure `float64`.
