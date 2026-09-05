"""Exact six-row extension of the frozen 0085 analytic operator witness."""

import contextlib
import io
from pathlib import Path
import runpy

import numpy as np

source = Path(__file__).resolve().parents[1] / "0085" / "jet_probe.py"
with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(source))
indices = old["indices"]
generator = np.zeros_like(old["angular"])
for axis in range(3):
    e = [int(i == axis) for i in range(3)]
    for j in range(3):
        for col, alpha in enumerate(indices):
            if sum(alpha) == 1:
                ell = alpha.index(1)
                direction = [int(i == ell) for i in range(3)]
                generator[axis, j * len(indices) + col] = -old["cross"](e, direction)[j]
rank = old["modular_rank"]
A = old["matrix"]
L = old["angular"]
print(f"constraint rank={rank(A)}; analytic universal upper bound=235")
print(f"physical-spin augmented rank={rank(np.vstack([A, L]))}")
print(f"generator-moment augmented rank={rank(np.vstack([A, generator]))}")
print(f"all-six moment augmented rank={rank(np.vstack([A, L, generator]))}")
print("All ranks are exact over F_101; no floating sign or tolerance.")
