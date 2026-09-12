"""A1 archive: H5 d(t)/a exclusion certificate. Run: python3 run_a1.py. CPython 3.12.2, numpy."""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '0108-cipher-radical', 'receipts', 'poc2-filament'))
from run_poc2 import rk4

G, a, T, dt = 1.0, 0.05, 4.088, 0.004
s = np.array([0.773723, -1e-9, 1.185226, 1e-9])
ds = []
for _ in range(int(T / dt)):
    s = rk4(s, dt)
    ds.append(np.sqrt((s[0] - s[2]) ** 2 + (s[1] - s[3]) ** 2) / a)
ds = np.array(ds)
print(f"min(d/a)={ds.min():.2f} mean={ds.mean():.2f} frac-below-4={(ds < 4).mean():.3f}")
