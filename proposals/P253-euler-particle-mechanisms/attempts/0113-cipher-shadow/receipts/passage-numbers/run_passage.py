"""Passage diagnostics archive for 0113 §3 (numbers only, no transfer claim).

Run: python3 run_passage.py from receipts/. Imports run_poc2 from the 0108 receipt
(same bytes drift reproduced). Env: CPython 3.12.2, numpy 1.26.4.
"""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..',
                                '0108-cipher-radical', 'receipts', 'poc2-filament'))
from run_poc2 import rk4

G, a = 1.0, 0.05
s = np.array([0.773723, -1e-9, 1.185226, 1e-9])
dt = 0.004
dmin, strainmax, Rmin = 1e9, 0.0, 1e9
for _ in range(int(4.088 / dt)):
    s = rk4(s, dt)
    R1, Z1, R2, Z2 = s
    d = np.sqrt((R1 - R2) ** 2 + (Z1 - Z2) ** 2)
    dmin = min(dmin, d)
    Rmin = min(Rmin, R1, R2)
    strainmax = max(strainmax, G / max(d * d, 1e-12))
Tcore = 2 * np.pi ** 2 * a * a / G
print(f"d_min={dmin:.4f} d/a={dmin/a:.1f} Rmin={Rmin:.4f} strainmax={strainmax:.3f}")
print(f"T_leap=4.088 T_core={Tcore:.4f} ratio={4.088/Tcore:.1f}")
print(f"strain/core-vort={strainmax/(G/(np.pi*a*a)):.3f}")
