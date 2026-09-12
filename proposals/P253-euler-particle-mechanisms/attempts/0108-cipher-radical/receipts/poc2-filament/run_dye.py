"""PoC-2 filamentation-audit proxy: passive dye support diameter over one leapfrog period.

Standby idea → artifact (shepherd bounded task). S9 lesson: return-map/Floquet closure
is blind to support filamentation; freeze the pair (D(t), N) with N trivially conserved
for passive dye and D(t) the sensitive member.
Run: python3 run_dye.py from receipts/. Env: CPython 3.12.2, numpy 1.26.4.
Scope: REDUCED axisymmetric model (m>=1 invisible — stated cap, Euler check will see it).
"""
import time
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_poc2 import rk4, mutual

G, a = 1.0, 0.05
T = 4.088
DT = 0.004


def field(s, z, rings):
    v = np.zeros(2)
    for Rr, Zr in rings:
        m = mutual(s, z, Rr, Zr)
        v[0] += m[0]
        v[1] += m[2]
    return v


def main():
    t0 = time.time()
    rng = np.random.default_rng(7)
    # dye disks radius 3a around each ring cross-section (frozen seed geometry)
    seeds = []
    for Rr in (0.773723, 1.185226):
        r = 3 * a * np.sqrt(rng.random(256))
        ph = rng.random(256) * 2 * np.pi
        seeds.append(np.stack([np.full(256, Rr) + r * np.cos(ph), r * np.sin(ph)], axis=1))
    dye = np.concatenate(seeds, axis=0)
    s = np.array([0.773723, -1e-9, 1.185226, 1e-9])
    c0 = dye.mean(axis=0)
    D0 = np.abs(dye - c0).max() * 2

    def advect(d, rings, dt):
        k1 = np.array([field(p[0], p[1], rings) for p in d])
        k2 = np.array([field(*(p + 0.5 * dt * k), rings) for p, k in zip(d, k1)])
        return d + dt * k2

    n = int(T / DT)
    for _ in range(n):
        s = rk4(s, DT)
        rings = [(s[0], s[1]), (s[2], s[3])]
        dye = advect(dye, rings, DT)
    c1 = dye.mean(axis=0)
    D1 = np.abs(dye - c1).max() * 2
    print(f"D0={D0:.4f} D1={D1:.4f} D(T)/D(0)={D1/D0:.4f} N={len(dye)} (conserved)")
    print(f"centroid drift={(c1-c0).tolist()} (translation-invariant audit: D is centroid-subtracted)")
    print(f"elapsed {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
