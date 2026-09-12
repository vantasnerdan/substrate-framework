#!/usr/bin/env python3
"""P253/0122 contour-fitted mesh (beacon, attempt-local).

Extracts the src=0 contour from the trust member, then red-refines
(base 40x20 tensor) twice inside a band around it. Output: adapted
(node, connectivity) npz + report. No solver changes here.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np


def extract_contour(fu, mu, c, ctr=(1.0, 0.0), nang=72):
    pts = []
    for th in np.linspace(0, math.pi, nang):
        d = np.array([math.cos(th), math.sin(th)])

        def s(t):
            r, z = ctr[0] + t * d[0], ctr[1] + t * d[1]
            if r <= 1e-9 or r >= 6.0 or z < 0 or z >= 3.0:
                return None
            try:
                v = float(fu(np.array([[r], [z]]))[0])
            except ValueError:
                return None
            return v - c * r**2 / 2 - mu

        s0 = s(0.0)
        if s0 is None or s0 <= 0:
            continue
        lo, hi = 0.0, 0.6
        while True:
            h = s(hi)
            if h is None or hi >= 2.0:
                break
            if h <= 0:
                break
            hi *= 1.5
        if hi >= 2.0 or s(hi) is None:
            continue
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            sm = s(mid)
            if sm is None:
                break
            if sm > 0:
                lo = mid
            else:
                hi = mid
        else:
            pts.append([ctr[0] + lo * d[0], ctr[1] + lo * d[1]])
            continue
    return np.array(pts)


def dist_to_contour(p, cpts):
    d = np.sqrt(((p[:, None, :] - cpts[None, :, :]) ** 2).sum(axis=2))
    return d.min(axis=1)


def red_refine_band(p, t, cpts, width):
    """One red-refinement pass on triangles with centroid in band."""
    tri = p[t.T]
    cent = tri.mean(axis=1)
    mark = dist_to_contour(cent, cpts) < width
    pm = p.tolist()
    key = {}
    out_t = []

    def mid(a, b):
        k = (min(a, b), max(a, b))
        if k not in key:
            key[k] = len(pm)
            pm.append([round(float(v), 12)
                       for v in (p[a] + p[b]) / 2])
        return key[k]

    for j in range(t.shape[1]):
        i0, i1, i2 = (int(v) for v in t[:, j])
        if not mark[j]:
            out_t.append([i0, i1, i2])
            continue
        m01, m12, m20 = mid(i0, i1), mid(i1, i2), mid(i2, i0)
        out_t += [[i0, m01, m20], [i1, m12, m01], [i2, m20, m12],
                  [m01, m12, m20]]
    return np.array(pm), np.array(out_t).T


def main() -> None:
    import sys
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1
    att = "proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member"
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    d = np.load(f"{att}/member-trust-r3.npz")
    cpts = extract_contour(fu, mu, c)
    cpts = np.vstack([cpts, cpts * np.array([1.0, -1.0])])  # R3 mirror
    print(f"contour: {len(cpts)} points", flush=True)
    p = mesh.p.T.copy()
    t = mesh.t.copy()
    for level, width in ((1, 0.30), (2, 0.15)):
        p, t = red_refine_band(p, t, cpts, width)
        print(f"level {level}: nodes={p.shape[0]} elems={t.shape[1]}",
              flush=True)
    out = "proposals/P253-euler-particle-mechanisms/attempts/0122-beacon-fitted/"
    np.savez(f"{out}fitted-mesh.npz", p=p, t=t,
             mu=mu, c=c, kap=float(d["kap"]), rbar=float(d["rbar"]))
    print("saved fitted-mesh.npz", flush=True)


if __name__ == "__main__":
    main()
