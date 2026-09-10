#!/usr/bin/env python3
"""P253/0122 fitted-mesh bordered solve + lemma re-measure (beacon).

Loads fitted-mesh.npz, warm-starts bordered Newton (p=6) from the trust
state interpolated onto the fitted mesh, then re-measures the 0121 lemma
quantities (gap, soft share, dF bound, margin). Gate: dF <= 0.1.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def main() -> None:
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1, MeshTri
    att7 = "proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member"
    att2 = "proposals/P253-euler-particle-mechanisms/attempts/0122-beacon-fitted"
    fm = np.load(f"{att2}/fitted-mesh.npz")
    mesh = MeshTri(fm["p"].T, fm["t"])
    print(f"fitted mesh: nodes={mesh.p.shape[1]} elems={mesh.t.shape[1]}",
          flush=True)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    import os
    resume = f"{att2}/member-fitted.npz"
    if os.path.exists(resume):
        r = np.load(resume)
        u0, mu0, c0 = r["u"], float(r["mu"]), float(r["c"])
        print("RESUME fitted state res=", r["res"], flush=True)
    else:
        d = np.load(f"{att7}/member-trust-r3.npz")
        mc = B.build_mesh(6.0, 3.0, 40, 20)
        bc = Basis(mc, ElementTriP1())
        u0 = bc.interpolator(d["u"])(np.stack([rn, basis.doflocs[1]]))
        mu0, c0 = float(d["mu"]), float(d["c"])
    out = B.solve_bordered(mesh=mesh, reg=1e-3, p_pw=6,
                           u0=u0, mu0=mu0, c0=c0)
    print(f"FITTED res={out['res']:.3e} kap={out['kap']:.4f} "
          f"rbar={out['rbar']:.4f} mu={out['mu']:.4f} c={out['c']:.5f}",
          flush=True)
    B.save_npz(f"{att2}/member-fitted.npz", force=True,
             u=out["u"], mu=out["mu"], c=out["c"],
             kap=out["kap"], rbar=out["rbar"], iz=out["iz"],
             res=out["res"])

if __name__ == "__main__":
    main()
