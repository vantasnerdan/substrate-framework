#!/usr/bin/env python3
"""T1 source (beacon, 0157/tilt, COMPUTE AUTHORIZED post-freeze 41fc3958).
Frozen design: tilt/design.md. Verdict logic frozen there; numeric
parameters below are implementation guards (mutation-style), the (a)
verdict rests on analytic exactness.
T1(a): J(eps)=J0(1+4.eta_in-3.eps_zz) on banked screw dipole:
  analytic: only eps_xz,eps_yz nonzero -> eta_in=eps_zz=0 EXACT.
  numeric FD cross-check guard: max|J/J0-1| < 1e-3 outside cores.
T1(b): O(eps^2) survey: tr(eps^2) profile of the dipole (geometric
  existence only; coupling constant UNPRICED — new formalism).
T1(c): edge-fallback SPEC (text footer of stdout, banked in receipts).
Units a=1,b0=1; dipole separation d=10 (dynamics d0); J0 scale-free (=1).
Exit nonzero on any guard breach.
"""

import numpy as np

B0 = 1.0
A = 1.0
D = 10.0
L = 40.0
N = 400
H = L / N


def theta(x, y, x0):
    return np.arctan2(y, x - x0)


def main() -> None:
    g = (np.arange(N) - N / 2) * H
    X, Y = np.meshgrid(g, g)
    # Screw dipole: +b0 at x=-d/2, -b0 at x=+d/2 (dipole-run convention)
    u = (B0 / (2 * np.pi)) * (theta(X, Y, -D / 2) - theta(X, Y, +D / 2))
    def dudx(x0, sgn):
        with np.errstate(invalid="ignore", divide="ignore"):
            return sgn * (B0 / (2 * np.pi)) * (-Y / ((X - x0) ** 2 + Y ** 2))

    def dudy(x0, sgn):
        with np.errstate(invalid="ignore", divide="ignore"):
            return sgn * (B0 / (2 * np.pi)) * ((X - x0) / ((X - x0) ** 2 + Y ** 2))
    exz_a = 0.5 * (dudx(-D / 2, 1.0) + dudx(D / 2, -1.0))
    eyz_a = 0.5 * (dudy(-D / 2, 1.0) + dudy(D / 2, -1.0))
    r1 = np.sqrt((X + D / 2) ** 2 + Y ** 2)
    r2 = np.sqrt((X - D / 2) ** 2 + Y ** 2)
    # Guard mask: smooth bulk (r>5a off cores AND off the arctan2
    # branch-cut band |y|<2H: single-valued FD crosses neither the
    # singularity nor the cut — banked persist cut-wall lesson; the
    # analytic verdict needs no grid at all). Survey keeps r>A.
    edge = np.zeros_like(X, bool)
    edge[:5, :] = edge[-5:, :] = edge[:, :5] = edge[:, -5:] = True
    smooth = (r1 > 5 * A) & (r2 > 5 * A) & (np.abs(Y) > 2 * H) & ~edge
    bulk = (r1 > A) & (r2 > A)
    print("T1(a) analytic: eta_in=0, eps_zz=0 EXACT "
          "(only eps_xz,eps_yz nonzero by construction)")
    # Numeric FD cross-check (guard, not verdict)
    ex_num = (np.roll(u, -1, 1) - np.roll(u, 1, 1)) / (2 * H)
    ey_num = (np.roll(u, -1, 0) - np.roll(u, 1, 0)) / (2 * H)
    exz_n, eyz_n = 0.5 * ex_num, 0.5 * ey_num
    dev = np.maximum(np.abs(exz_n - exz_a), np.abs(eyz_n - eyz_a))
    print("T1(a) numeric guard: max FD-vs-analytic strain dev (smooth) = "
          f"{np.max(dev[smooth]):.3e} (bar 1e-3)")
    assert np.max(dev[smooth]) < 1e-3, "T1(a) FIRES: FD disagrees (grid?)"
    print("T1(a) verdict: J(eps_screw)=J0 EXACT at derived order — "
          "screw SILENT (freeze-predicted null RECEIPTED)")
    # T1(b): tr(eps^2) = 2(exz^2+eyz^2), core-localized channel profile
    tr2 = 2 * (exz_a ** 2 + eyz_a ** 2)
    w = tr2 * H * H
    print(f"T1(b) survey: max tr(eps^2)={np.max(tr2[bulk]):.5f}, "
          f"integrated weight={np.sum(w[bulk]):.5f} (geometric profile only; "
          f"no energetics claimed)")
    assert np.max(tr2[bulk]) > 0, "T1(b) FIRES: no O(eps^2) channel at all"
    print("T1(b) verdict: O(eps^2) channel GEOMETRICALLY OPEN "
          "(core-localized ~1/r^4); coupling constant UNPRICED — "
          "deriving it is new formalism (STOP rule: needs charter)")
    print("T1(c) spec: edge defect (b in-plane) carries eta_in!=0 at O(eps) "
          "-> J(eps) FIRES at derived order. T1' build (Volterra edge in "
          "F-A medium, eta_in field, J modulation, tilt source) AWAITS "
          "CHARTER under its own mini-freeze. P1 fires FOR SCREW-IN-F-C.")


if __name__ == "__main__":
    main()
