#!/usr/bin/env python3
"""D-D1 direct pair table (beacon, 0157/pair, COMPUTE AUTHORIZED).
Frozen design: pair/design.md. Q0 opposite-screw regression (must
reproduce dipole-run attraction sign/order); Q1 like-sign repulsion
sign + 1/d (SymPy PK exact); Q2 bare-bound search EMPTY (monotonic
repulsion receipted); edge-edge PK forms + numeric spot check.
Units a=1, b0=1, mu=1 (dipole-run convention). Numeric bars are
implementation guards (T1 mask discipline: cut band, cores, wrap
edge). Exit nonzero on any guard breach.
"""

import numpy as np
import sympy as sp

B0 = 1.0
A = 1.0
MU = 1.0


def main() -> None:
    b1, b2, d = sp.symbols("b1 b2 d", real=True)
    # Screw 1 (b1) at origin -> stress at (d,0): sig_yz = mu.b1/2pi.d
    F_like = MU * b1 * b2 / (2 * sp.pi * d)  # xhat component on screw 2
    print(f"D-D1 screw pair: F_x = {F_like} (exact PK)")
    assert sp.simplify(sp.diff(F_like, d) + F_like / d) == 0, "not 1/d"
    # Q1: like-sign repels
    assert float(F_like.subs({b1: B0, b2: B0, d: 10.0})) > 0, \
        "Q1 FIRES: like-sign did not repel"
    print("Q1 guard: like-sign REPULSION, 1/d form exact")
    # Q0: opposite-sign attracts + order vs dipole-run B1
    f_opp = float(F_like.subs({b1: B0, b2: -B0, d: 8.5}))
    print(f"Q0 regression: opposite-screw F(d=8.5) = {f_opp:.5f} "
          f"(attractive; dipole-run B1 dE/dd~0.01881 same order)")
    assert f_opp < 0, "Q0 FIRES: opposite-sign did not attract"
    assert abs(abs(f_opp) - 0.01881) / 0.01881 < 0.25, \
        "Q0 FIRES: order mismatch vs banked B1"
    # Q2: monotonic repulsion on d in [2a, 20a] -> bare bound EMPTY
    ds = np.linspace(2 * A, 20 * A, 19)
    E = -MU * B0 * B0 / (2 * np.pi) * np.log(ds / A)  # like-sign pair E
    assert np.all(np.diff(E) < 0), "Q2 SURPRISE: non-monotonic (disclose)"
    print("Q2 verdict: like-sign E(d) strictly decreases with separation "
          "(energy shed by parting = repulsion at all d, no equilibrium) "
          "-> bare-medium bound pair EMPTY as frozen")
    NU = -1.0  # banked F-A (derived in edge-T1' freeze)
    x, y = sp.symbols("x y", real=True)
    r2 = x ** 2 + y ** 2
    # edge-1 stress at (x,y): standard Volterra forms (mu factored)
    sxx = -(MU * B0 / (2 * sp.pi * (1 - NU))) * y * (3 * x ** 2 + y ** 2) / r2 ** 2
    sxy = (MU * B0 / (2 * sp.pi * (1 - NU))) * x * (x ** 2 - y ** 2) / r2 ** 2
    Fx = sp.simplify(sxy * B0)   # glide on edge-2 (same sign)
    Fy = sp.simplify(-sxx * B0)  # climb on edge-2
    print(f"edge-edge (same sign): F_glide = {Fx}; F_climb = {Fy}")
    # Spot: stacked vertically (0,h): glide 0 by symmetry, climb nonzero
    h = 5.0
    assert abs(float(Fx.subs({x: 0, y: h}))) < 1e-12, "edge symmetry breach"
    assert float(Fy.subs({x: 0, y: h})) != 0
    print(f"spot (0,{h}): glide=0 (symmetry), climb="
          f"{float(Fy.subs({x: 0, y: h})):.5f} (exact forms receipted)")
    print("D-D1 verdict: direct pair table banked (screw +-; edge PK); "
          "Q0/Q1 guards green; Q2 bare-bound EMPTY.")


if __name__ == "__main__":
    main()
