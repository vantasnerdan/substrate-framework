#!/usr/bin/env python3
"""Edge-T1' compute (beacon, 0157/tilt, chartered post-mini-freeze).
Frozen design: tilt/edge-t1prime-design.md. Verdicts per E0-E3.
Single Volterra edge (line z, Burgers B0.xhat), banked nu=-1
(DERIVED from banked F-A constants in the freeze; kappa=0 footnote
travels). Observables: eta_in analytic + FD guard (E0, bar 1e-3);
E1 nonzero check; E2 angular-mean scope verdict (expect 0);
dJ/J0 = 4.eta_in modulation stats (peak, pattern); E3 core-shift
STAGED (recorded, not fired). Masks: cut band |y|<2H (arctan cut),
core r>A excluded (bulk) / r>5A (smooth guard), 5-cell wrap edge.
Units a=1, b0=1. Exit nonzero on E0/E1 breach.
"""

import numpy as np

B0 = 1.0
A = 1.0
NU = -1.0  # derived in freeze from banked mu=K/10, lam=-K/15
L = 40.0
N = 400
H = L / N


def main() -> None:
    g = (np.arange(N) - N / 2) * H
    X, Y = np.meshgrid(g, g)
    r2 = X ** 2 + Y ** 2
    th = np.arctan2(Y, X)
    fac = B0 / (2 * np.pi)
    ux = fac * (th + X * Y / (2 * (1 - NU) * r2))
    uy = -fac * ((1 - 2 * NU) * np.log(np.sqrt(r2)) / (2 * (1 - NU))
                 + (X ** 2 - Y ** 2) / (4 * (1 - NU) * r2))
    r = np.sqrt(r2)
    # Analytic eta_in (dipolar, zero mean)
    eta_a = -B0 * (1 - 2 * NU) * Y / (4 * np.pi * (1 - NU) * r2)
    edge = np.zeros_like(X, bool)
    edge[:5, :] = edge[-5:, :] = edge[:, :5] = edge[:, -5:] = True
    smooth = (r > 5 * A) & (np.abs(Y) > 2 * H) & ~edge
    bulk = (r > A) & ~edge
    # E0: FD guard on exx+eyy vs analytic 2.eta
    with np.errstate(invalid="ignore", divide="ignore"):
        exx = (np.roll(ux, -1, 1) - np.roll(ux, 1, 1)) / (2 * H)
        eyy = (np.roll(uy, -1, 0) - np.roll(uy, 1, 0)) / (2 * H)
    dev = np.abs((exx + eyy) / 2 - eta_a)
    print(f"E0 guard: max FD-vs-analytic eta_in dev (smooth) = "
          f"{np.nanmax(dev[smooth]):.3e} (bar 1e-3)")
    assert np.nanmax(dev[smooth]) < 1e-3, "E0 FIRES: grid/cut bug"
    # E1: dilatation genuinely present (implementation-grade existence)
    peak = np.nanmax(np.abs(eta_a[bulk]))
    print(f"E1: peak |eta_in| (bulk) = {peak:.5f} (must be >> discretization)")
    assert peak > 1e-3, "E1 FIRES: eta_in absent (implementation kill)"
    # E2: angular-mean scope verdict on ring r=8 (expect 0: no monopole)
    ring = np.abs(r - 8.0) < H
    m = float(np.mean(eta_a[ring]))
    print(f"E2: angular mean eta_in at r=8: {m:.3e} "
          f"(scale {peak:.3e}; expect ~0 -> no far-field monopole)")
    assert abs(m) / peak < 0.05, \
        "E2 SURPRISE: nonzero monopole (disclose, contradicts dipole form)"
    print("E2 verdict: far-field tilt-monopole claim KILLED IN SCOPE; "
          "near-field modulation stands")
    # Modulation map: dJ/J0 = 4.eta_in (eps_zz=0 plane strain)
    dJ = 4 * eta_a
    print(f" modulation: peak |dJ/J0|={np.nanmax(np.abs(dJ[bulk])):.5f} "
          f"at core flank; dipolar sin(th)/r, zero mean")
    # E3 staged: core-shift scale ~ b (Burgers scale) vs B1 bars — recorded
    print("EDGE-T1' verdict: edge FIRES J at derived order (P1-redirect "
          "SUCCESS); near-field tilt-modulation staged for gate-3.")


if __name__ == "__main__":
    main()
