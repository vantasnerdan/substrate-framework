#!/usr/bin/env python3
"""D-D2 edge pair in imposed tilt wave (beacon, 0157/pair, COMPUTE AUTHORIZED).
Frozen: pair/design.md. Banked inputs only: J(eps) (edge-T1'),
w = 1/2.J(eps).th^2 + W_FA + J0/ell^3 scale (t2-design.md),
D-D1 edge PK forms, E2 monopole kill, RC7 additivity, T2 F3
estimate-level Peierls (upgraded here to quadrature).
Mechanism (frozen intent): single-edge self-shift
U(x1)/L = (J0/2ell^3).int(4.eta_in-3.eps_zz)(x-x1).th(x)^2 d2x;
pair bound (adjacent troughs, d*=pi/q) iff 2.DU > F_rep(d*).d*
(Q3). Cross O(eta^2) OUT per STOP. Glide trap from diagonal
prestress is identically 0 (T2 F1 travels) — this self-shift
channel is the only derived-order trap; a verdict EITHER way
satisfies Q3. Exit nonzero on guard breach.
"""

import numpy as np

MU = 1.0
B0 = 1.0
A = 1.0
NU = -1.0
BETA = 6.822e-03   # (J0/ell^3)/mu, T2 receipt
TH = 0.1            # fiducial tilt amplitude
H = 60.0            # half-box (long-wave: box >> wavelength)
N = 600             # grid (dx = 0.2a; core masked)


def edge_fields(X, Y):
    """Volterra edge (b=xhat) strain, nu=-1. D-D1 forms, cited."""
    r2 = X ** 2 + Y ** 2
    fac = B0 / (2 * np.pi * (1 - NU))
    exx = -fac * Y * (3 * X ** 2 + Y ** 2) / r2 ** 2
    eyy = fac * Y * (X ** 2 - Y ** 2) / r2 ** 2
    ezz = np.zeros_like(X)  # plane strain
    eta_in = (exx + eyy) / 2.0  # -fac.Y/r^2: ODD in Y (parity null below)
    return eta_in, ezz


def trap_barrier(q, p, th, H, N, mask_core=2.0):
    """Barrier DU/L for one edge in oblique wave th.cos(qx+py).
    PARITY: M odd in Y -> normal wave (p=0) gives EXACT null; p != 0
    opens the dU/dy-moment channel (receipted mechanism). Units: J0/ell^3
    = BETA.MU. Scans x1 over one x-period pi/q at y1 = 0."""
    x = np.linspace(-H, H, N)
    dx = x[1] - x[0]
    X, Y = np.meshgrid(x, x)
    eta_in, ezz = edge_fields(X, Y)
    M = 4.0 * eta_in - 3.0 * ezz          # J-modulation pattern
    r = np.sqrt(X ** 2 + Y ** 2)
    ok = (r > mask_core * A) & (np.abs(X) < H - 5) & (np.abs(Y) < H - 5)
    M = np.where(ok, M, 0.0)
    period = np.pi / q
    xs = np.linspace(0, period, 25, endpoint=False)
    U = []
    for x1 in xs:
        th2 = th ** 2 * np.cos(q * (X + x1) + p * Y) ** 2
        U.append(0.5 * BETA * MU * np.sum(M * th2) * dx * dx)
    U = np.array(U)
    return U.max() - U.min()


def main() -> None:
    q = 2 * np.pi / 20.0  # fiducial wave (qa=0.31; T2 F3 caveat travels)
    p = q                 # 45-degree oblique: opens parity channel
    # Parity receipt: normal-incidence trap must be EXACT null
    DUp0 = trap_barrier(q, 0.0, TH, H, N)
    print(f"parity receipt: p=0 barrier = {DUp0:.3e} (exact null expected)")
    assert DUp0 < 1e-12, "PARITY FIRES: normal wave traps (mechanism wrong)"
    DU = trap_barrier(q, p, TH, H, N)
    # Guard G2: zero tilt -> zero trap
    assert trap_barrier(q, p, 0.0, H, N) == 0.0, "G2 FIRES: trap at zero tilt"
    # Guard G3: Th^2 scaling (structure of w = 1/2.J.th^2)
    assert abs(trap_barrier(q, p, 2 * TH, H, N) / DU - 4.0) < 1e-9, \
        "G3 FIRES: not th^2"
    # Guard G4: box stability (far-field oscillatory tail must not govern)
    DUbig = trap_barrier(q, p, TH, 1.5 * H, int(1.5 * N))
    print(f"D-D2 single-edge trap: DU/L = {DU:.3e}, box1.5x = {DUbig:.3e}")
    assert abs(DUbig - DU) / DU < 0.2, \
        "G4 FIRES: box-dominated barrier (disclose, hold lane)"
    print(f"G2/G3/G4 green; period pi/q = {np.pi / q:.2f}a")
    # Q3: pair at adjacent troughs d* = pi/q vs glide repulsion (D-D1)
    dstar = np.pi / q
    F_rep = 0.25 / (np.pi * dstar)  # D-D1 same-sign glide, mu=b=1
    work = F_rep * dstar            # q-independent: 0.25/pi
    print(f"Q3: 2.DU = {2 * DU:.3e} vs unbind work F_rep.d* = {work:.3e}")
    ratio = 2 * DU / work
    verdict = "BOUND" if ratio > 1 else "UNBOUND at fiducial"
    print(f"Q3 verdict: ratio = {ratio:.3e} -> {verdict}")
    print(f"parametric reopen: bound iff beta.th^2.(qa)^n.O(1) ~ 1 "
          f"(fiducial beta.th^2.qa = {BETA * TH**2 * q * A:.2e}); "
          "T2 F3 averaging-suppression caveat travels.")
    print(f"D-D2 verdict: tilt-wave pair {verdict} (oblique parity channel; "
          "normal-incidence null exact); inequality is the verdict.")


if __name__ == "__main__":
    main()
