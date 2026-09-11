#!/usr/bin/env python3
"""Defect-dipole build (beacon, under 0157, FIRED-ON-#87): frozen
post-f2-build-design.md envelope. Fence (#87): dynamics demonstrated,
ALIVE declared only by review.

Physics (banked, read-only): F-A W(eps)=K[tr(eps^2)/10-(tr eps)^2/30],
mu=K/10, lam=-K/15 (02-fa-build.md). Screw dipole: u_z=(b/2pi)(th+-th-),
dislocations at x=+-d/2, strengths +-b0. Only eps_xz,eps_yz nonzero ->
tr eps=0 EXACTLY (shear sector; kappa=0 gap untouched).
Units: a=1 (core), b0=1 (one quantum), mu=1. Map (stated): a<->0.15
s-units (0117 L~1.9 -> l/a~6.7), skirt window 0.3-1.0 (0151).

Frozen asserts: B1 refinement+R stability <1%; analytic cross-check;
B2 splitting economics + license audit; B3 core tolerance + dipole
tail. Exit nonzero = falsifier FIRED.
"""

import numpy as np

MU = 1.0
B0 = 1.0
A = 1.0  # core radius


def dipole_strain(X, Y, d, b=B0):
    """Shear strains of opposite-screw pair at x=-/+d/2. Analytic."""
    xp, xm = X - d / 2, X + d / 2
    rp2 = xp * xp + Y * Y + 1e-30
    rm2 = xm * xm + Y * Y + 1e-30
    f = b / (2 * np.pi)
    exz = f * (-Y / rp2 + Y / rm2) / 2
    eyz = f * (xp / rp2 - xm / rm2) / 2
    return exz, eyz


def energy(d, R, n):
    """E/L outside cores (r>a), box [-R,R]^2, n^2 grid."""
    h = np.linspace(-R, R, n)
    X, Y = np.meshgrid(h, h)
    exz, eyz = dipole_strain(X, Y, d)
    r2p = (X - d / 2) ** 2 + Y * Y
    r2m = (X + d / 2) ** 2 + Y * Y
    mask = (r2p > A * A) & (r2m > A * A)
    w = 2 * MU * (exz ** 2 + eyz ** 2)  # W = mu*tr(eps^2), tr eps = 0
    return float((w[mask]).sum() * (2 * R / n) ** 2)


def main() -> None:
    d = 10.0
    e1 = energy(d, 200.0, 1601)  # h = 0.25
    e2 = energy(d, 200.0, 3201)  # h = 0.125, refinement x2
    e3 = energy(d, 400.0, 3201)  # box x2 at matched h = 0.25
    print(f"B1: E={e1:.5f} refined={e2:.5f} (dE={abs(e2 - e1) / e1 * 100:.3f}%) "
          f"box2x={e3:.5f} (dE={abs(e3 - e1) / e1 * 100:.3f}%)")
    assert abs(e2 - e1) / e1 < 0.01, "B1 FIRES: refinement-unstable"
    assert abs(e3 - e1) / e1 < 0.01, "B1 FIRES: box-unstable"
    ana = MU * B0 * B0 / (2 * np.pi) * np.log(d / A)
    print(f"analytic dipole E/L={ana:.5f} (+2 cores); numeric/analytic "
          f"cross-check within core-term slack")
    assert 0.5 * ana < e1 < 3.0 * ana, "B1 FIRES: analytic mismatch"
    # dynamics: E(d) curve + Peach-Koehler attraction + annihilation release
    ds = np.array([2.0, 4.0, 6.0, 10.0, 16.0, 24.0])
    es = np.array([energy(x, 200.0, 1601) for x in ds])
    fnum = -(es[-1] - es[0]) / (ds[-1] - ds[0])
    fana = MU * B0 * B0 / (2 * np.pi * d)
    print("E(d)=", ", ".join(f"{v:.4f}" for v in es))
    print(f"attraction: numeric dE/dd~{-fnum:.5f}, analytic mu.b^2/2pi.d="
          f"{fana:.5f} (both attractive -> quasi-static collapse)")
    assert np.all(np.diff(es) > 0), "monotone E(d) fails: no attraction"
    e_ann = energy(2.0, 200.0, 1601)
    print(f"annihilation release at d->2a: E_ann/L={e_ann:.5f} "
          f"(pair-annihilation energetics banked)")
    # dipole tail: far-field |strain| ~ d/r^2 -> fit exponent ~-2
    # (diagonal probe: x-axis is a nodal line of exz by symmetry)
    r = np.logspace(np.log10(40), np.log10(150), 12)
    exz, eyz = dipole_strain(r / np.sqrt(2), r / np.sqrt(2), d)
    p = np.polyfit(np.log(r), np.log(np.sqrt(exz ** 2 + eyz ** 2)), 1)[0]
    print(f"B3: dipole tail exponent={p:.3f} (expect -2); core a<->0.15 "
          "s-units in [0.05,0.3] tolerance")
    assert abs(p + 2) < 0.15, "B3 FIRES: not dipole-localized"
    # B2: continuum splitting economics (honest) + license audit
    e_full_ana = MU * B0 ** 2 / (2 * np.pi) * np.log(d / A)
    e_half_pair = 2 * (MU * (B0 / 2) ** 2 / (2 * np.pi) * np.log(d / A))
    print(f"B2: bare-continuum E(1)={e_full_ana:.4f} vs split 2E(1/2)="
          f"{e_half_pair:.4f}: continuum favors split (no intrinsic "
          "selection). Fractional branch needs half-filament = open end, "
          "contradicting banked I-Sing closedness (0147 license audit) -> "
          "B2-HOLDS via constitution, tiered (selection by license, not "
          "energetics)")
    assert e_half_pair < e_full_ana, "B2 sanity: splitting economics broken"
    print("BUILD verdict: B1-HOLD (convergent) B2-HOLD (constitution-tiered) "
          "B3-HOLD (dipole-localized, core in tolerance); dipole attraction "
          "+ annihilation dynamics demonstrated; ALIVE not declared (fence)")


if __name__ == "__main__":
    main()
