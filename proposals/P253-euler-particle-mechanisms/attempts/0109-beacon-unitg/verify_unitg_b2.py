#!/usr/bin/env python3
"""P253/0109 Unit-G witness + B2 edge-transfer exposing checks (beacon).

Scope: finite-dimensional row algebra and frozen-column identities ONLY.
Imports canonical framework APIs; duplicates nothing. Does NOT prove cutoff
convergence, witness nondegeneracy on the actual Cao member, IFT hypotheses,
the accessible finite-Cao transfer, or any P2/particle claim. Those remain
pinned in witness-status.md / b2-edge-transfer.md.

New predicates beyond P253/0107's 10 (whose receipt is cited, not rerun):
  B1 column square, B2 closed-form exp instance, B3 edge shear, B4 sqrt(R/Z)
  amplification, B5 order-zero symmetrizer, B6 K-center value,
  G1 block lower-triangular determinant instance, G2 Cao-jet import smoke.
"""

from __future__ import annotations

import sympy as sp
from substrate_framework.euler_p2_principal import charged_column_slow_block
from substrate_framework.euler_cao_schur import cao_thin_ring_schur_jet


def check(label: str, predicate: bool) -> None:
    if not predicate:
        raise AssertionError(label)
    print(f"PASS {label}")


def main() -> None:
    x, Z, R, t = sp.symbols("x Z R t", real=True)
    xp, Zp, Rp = sp.symbols("xp Zp Rp", positive=True)

    # B1: C^2 = -x^2 Z R I (0107 derivation (47); framework block constructor).
    C = charged_column_slow_block(x, Z, R)
    check("B1 column square", sp.simplify(C**2 + x**2 * Z * R * sp.eye(2)) == sp.zeros(2))

    # B5: H0 = diag(R,Z) symmetrizes the frozen column: H0 C skew (51a).
    H0 = sp.diag(R, Z)
    HC = H0 * C
    check("B5 order-zero symmetrizer", sp.simplify(HC + HC.T) == sp.zeros(2))

    # B2: closed-form exp instance (51c) against SymPy matrix exponential.
    vals = {x: sp.Rational(3, 10), Z: sp.Rational(1, 2), R: sp.Rational(11, 10),
            t: sp.Rational(7, 10)}
    Cn = C.subs(vals)
    nu = abs(vals[x]) * sp.sqrt(vals[Z] * vals[R])
    closed = sp.cos(nu * vals[t]) * sp.eye(2) + sp.sin(nu * vals[t]) / nu * Cn
    exact = sp.Matrix((Cn * vals[t]).exp())  # Jordan oracle at same t, independent of closed form
    diff = sp.N((sp.Matrix(exact) - closed).norm(), 30)
    check("B2 frozen-column exp closed form", diff < sp.Float(10) ** -25)

    # B3: Z = 0 edge shear (51d): nilpotent square, exp = I + tC.
    C0 = charged_column_slow_block(x, 0, R)
    check("B3 edge nilpotent square", sp.simplify(C0**2) == sp.zeros(2))
    shr = sp.eye(2) + vals[t] * C0.subs({x: vals[x], R: vals[R]})
    check("B3 edge shear entry", sp.simplify(shr[1, 0] - vals[t] * vals[x] * vals[R]) == 0)

    # B4: at t_s = pi/(2 nu), lower entry magnitude = sqrt(R/Z) (R,Z > 0).
    Cp = charged_column_slow_block(xp, Zp, Rp)
    nup = xp * sp.sqrt(Zp * Rp)  # |x| = x for xp > 0
    ts = sp.pi / (2 * nup)
    lower = (sp.sin(nup * ts) / nup * Cp)[1, 0]
    check("B4 edge amplification sqrt(R/Z)",
          sp.simplify(lower - sp.sqrt(Rp / Zp)) == 0)

    # B6: smooth-center classifier K = 8 zeta/Omega -> 16 at zeta = 2 Omega.
    Om = sp.symbols("Omega", positive=True)
    check("B6 center K equals 16", sp.simplify(8 * (2 * Om) / Om) == 16)

    # G1: 9x9 staged block determinant instance: det = detH detC detG.
    rng = [[3, 1, 4], [1, 5, 9], [2, 6, 5]]
    H = sp.Matrix(rng)
    Cc = sp.Matrix([[2, 7, 1], [8, 2, 8], [1, 8, 2]])
    G = sp.Matrix([[4, 5, 9], [0, 4, 5], [2, 3, 5]])
    KI = sp.Matrix([[1, 0, 2], [0, 1, 1], [3, 1, 0]])
    KC = sp.Matrix([[0, 2, 1], [1, 0, 3], [2, 2, 1]])
    M = sp.Matrix(sp.BlockMatrix([[H, sp.zeros(3), sp.zeros(3)],
                                  [sp.zeros(3), Cc, sp.zeros(3)],
                                  [KI, KC, G]]))
    check("G1 staged block determinant",
          M.det() == H.det() * Cc.det() * G.det() != 0)

    # G2: Cao thin-ring jet import smoke pins the P1 carrier dependency.
    jet = cao_thin_ring_schur_jet(1, 10, 5, 1)
    check("G2 Cao jet import smoke",
          jet.physical_axial_impulse.is_positive
          and jet.translation_speed.is_positive
          and jet.physical_schur_determinant != 0)

    print("ALL 9 UNITG-B2 CHECKS PASSED")


if __name__ == "__main__":
    main()
