"""Exact shared-mean canonical momentum and physical field-map oracle."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0102-canonical-coefficient-join")
    rho, P, B, Hq, kappa = s.symbols("rho P B Hq kappa", positive=True)
    C, N, t = s.symbols("C N t", real=True)
    q, qd, ud, p, r, vp = s.symbols("q qd ud p r vp", real=True)
    sp, sm = s.symbols("sp sm", real=True)
    j = B**2 / P

    # Independent time-reversal reactions, with their COMMON macro momentum.
    plus = B * sp * qd - (P * sp**2 + 2 * N * sp * q + Hq * q**2) / 2
    minus = -B * sm * qd - (P * sm**2 + 2 * N * sm * q + Hq * q**2) / 2
    paired = s.expand(((plus + minus) / 2).subs({sp: p + r, sm: r - p}))
    expected = B * p * qd - (P * (p**2 + r**2) + 2 * N * r * q + Hq * q**2) / 2
    ledger.check("independent paired momenta separate kinetic and static reactions",
                 s.simplify(paired - expected) == 0)
    ledger.check("even reaction gives full potential Schur complement",
                 s.simplify(paired.subs(r, -N * q / P)
                            - (B * p * qd - P * p**2 / 2 - (Hq - N**2 / P) * q**2 / 2)) == 0)

    # Explicit canonical V sector before elimination. Keeping its positive
    # one-form Bp term but omitting the negative H cross would manufacture spin.
    before = rho * vp * ud + t * C * vp * qd + t * B * p * ud - rho * vp**2 / 2 - t * B * p * vp
    vsol = ud + t * (C * qd - B * p) / rho
    after = s.expand(before.subs(vp, vsol))
    expected_mean = rho * ud**2 / 2 + C * t * ud * qd + t**2 * (C * qd - B * p) ** 2 / (2 * rho)
    ledger.check("full common macro momentum stationary equation",
                 s.simplify(s.diff(before, vp).subs(vp, vsol)) == 0)
    ledger.check("L-driven affine pairing cancels in complete mean reduction",
                 s.simplify(after - expected_mean) == 0)
    ledger.check("surviving mixed row is generator moment not induced-velocity spin",
                 s.simplify(s.diff(after, ud, qd) - C * t) == 0)
    ledger.check("retained Bs Udot mutation is exposed",
                 s.simplify(after.coeff(ud, 1).coeff(p, 1)) == 0)

    lagrangian = after + B * p * qd - P * p**2 / 2 - kappa * q**2 / 2
    psol = B * (1 - t**2 * C / rho) * qd / (P - t**2 * B**2 / rho)
    reduced = s.factor(lagrangian.subs(p, psol))
    jeff = t**2 * C**2 / rho + B**2 * (1 - t**2 * C / rho) ** 2 / (P - t**2 * B**2 / rho)
    ledger.check("complete reaction equation retains mean square",
                 s.simplify(s.diff(lagrangian, p).subs(p, psol)) == 0)
    ledger.check("full gradient inertia before moment completion",
                 s.simplify(s.diff(reduced, qd, qd) - jeff) == 0)
    ledger.check("general small-gradient inertia correction is squared mismatch",
                 s.simplify(s.series(jeff, t, 0, 3).removeO() - j - t**2 * (C - j) ** 2 / rho) == 0)
    ledger.check("matched moment gives exact derived inertia not merely its jet",
                 s.simplify(jeff.subs(C, j) - j) == 0)
    ledger.check("matched moment leaves physical reaction unchanged",
                 s.simplify(psol.subs(C, j) - B * qd / P) == 0)
    ledger.check("matched moment makes the retained hybrid residual zero",
                 s.simplify((C * qd - B * psol).subs(C, j)) == 0)

    # Independent physical coordinate pullback, including the U-gradient mass.
    phid = s.symbols("phid", real=True)
    matched = s.simplify(reduced.subs(C, j).subs(qd, phid - t * ud))
    expected_kinetic = rho * ud**2 / 2 + j * phid**2 / 2 - j * t**2 * ud**2 / 2
    ledger.check("physical-angle kinetic pullback retains negative macro gradient mass",
                 s.simplify(matched - expected_kinetic + kappa * q**2 / 2) == 0)
    ledger.check("leading physical mixed mass is exactly zero after completion",
                 s.diff(matched, ud, phid) == 0)
    b = (C - j) / 2
    g = -kappa / 2
    ell = s.simplify(g - kappa * b / j)
    ledger.check("physical coupling is computed from generator moment",
                 s.simplify(ell + kappa * C / (2 * j)) == 0)
    ledger.check("completed pair has nonzero physical optical transfer",
                 s.simplify(ell.subs(C, j) + kappa / 2) == 0)
    ledger.check("zero generator moment mutation reproduces the old cancellation",
                 ell.subs(C, 0) == 0)

    # Six response directions, each with its own support: zero response KKS.
    Bc, amplitude, Pactual = s.symbols("Bc amplitude Pactual", nonzero=True)
    actual_B = amplitude * Bc
    actual_j = actual_B**2 / Pactual
    raw_q_moments = s.Matrix(s.symbols("mq:6"))
    raw_s_moments = s.Matrix(s.symbols("ms:6"))
    moment_matrix = s.Matrix.hstack(raw_q_moments, raw_s_moments, s.eye(6))
    target_q = s.Matrix([actual_j, 0, 0, 0, 0, 0])
    target_s = s.Matrix([0, 0, 0, actual_B, 0, 0])
    q_column = s.Matrix([1, 0, *(target_q - raw_q_moments)])
    s_column = s.Matrix([0, 1, *(target_s - raw_s_moments)])
    extended_kks = s.diag(s.Matrix([[0, actual_B], [-actual_B, 0]]), s.zeros(6))
    ledger.check("S has all six prescribed moments before computing j",
                 moment_matrix * s_column == target_s)
    ledger.check("Q generator moment is an explicit triangular six-row completion",
                 moment_matrix * q_column == target_q)
    ledger.check("six response corrections preserve the actual raw KKS pairing",
                 (q_column.T * extended_kks * s_column)[0] == actual_B)

    # Complete finite non-diagonal operator instance of the Woodbury identity.
    PP = s.Matrix([[3, 1], [1, 2]])
    DD = s.Matrix([[2, 1]])
    JJ = (DD * PP.inv() * DD.T)[0]
    pvector = PP.inv() * DD.T * qd
    residual = JJ * qd - (DD * pvector)[0]
    equation = PP * pvector - DD.T * qd - (t**2 / rho) * DD.T * residual
    ledger.check("full non-diagonal reaction operator has the same matched solution",
                 equation == s.zeros(2, 1) and residual == 0)

    kg, jg, tg = s.symbols("kg jg tg", positive=True)
    curvature, mp = s.symbols("curvature mp", real=True)
    normal = curvature - kappa * mp / j
    shifted = normal.subs({curvature: curvature + tg**2 * kg, mp: mp + tg**2 * jg})
    ledger.check("positive compact gradient keeps its entire inertia cost",
                 s.simplify(shifted - normal - tg**2 * (kg - kappa * jg / j)) == 0)
    crate, erate, a_rate = s.symbols("crate erate a_rate", real=True)
    rate_lagrangian = p * (B * qd + crate * erate + B * a_rate) - P * p**2 / 2
    rate_reaction = (B * qd + crate * erate + B * a_rate) / P
    rate_reduced = s.expand(rate_lagrangian.subs(p, rate_reaction))
    ledger.check("STF canonical rate source retains its complete reaction square",
                 s.simplify(rate_reduced - (B * qd + crate * erate + B * a_rate) ** 2 / (2 * P)) == 0)
    ledger.check("STF rate source diagonal inertia survives time-reversal pairing",
                 s.simplify(s.diff(rate_reduced, erate, erate) - crate**2 / P) == 0)
    strain_raw = s.Matrix(3, 3, s.symbols("strain:9"))
    strain = (strain_raw + strain_raw.T) / 2 - s.trace(strain_raw) * s.eye(3) / 3
    mixed_axial = s.Matrix([sum(s.LeviCivita(i, jj, kk) * strain[jj, kk]
                                   for jj in range(3) for kk in range(3)) for i in range(3)])
    ledger.check("isotropic STF-to-axial mixed source vanishes without deleting its norm",
                 mixed_axial == s.zeros(3, 1))
    geom, betad = s.symbols("geom betad", real=True)
    physical_spin = geom * betad + j * (phid - betad)
    ledger.check("actual affine spin is retained in the canonical observation map",
                 s.simplify(physical_spin - j * phid - (geom - j) * betad) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
