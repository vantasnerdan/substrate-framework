"""Exact complete-action jet and geometric-gradient-lift receipts."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0047-geometric-gradient-lift")
    h, mass_h, pairing = sp.symbols("h t B", positive=True)
    mixed, a, b2, parameter, frequency = sp.symbols("g a b2 epsilon nu", real=True)
    u, v, w, s11, s12, s22 = sp.symbols("u v w S11 S12 S22", real=True)
    h0 = sp.Matrix([[h, mixed], [mixed, mass_h]])
    h2 = sp.Matrix([[u, v], [v, w]])
    symplectic = sp.Matrix([[0, 1], [-1, 0]])
    symmetric = sp.Matrix([[s11, s12], [s12, s22]])
    energy = h0 + sp.I * parameter * a * symplectic + parameter**2 * h2
    omega = (
        pairing * symplectic
        + sp.I * parameter * symmetric
        + parameter**2 * b2 * symplectic
    )
    polynomial = sp.expand((energy - sp.I * frequency * omega).det())
    determinant0 = h0.det()
    gap2 = determinant0 / pairing**2
    linear = 2 * pairing * a + sp.trace(h0.adjugate() * symmetric)
    potential2 = sp.trace(h0.adjugate() * h2) - a**2
    inertia2 = symmetric.det() - 2 * pairing * b2

    def zero(name: str, expression: sp.Expr) -> None:
        ledger.check(name, sp.simplify(expression) == 0)

    zero(
        "complete characteristic constant term",
        polynomial.coeff(parameter, 0) - determinant0 + pairing**2 * frequency**2,
    )
    zero(
        "complete characteristic retains chiral first derivative",
        polynomial.coeff(parameter, 1) - frequency * linear,
    )
    zero(
        "complete characteristic quadratic term",
        polynomial.coeff(parameter, 2) - potential2 - frequency**2 * inertia2,
    )
    normalized_c = sp.simplify((potential2 + gap2 * inertia2) / mass_h)
    stiffness, correction11, correction12, correction_a = sp.symbols(
        "T r11 r12 da", real=True
    )
    changed_c = normalized_c.subs(
        {u: u + stiffness + correction11, v: v + correction12, a: a + correction_a},
        simultaneous=True,
    )
    expected_change = (
        stiffness
        + correction11
        - 2 * mixed * correction12 / mass_h
        - (2 * a * correction_a + correction_a**2) / mass_h
    )
    zero(
        "gradient-cage energy survives full shape and KKS elimination",
        changed_c - normalized_c - expected_change,
    )
    zero(
        "added q gradient enters with exact unit coefficient",
        sp.diff(normalized_c, u) - 1,
    )
    zero(
        "momentum-gradient contribution remains in the formula",
        sp.diff(normalized_c, w) - h / mass_h,
    )
    zero(
        "KKS-gradient contribution remains in the formula",
        sp.diff(normalized_c, b2) + 2 * gap2 * pairing / mass_h,
    )
    omega0 = pairing * symplectic
    omega1 = sp.I * symmetric
    darboux1 = -omega0.inv() * omega1 / 2
    ledger.check(
        "first-order Darboux transformation is derived",
        (darboux1.conjugate().T * omega0 + omega0 * darboux1 + omega1).applyfunc(
            sp.simplify
        )
        == sp.zeros(2),
    )

    macro, bond = sp.symbols("kappa d", real=True)
    difference = sp.exp(sp.I * bond * macro / 2) - sp.exp(-sp.I * bond * macro / 2)
    zero(
        "centered neighboring-angle difference vanishes at zero wave",
        difference.subs(macro, 0),
    )
    zero(
        "centered difference gives geometric first-gradient amplitude",
        sp.diff(difference, macro).subs(macro, 0) - sp.I * bond,
    )
    zero(
        "centered difference has no quadratic amplitude term",
        sp.diff(difference, macro, 2).subs(macro, 0),
    )
    zero(
        "its energy starts at the squared geometric bond length",
        sp.diff(sp.conjugate(difference) * difference, macro, 2).subs(macro, 0) / 2
        - bond**2,
    )
    self_kks1 = sp.symbols("self_kks1", real=True)
    self_added = sp.conjugate(difference) * difference * sp.I * macro * self_kks1
    zero(
        "self KKS correction starts beyond gradient order",
        sp.diff(self_added, macro, 2).subs(macro, 0),
    )
    rho, carrier, lam, principal, error, geom = sp.symbols(
        "rho K lambda A C d", positive=True
    )
    lower = rho * geom**2 * ((1 + carrier / lam) * principal - error)
    ledger.check(
        "actual-EPS carrier lower bound grows strictly",
        sp.diff(lower, carrier).is_positive,
    )
    ledger.check(
        "omitting the chiral energy term changes the action coefficient",
        sp.simplify(normalized_c - normalized_c.subs(a, 0)) != 0,
    )
    ledger.check(
        "omitting KKS curvature changes the action coefficient",
        sp.simplify(normalized_c - normalized_c.subs(b2, 0)) != 0,
    )
    print("C_lab(e) =", sp.factor(normalized_c))
    print("Exact lift delta C =", expected_change)
    print(
        "Scope: full local action jet. Compact-support estimates and finite positive thresholds are proved in gradient-cage-proof.md."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
