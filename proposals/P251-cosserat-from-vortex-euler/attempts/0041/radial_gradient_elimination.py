"""Exact exposing check for the conjugate-radius gradient correction."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0041-radius-gradient-elimination")
    h, stiffness, p0, cx, cq = sp.symbols("h K P0 Cx Cq", positive=True)
    k, x, q, qdot, psidot = sp.symbols("k x q qdot Psidot", real=True)
    action = (
        -p0 * x * qdot - (h + cx * k**2) * x**2 / 2 - (stiffness + cq * k**2) * q**2 / 2
    )

    def zero(name: str, expression: sp.Expr) -> None:
        ledger.check(name, sp.simplify(expression) == 0)

    solution = sp.solve(sp.diff(action, x), x)[0]
    eliminated = sp.simplify(action.subs(x, solution))
    zero(
        "radius momentum solved without deleting its gradient",
        solution + p0 * qdot / (h + cx * k**2),
    )
    zero(
        "exact eliminated inertia is positive rational operator",
        sp.diff(eliminated, qdot, 2) - p0**2 / (h + cx * k**2),
    )
    frequency2 = (h + cx * k**2) * (stiffness + cq * k**2) / p0**2
    zero(
        "optical k squared coefficient contains radius term",
        sp.diff(frequency2, k, 2).subs(k, 0) / 2 - (h * cq + stiffness * cx) / p0**2,
    )
    normalizer = sp.sqrt(1 + cx * k**2 / h)
    zero(
        "collective field normalizes exact inertia",
        (p0**2 / (h + cx * k**2)) * normalizer**2 - p0**2 / h,
    )
    normalized_stiffness = (stiffness + cq * k**2) * normalizer**2
    zero(
        "normalized local gradient coefficient",
        sp.diff(normalized_stiffness, k, 2).subs(k, 0) / 2 - cq - stiffness * cx / h,
    )
    zero(
        "normalized fourth derivative remainder retained",
        normalized_stiffness
        - stiffness
        - (cq + stiffness * cx / h) * k**2
        - cq * cx * k**4 / h,
    )
    ledger.check(
        "omitting radius gradient changes optical dispersion",
        sp.simplify((h * cq + stiffness * cx) / p0**2 - h * cq / p0**2) != 0,
    )
    zero(
        "hexagon singular coefficient gains derived factor ten",
        (cq + stiffness * cx / h).subs({stiffness: 9 * h, cx: cq}) - 10 * cq,
    )
    print(
        "Scope: free-relative sector; the affine-cage collective field requires its own full matrix map."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
