"""Constrained-field normalization of the full point-limit gradient matrices."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0041-cage-gradient-matrix")
    rho, gamma, d, line = sp.symbols("rho Gamma d Tline", positive=True)
    inertia = 6 * sp.pi * rho * d**4
    stiffness = 27 * rho * gamma**2 / (8 * sp.pi)
    m0 = inertia * sp.Matrix(
        [
            [sp.Rational(11, 10), -sp.Rational(9, 10)],
            [-sp.Rational(9, 10), sp.Rational(11, 10)],
        ]
    )
    transform = sp.Matrix([[sp.Rational(2, 11), sp.Rational(9, 11)], [0, 1]])
    cp = line * sp.eye(2) / (3 * rho**2 * gamma**2 * d**2)
    ctheta = 3 * line * d**2 * sp.eye(2)
    ktheta = stiffness * sp.Matrix([[1, -1], [-1, 1]])
    kinetic = transform.T * m0 * transform
    gradient_inertia = transform.T * m0 * cp * m0 * transform
    angle_gradient = transform.T * ctheta * transform
    potential = transform.T * ktheta * transform

    def zero(name: str, expression: sp.Expr) -> None:
        ledger.check(name, sp.simplify(expression) == 0)

    zero("constant map diagonalizes kinetic term", kinetic[0, 1])
    zero("collective kinetic coefficient", kinetic[0, 0] - 2 * inertia / 55)
    zero("collective gap coefficient", potential[0, 0] - 4 * stiffness / 121)
    zero("twist-only leading gradient", angle_gradient[0, 0] - 12 * line * d**2 / 121)
    complete = sp.simplify(
        angle_gradient[0, 0] + potential[0, 0] * gradient_inertia[0, 0] / kinetic[0, 0]
    )
    zero(
        "complete constrained singular ratio",
        complete / angle_gradient[0, 0] - sp.Rational(1129, 220),
    )
    ledger.check(
        "radius-gradient inertia is nonzero",
        sp.simplify(gradient_inertia[0, 0]).is_positive,
    )
    ledger.check(
        "free-relative tenfold factor is not cage coefficient",
        sp.simplify(complete / angle_gradient[0, 0] - 10) != 0,
    )
    print(
        "Scope: derived leading-log point-limit matrix; finite-core coefficients remain exact integral expressions."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
