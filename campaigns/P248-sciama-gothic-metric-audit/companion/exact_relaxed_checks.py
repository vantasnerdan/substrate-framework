"""Exact relaxed-equation, density-weight, and weak-field checks for P248."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def run() -> int:
    ledger = CheckLedger("P248-RELAXED")
    t, x, y, z = sp.symbols("t x y z", real=True)
    coordinates = (t, x, y, z)
    potential = sp.Matrix(
        [
            x * t**2,
            -t * x**2 / 2,
            y * z**2,
            -z * y**2 / 2,
        ]
    )

    def box(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, x, 2) + sp.diff(expression, y, 2) + sp.diff(expression, z, 2) - sp.diff(expression, t, 2)

    divergence = sum(
        sp.diff(potential[index], coordinates[index]) for index in range(4)
    )
    divergence_of_box = sum(
        sp.diff(box(potential[index]), coordinates[index]) for index in range(4)
    )
    ledger.check(
        "commuting divergence with flat wave operator",
        sp.simplify(divergence_of_box - box(divergence)) == 0,
    )

    epsilon = sp.symbols("epsilon", real=True)
    h00, h11, h22, h33 = sp.symbols("h00 h11 h22 h33", real=True)
    eta = sp.diag(-1, 1, 1, 1)
    perturbation = sp.diag(h00, h11, h22, h33)
    metric = eta + epsilon * perturbation
    gothic = sp.sqrt(-metric.det()) * metric.inv()
    linear_gothic = gothic.applyfunc(
        lambda entry: sp.diff(entry, epsilon).subs(epsilon, 0)
    )
    raised_h = eta * perturbation * eta
    trace_h = sp.trace(eta * perturbation)
    expected = -raised_h + sp.Rational(1, 2) * trace_h * eta
    ledger.check(
        "weight-one gothic variable linearizes to minus trace reverse",
        (linear_gothic - expected).applyfunc(sp.simplify) == sp.zeros(4),
    )

    # A diagonal orientation-preserving coordinate rescaling supplies an
    # independent density-weight check.
    j1, j2, j3, j4 = sp.symbols("j1 j2 j3 j4", positive=True)
    jacobian = sp.diag(j1, j2, j3, j4)  # dx'/dx
    inverse_jacobian = jacobian.inv()
    base_metric = sp.diag(-2, 3, 5, 7)
    transformed_metric = inverse_jacobian.T * base_metric * inverse_jacobian
    transformed_gothic = sp.sqrt(-transformed_metric.det()) * transformed_metric.inv()
    base_gothic = sp.sqrt(-base_metric.det()) * base_metric.inv()
    expected_density = (
        jacobian * base_gothic * jacobian.T / jacobian.det()
    )
    ledger.check(
        "gothic inverse transforms as contravariant weight one",
        (transformed_gothic - expected_density).applyfunc(sp.simplify)
        == sp.zeros(4),
    )

    tetrad = sp.diag(2, 3, 5, 7)
    tetrad_metric = tetrad.T * eta * tetrad
    ledger.check(
        "tetrad determinant equals positive metric volume density",
        sp.sqrt(-tetrad_metric.det()) == tetrad.det(),
    )

    phi, speed = sp.symbols("Phi c", real=True, nonzero=True)
    newtonian_index = 1 - 2 * phi / speed**2
    ledger.check(
        "weak optical potential convention is explicit",
        sp.diff(newtonian_index, phi).subs(phi, 0) == -2 / speed**2,
    )

    # Order bookkeeping: theta starts at quadratic order, so first order sees
    # matter only and the next source contains the quadratic self-term once.
    order = sp.symbols("epsilon", real=True)
    matter0, matter1, self2 = sp.symbols("T0 T1 theta2", real=True)
    source = matter0 + order * matter1 + order**2 * self2
    ledger.check(
        "post-Minkowskian source orders are not duplicated",
        sp.expand(source).coeff(order, 0) == matter0
        and sp.expand(source).coeff(order, 1) == matter1
        and sp.expand(source).coeff(order, 2) == self2,
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(run())
