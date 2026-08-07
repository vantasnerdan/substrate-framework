#!/usr/bin/env python3
"""Independent algebraic rederivation of the P228 frame, gauge ray, and G."""

from __future__ import annotations

import sympy as sp

from substrate_framework.dimensional_sine_gordon import (
    dimensional_sine_gordon_coefficients,
)
from substrate_framework.emergent_spin2 import (
    sine_gordon_spin_two_coupling_ledger,
)
from substrate_framework.verification import CheckLedger


SPACETIME_COMPONENTS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 2),
    (2, 3),
    (3, 3),
)


def _direct_quadratic_gauge_constraint(momentum: tuple[int, ...]) -> sp.Matrix:
    """Construct the coefficient constraint without package kinetic helpers."""

    metric = sp.diag(-1, 1, 1, 1)
    variables = sp.symbols("r00 r01 r02 r03 r11 r12 r13 r22 r23 r33")
    tensor = sp.zeros(4)
    for variable, (row, column) in zip(
        variables,
        SPACETIME_COMPONENTS,
        strict=True,
    ):
        tensor[row, column] = variable
        tensor[column, row] = variable
    raised = metric * tensor * metric
    covector = sp.Matrix(momentum)
    vector = metric * covector
    norm = (covector.T * vector)[0]
    trace = sp.trace(metric * tensor)
    divergence = sp.Matrix(
        [
            sum(covector[row] * raised[row, column] for row in range(4))
            for column in range(4)
        ]
    )
    invariants = (
        norm
        * sum(
            tensor[row, column] * raised[row, column]
            for row in range(4)
            for column in range(4)
        ),
        (divergence.T * metric * divergence)[0],
        sum(divergence[column] * covector[column] * trace for column in range(4)),
        norm * trace**2,
    )
    gauge = sp.Matrix(
        [
            [
                (momentum[row] if column == gauge_index else 0)
                + (momentum[column] if row == gauge_index else 0)
                for gauge_index in range(4)
            ]
            for row, column in SPACETIME_COMPONENTS
        ]
    )
    residuals = []
    for invariant in invariants:
        residuals.append(sp.hessian(invariant, variables) * gauge)
    return sp.Matrix(
        [
            [residual[row, column] for residual in residuals]
            for row in range(10)
            for column in range(4)
        ]
    )


def main() -> int:
    checks = CheckLedger("P228-INDEPENDENT-ALGEBRA")
    golden_ratio = (1 + sp.sqrt(5)) / 2
    normalization = sp.sqrt(1 + golden_ratio**2)
    directions = tuple(
        sp.Matrix(vector) / normalization
        for vector in (
            (0, 1, golden_ratio),
            (0, -1, golden_ratio),
            (1, golden_ratio, 0),
            (-1, golden_ratio, 0),
            (golden_ratio, 0, 1),
            (golden_ratio, 0, -1),
        )
    )
    projection = sp.Matrix(
        [
            [
                n[0] ** 2,
                n[1] ** 2,
                n[2] ** 2,
                2 * n[0] * n[1],
                2 * n[0] * n[2],
                2 * n[1] * n[2],
            ]
            for n in directions
        ]
    )
    hxx, hyy, hzz, hxy, hxz, hyz = sp.symbols("hxx hyy hzz hxy hxz hyz")
    tensor = sp.Matrix([[hxx, hxy, hxz], [hxy, hyy, hyz], [hxz, hyz, hzz]])
    components = sp.Matrix([hxx, hyy, hzz, hxy, hxz, hyz])
    channel_norm = sp.expand((projection * components).dot(projection * components))
    expected_norm = sp.Rational(2, 5) * sp.trace(tensor) ** 2 + sp.Rational(
        4, 5
    ) * sp.trace(tensor * tensor)
    checks.check(
        "direct icosahedral algebra reproduces the rank and tight-frame law",
        projection.rank() == 6
        and sp.simplify(channel_norm - expected_norm) == 0
        and sp.simplify(
            sum((n * n.T for n in directions), sp.zeros(3)) - 2 * sp.eye(3)
        )
        == sp.zeros(3),
    )

    ell, coordinate = sp.symbols("ell s", positive=True, real=True)
    kink_gradient_squared = 4 / (ell**2 * sp.cosh(coordinate / ell) ** 2)
    checks.check(
        "direct kink integration gives the wall translation normalization",
        sp.integrate(
            kink_gradient_squared,
            (coordinate, -sp.oo, sp.oo),
        )
        == 8 / ell,
    )

    coefficient_constraints = sp.Matrix.vstack(
        *(
            _direct_quadratic_gauge_constraint(momentum)
            for momentum in ((1, 2, 3, 5), (2, -1, 1, 3))
        )
    )
    nullspace = coefficient_constraints.nullspace()
    ray = sp.simplify(nullspace[0] / nullspace[0][0])
    checks.check(
        "an independent invariant Hessian solve reproduces the unique gauge ray",
        coefficient_constraints.rank() == 3
        and len(nullspace) == 1
        and tuple(ray) == (1, -2, 2, -1),
    )

    kx, ky, kz = 1, 2, 3
    constraints = sp.Matrix(
        [
            [kx, 0, 0, ky, kz, 0],
            [0, ky, 0, kx, 0, kz],
            [0, 0, kz, 0, kx, ky],
            [1, 1, 1, 0, 0, 0],
        ]
    )
    basis = sp.Matrix.hstack(*constraints.nullspace())
    projected_norm = sp.simplify(basis.T * sp.diag(1, 1, 1, 2, 2, 2) * basis)
    checks.check(
        "direct admissibility algebra leaves two positive tensor modes",
        constraints.rank() == 4
        and len(constraints.nullspace()) == 2
        and projected_norm.is_positive_definite is True,
    )

    lam, tension, onsite = sp.symbols("lambda T mu", positive=True)
    speed_squared = tension / lam
    direct_spacetime_normalization = sp.Rational(4, 5) * 8 * onsite
    direct_kappa = sp.simplify(1 / (4 * direct_spacetime_normalization))
    direct_newton = sp.simplify(speed_squared**2 * direct_kappa / (8 * sp.pi))
    package = sine_gordon_spin_two_coupling_ledger(
        dimensional_sine_gordon_coefficients(lam, tension, onsite)
    )
    checks.check(
        "independent normalization arithmetic reproduces package kappa and G",
        direct_spacetime_normalization == 32 * onsite / 5
        and direct_kappa == 5 / (128 * onsite)
        and direct_newton == 5 * tension**2 / (1024 * sp.pi * lam**2 * onsite)
        and package.einstein_coupling == direct_kappa
        and package.newton_constant == direct_newton,
    )
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
