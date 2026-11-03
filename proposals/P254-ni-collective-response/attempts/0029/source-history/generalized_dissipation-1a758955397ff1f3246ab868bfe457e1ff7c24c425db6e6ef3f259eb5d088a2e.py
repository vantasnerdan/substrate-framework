"""Exact generalized-coordinate power balance and declared dissipation.

An externally supplied scalar power constrains only the contraction of a
generalized force with the generalized velocity.  It does not normally select
the force components.  This module exposes that affine balance family and a
separate Rayleigh construction in which the damping matrix is declared.

Nothing here derives a radiated power, self-field, regularization, causal
response, source dynamics, radiation-reaction coefficient, or physical
generalized-coordinate metric.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, Iterable

import sympy as sp


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return sp.simplify(expression)


def _exact_nonnegative(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be explicitly nonnegative")
    return expression


def _exact_vector(values: Iterable[Any], name: str) -> sp.ImmutableMatrix:
    entries = tuple(_exact_real(value, f"{name}[{index}]") for index, value in enumerate(values))
    if not entries:
        raise ValueError(f"{name} must contain at least one component")
    return sp.ImmutableMatrix(entries)


def _exact_square_matrix(
    values: Any,
    size: int,
    name: str,
) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(values)
    if matrix.shape != (size, size):
        raise ValueError(f"{name} must have shape ({size}, {size})")
    entries = [
        _exact_real(matrix[row, column], f"{name}[{row},{column}]")
        for row in range(size)
        for column in range(size)
    ]
    return sp.ImmutableMatrix(size, size, entries)


def _require_symmetric(matrix: sp.ImmutableMatrix, name: str) -> None:
    if matrix != matrix.T:
        raise ValueError(f"{name} must be symmetric")


def _require_positive_definite(matrix: sp.ImmutableMatrix, name: str) -> None:
    _require_symmetric(matrix, name)
    for size in range(1, matrix.rows + 1):
        minor = sp.simplify(matrix[:size, :size].det())
        if minor.is_positive is not True:
            raise ValueError(f"{name} must be explicitly positive definite")


def _require_positive_semidefinite(matrix: sp.ImmutableMatrix, name: str) -> None:
    _require_symmetric(matrix, name)
    indices = range(matrix.rows)
    for size in range(1, matrix.rows + 1):
        for selected in combinations(indices, size):
            minor = sp.simplify(matrix.extract(selected, selected).det())
            if minor.is_nonnegative is not True:
                raise ValueError(f"{name} must be explicitly positive semidefinite")


@dataclass(frozen=True)
class MetricPowerBalance:
    """A declared-metric particular solution of ``force.rates=-power``.

    The returned force uniquely minimizes ``force.T*metric**-1*force`` among
    balanced forces.  Every other balanced force is this particular force plus
    a covector whose ordinary contraction with ``rates`` vanishes.  The metric
    is additional model data and may carry the coordinate-unit conversions.
    """

    power: sp.Expr
    rates: sp.ImmutableMatrix
    coordinate_metric: sp.ImmutableMatrix
    squared_rate: sp.Expr
    particular_force: sp.ImmutableMatrix
    balance_residual: sp.Expr


@dataclass(frozen=True)
class RayleighDissipation:
    """Exact ledger for a declared symmetric positive-semidefinite matrix."""

    rates: sp.ImmutableMatrix
    damping_matrix: sp.ImmutableMatrix
    rayleigh_function: sp.Expr
    generalized_force: sp.ImmutableMatrix
    dissipated_power: sp.Expr
    energy_rate_without_external_work: sp.Expr


def power_balance_residual(
    power: Any,
    rates: Iterable[Any],
    generalized_force: Iterable[Any],
) -> sp.Expr:
    """Return ``force.rates + power`` for exact real components."""

    supplied_power = _exact_nonnegative(power, "power")
    velocity = _exact_vector(rates, "rates")
    force = _exact_vector(generalized_force, "generalized_force")
    if force.rows != velocity.rows:
        raise ValueError("generalized_force and rates must have equal length")
    return sp.simplify((force.T * velocity)[0] + supplied_power)


def metric_power_balance(
    power: Any,
    rates: Iterable[Any],
    coordinate_metric: Any | None = None,
) -> MetricPowerBalance:
    """Return the minimum declared-metric force with work ``-power``.

    For nonzero rate vector ``u`` and positive-definite ``G``, the result is
    ``Q0=-power*G*u/(u.T*G*u)``.  It is only a particular allocation: the full
    family is ``Q0+z`` for arbitrary ``z`` satisfying ``z.T*u=0``.
    """

    supplied_power = _exact_nonnegative(power, "power")
    velocity = _exact_vector(rates, "rates")
    if coordinate_metric is None:
        metric = sp.ImmutableMatrix.eye(velocity.rows)
    else:
        metric = _exact_square_matrix(
            coordinate_metric,
            velocity.rows,
            "coordinate_metric",
        )
    _require_positive_definite(metric, "coordinate_metric")
    squared_rate = sp.simplify((velocity.T * metric * velocity)[0])
    if squared_rate.is_positive is not True:
        raise ValueError("rates must be explicitly nonzero in the coordinate metric")
    force = sp.ImmutableMatrix(
        [sp.simplify(value) for value in (-supplied_power * metric * velocity / squared_rate)]
    )
    residual = power_balance_residual(supplied_power, velocity, force)
    if residual != 0:
        raise AssertionError("constructed force does not satisfy power balance")
    return MetricPowerBalance(
        power=supplied_power,
        rates=velocity,
        coordinate_metric=metric,
        squared_rate=squared_rate,
        particular_force=force,
        balance_residual=residual,
    )


def scalar_power_balance_force(power: Any, rate: Any) -> sp.Expr:
    """Return the unique one-component balance force ``-power/rate``.

    The rate must be explicitly nonzero.  At zero rate, positive power is
    inconsistent with instantaneous work balance, while zero power leaves the
    force unconstrained; neither case has a unique quotient.
    """

    result = metric_power_balance(power, (rate,))
    return sp.simplify(result.particular_force[0])


def rayleigh_dissipation(
    rates: Iterable[Any],
    damping_matrix: Any,
) -> RayleighDissipation:
    """Construct ``R=u.T*D*u/2``, ``Q=-D*u``, and ``P=u.T*D*u``.

    ``D`` must be exact, symmetric, and explicitly positive semidefinite.  It
    is declared constitutive data; the function does not infer it from a power
    curve or from field radiation.
    """

    velocity = _exact_vector(rates, "rates")
    damping = _exact_square_matrix(
        damping_matrix,
        velocity.rows,
        "damping_matrix",
    )
    _require_positive_semidefinite(damping, "damping_matrix")
    power = sp.simplify((velocity.T * damping * velocity)[0])
    rayleigh = sp.simplify(power / 2)
    force = sp.ImmutableMatrix([sp.simplify(value) for value in (-damping * velocity)])
    energy_rate = sp.simplify((force.T * velocity)[0])
    if sp.simplify(energy_rate + power) != 0:
        raise AssertionError("Rayleigh force and dissipated power disagree")
    return RayleighDissipation(
        rates=velocity,
        damping_matrix=damping,
        rayleigh_function=rayleigh,
        generalized_force=force,
        dissipated_power=power,
        energy_rate_without_external_work=energy_rate,
    )


def energy_rate_with_external_work(
    dissipation: RayleighDissipation,
    external_force: Iterable[Any],
) -> sp.Expr:
    """Return ``external_force.rates-dissipated_power`` exactly."""

    supplied = _exact_vector(external_force, "external_force")
    if supplied.rows != dissipation.rates.rows:
        raise ValueError("external_force and rates must have equal length")
    return sp.simplify(
        (supplied.T * dissipation.rates)[0] - dissipation.dissipated_power
    )
