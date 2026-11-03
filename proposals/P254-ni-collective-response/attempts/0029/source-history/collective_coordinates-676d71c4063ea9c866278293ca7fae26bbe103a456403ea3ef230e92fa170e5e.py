"""Exact conditional pullbacks to one collective coordinate.

The helpers in this module derive a reduced kinetic metric from a declared
field-profile family.  They do not select that family, prove its physical
realization, identify coefficients belonging to separately declared models,
or turn an unstable stationary point into a stable mode or event threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive_quantity(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_positive is False or (
        expression.is_number and expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be positive")
    return expression


def _real_quantity(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is False or (
        expression.is_number and expression.is_real is not True
    ):
        raise ValueError(f"{name} must be real")
    return expression


@dataclass(frozen=True)
class CollectiveCoordinateMetric:
    """Pulled-back kinetic metric for a declared field-profile family."""

    profile_derivative: sp.Expr
    geometric_integral: sp.Expr
    inertia: sp.Expr


@dataclass(frozen=True)
class CollectiveCoordinateLinearization:
    """Rest-stationary linearization of a one-coordinate Lagrangian."""

    inertia: sp.Expr
    curvature: sp.Expr
    spectral_ratio: sp.Expr
    characteristic_root_squared: sp.Expr
    characteristic_roots: tuple[sp.Expr, sp.Expr]
    stability: str
    stable_angular_frequency: sp.Expr | None
    instability_rate: sp.Expr | None


@dataclass(frozen=True)
class CollectiveCoordinateReparameterization:
    """Exact metric and potential transformation under ``q=g(Q)``."""

    inverse_map: sp.Expr
    jacobian: sp.Expr
    second_jacobian: sp.Expr
    transformed_inertia: sp.Expr
    transformed_potential: sp.Expr
    transformed_hessian: sp.Expr
    hessian_chain_term: sp.Expr
    stationary_hessian: sp.Expr
    stationary_spectral_ratio: sp.Expr


@dataclass(frozen=True)
class CollectiveCoordinateDimensionLedger:
    """Energy/length/time exponents for the one-coordinate pullback."""

    base_dimensions: tuple[str, str, str]
    quantity_names: tuple[str, ...]
    dimension_matrix: sp.ImmutableMatrix


def collective_coordinate_metric(
    inertia_density: Any,
    profile: Any,
    spatial_coordinate: sp.Symbol,
    collective_coordinate: sp.Symbol,
    *,
    lower_bound: Any = -sp.oo,
    upper_bound: Any = sp.oo,
) -> CollectiveCoordinateMetric:
    r"""Return ``lambda*integral((partial_q phi)**2 dx)``.

    The caller supplies the profile, domain, and coordinate convention.  The
    returned integral can be zero, infinite, or unevaluated; downstream claims
    of a positive finite inertia must check those properties explicitly.
    """

    coefficient = _positive_quantity(inertia_density, "inertia_density")
    field_profile = _real_quantity(profile, "profile")
    derivative = sp.simplify(sp.diff(field_profile, collective_coordinate))
    geometric_integral = sp.integrate(
        derivative**2,
        (spatial_coordinate, sp.sympify(lower_bound), sp.sympify(upper_bound)),
    )
    return CollectiveCoordinateMetric(
        profile_derivative=derivative,
        geometric_integral=geometric_integral,
        inertia=sp.simplify(coefficient * geometric_integral),
    )


def reduced_collective_lagrangian(
    inertia: Any,
    potential: Any,
    velocity: Any,
) -> sp.Expr:
    """Return the declared reduced form ``M*qdot**2/2-U``."""

    metric = _positive_quantity(inertia, "inertia")
    energy = _real_quantity(potential, "potential")
    rate = _real_quantity(velocity, "velocity")
    return sp.simplify(metric * rate**2 / 2 - energy)


def reduced_collective_euler_lagrange(
    inertia: Any,
    potential: Any,
    coordinate: sp.Symbol,
    velocity: Any,
    acceleration: Any,
) -> sp.Expr:
    r"""Return ``M*qddot + M_q*qdot**2/2 + U_q``.

    This is the Euler--Lagrange left-hand side of
    ``L=M(q)*qdot**2/2-U(q)``.  A coordinate-dependent metric therefore adds
    the displayed connection term.
    """

    metric = _positive_quantity(inertia, "inertia")
    energy = _real_quantity(potential, "potential")
    rate = _real_quantity(velocity, "velocity")
    acceleration_value = _real_quantity(acceleration, "acceleration")
    return sp.simplify(
        metric * acceleration_value
        + sp.diff(metric, coordinate) * rate**2 / 2
        + sp.diff(energy, coordinate)
    )


def stationary_collective_linearization(
    inertia_at_stationary_point: Any,
    potential_curvature: Any,
) -> CollectiveCoordinateLinearization:
    r"""Classify ``M*delta_qddot+U''*delta_q=0`` for positive ``M``."""

    metric = _positive_quantity(
        inertia_at_stationary_point,
        "inertia_at_stationary_point",
    )
    curvature = _real_quantity(potential_curvature, "potential_curvature")
    ratio = sp.simplify(curvature / metric)
    root_squared = sp.simplify(-ratio)
    root = sp.sqrt(root_squared)

    if curvature.is_positive is True:
        stability = "stable"
        stable_frequency = sp.sqrt(ratio)
        instability_rate = None
    elif curvature.is_zero is True:
        stability = "neutral"
        stable_frequency = sp.S.Zero
        instability_rate = sp.S.Zero
    elif curvature.is_negative is True:
        stability = "unstable"
        stable_frequency = None
        instability_rate = sp.sqrt(root_squared)
    else:
        stability = "undetermined"
        stable_frequency = None
        instability_rate = None

    return CollectiveCoordinateLinearization(
        inertia=metric,
        curvature=curvature,
        spectral_ratio=ratio,
        characteristic_root_squared=root_squared,
        characteristic_roots=(-root, root),
        stability=stability,
        stable_angular_frequency=stable_frequency,
        instability_rate=instability_rate,
    )


def transform_collective_coordinate(
    inertia: Any,
    potential: Any,
    coordinate: sp.Symbol,
    new_coordinate: sp.Symbol,
    inverse_map: Any,
) -> CollectiveCoordinateReparameterization:
    r"""Transform a reduced model under the local inverse ``q=g(Q)``.

    At a stationary point, the returned ``stationary_hessian`` omits the
    gradient-proportional chain term and its ratio to the transformed inertia
    equals the old ``U_qq/M`` evaluated at ``q=g(Q)``.
    """

    metric = _positive_quantity(inertia, "inertia")
    energy = _real_quantity(potential, "potential")
    mapping = _real_quantity(inverse_map, "inverse_map")
    jacobian = sp.diff(mapping, new_coordinate)
    if jacobian.is_zero is True:
        raise ValueError("inverse_map must be locally invertible")
    second_jacobian = sp.diff(mapping, new_coordinate, 2)
    transformed_metric = sp.simplify(
        metric.subs(coordinate, mapping) * jacobian**2
    )
    transformed_potential = sp.simplify(energy.subs(coordinate, mapping))
    old_gradient = sp.diff(energy, coordinate).subs(coordinate, mapping)
    old_hessian = sp.diff(energy, coordinate, 2).subs(coordinate, mapping)
    chain_term = sp.simplify(old_gradient * second_jacobian)
    stationary_hessian = sp.simplify(old_hessian * jacobian**2)
    transformed_hessian = sp.simplify(stationary_hessian + chain_term)
    return CollectiveCoordinateReparameterization(
        inverse_map=mapping,
        jacobian=jacobian,
        second_jacobian=second_jacobian,
        transformed_inertia=transformed_metric,
        transformed_potential=transformed_potential,
        transformed_hessian=transformed_hessian,
        hessian_chain_term=chain_term,
        stationary_hessian=stationary_hessian,
        stationary_spectral_ratio=sp.simplify(
            stationary_hessian / transformed_metric
        ),
    )


def capillary_barrier_top_linearization(
    pressure: Any,
    inertia_at_barrier_top: Any,
) -> CollectiveCoordinateLinearization:
    """Return the unstable capillary-top linearization for positive ``P,M``."""

    drive = _positive_quantity(pressure, "pressure")
    return stationary_collective_linearization(
        inertia_at_barrier_top,
        -2 * sp.pi * drive,
    )


def collective_coordinate_dimension_ledger() -> CollectiveCoordinateDimensionLedger:
    """Return the exact ``(E,L,T)`` dimension ledger for a radius pullback."""

    names = (
        "inertia_density",
        "profile_derivative",
        "spatial_measure",
        "geometric_integral",
        "collective_inertia",
        "potential_curvature",
        "spectral_ratio",
        "rate",
    )
    matrix = sp.Matrix(
        [
            [1, 0, 0, 0, 1, 1, 0, 0],
            [-1, -1, 1, -1, -2, -2, 0, 0],
            [2, 0, 0, 0, 2, 0, -2, -1],
        ]
    )
    return CollectiveCoordinateDimensionLedger(
        base_dimensions=("E", "L", "T"),
        quantity_names=names,
        dimension_matrix=sp.ImmutableMatrix(matrix),
    )
