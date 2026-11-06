"""Exact canonical Einstein-scalar ledgers and a coupled FLRW solution.

The declared mostly-plus action is

``S = integral sqrt(-g) * (R/(2*kappa) - (d phi)^2/2 - V(phi))``.

The module exposes the scalar stress fixed by that action and a complete exact
nonvacuum solution for ``V=0``.  It does not identify the scalar with a
refractive index or breather, select a physical gravitational coupling, or
claim that nonzero stress alone determines a geometry.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _exact_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if result.is_real is not True:
        raise ValueError(f"{name} must be declared real")
    return sp.simplify(result)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    result = _exact_real(value, name)
    if result.is_positive is not True:
        raise ValueError(f"{name} must be declared positive")
    return result


def _exact_matrix(value: Any, rows: int, columns: int, name: str) -> sp.Matrix:
    result = sp.Matrix(value)
    if result.shape != (rows, columns):
        raise ValueError(f"{name} must have shape {rows} by {columns}")
    for entry in result:
        _exact_real(entry, name)
    return result.applyfunc(sp.simplify)


@dataclass(frozen=True)
class ScalarStressLedger:
    """Covariant stress data for ``L=-epsilon*(d phi)^2/2-V``."""

    metric_covariant: sp.Matrix
    metric_contravariant: sp.Matrix
    gradient_covariant: sp.Matrix
    kinetic_coefficient: sp.Expr
    potential: sp.Expr
    gradient_norm_squared: sp.Expr
    covariant: sp.Matrix
    mixed: sp.Matrix
    trace: sp.Expr


@dataclass(frozen=True)
class MasslessScalarFlatFLRW:
    """Expanding or contracting scalar branch on the domain ``time>0``."""

    gravitational_coupling: sp.Expr
    time: sp.Expr
    reference_time: sp.Expr
    reference_scale_factor: sp.Expr
    reference_scalar: sp.Expr
    branch: sp.Expr
    scale_factor: sp.Expr
    scalar_field: sp.Expr
    hubble_rate: sp.Expr
    scalar_velocity: sp.Expr
    energy_density: sp.Expr
    pressure: sp.Expr
    friedmann_residual: sp.Expr
    spatial_einstein_residual: sp.Expr
    scalar_equation_residual: sp.Expr
    continuity_residual: sp.Expr
    ricci_scalar: sp.Expr
    kretschmann_scalar: sp.Expr


def minimally_coupled_scalar_stress(
    metric_covariant: Any,
    gradient_covariant: Any,
    potential: Any,
    *,
    kinetic_coefficient: Any = 1,
) -> ScalarStressLedger:
    """Return the exact scalar stress for a declared four-dimensional metric.

    The action density is
    ``-epsilon*g^ab*phi_a*phi_b/2 - V``.  The canonical healthy convention has
    ``epsilon=1``; a negative value is retained only for explicit ghost
    countermodels.
    """

    metric = _exact_matrix(metric_covariant, 4, 4, "metric_covariant")
    if metric != metric.T:
        raise ValueError("metric_covariant must be symmetric")
    determinant = sp.simplify(metric.det())
    if determinant == 0:
        raise ValueError("metric_covariant must be invertible")
    inverse = metric.inv().applyfunc(sp.simplify)
    gradient = _exact_matrix(gradient_covariant, 4, 1, "gradient_covariant")
    scalar_potential = _exact_real(potential, "potential")
    epsilon = _exact_real(kinetic_coefficient, "kinetic_coefficient")
    if epsilon == 0:
        raise ValueError("kinetic_coefficient must be nonzero")

    norm = sp.simplify((gradient.T * inverse * gradient)[0])
    stress = (
        epsilon * gradient * gradient.T
        - metric * (epsilon * norm / 2 + scalar_potential)
    ).applyfunc(sp.simplify)
    mixed = (inverse * stress).applyfunc(sp.simplify)
    trace = sp.simplify(sp.trace(mixed))
    expected_trace = sp.simplify(-epsilon * norm - 4 * scalar_potential)
    if sp.simplify(trace - expected_trace) != 0:
        raise AssertionError("scalar stress trace does not close")

    return ScalarStressLedger(
        metric_covariant=metric,
        metric_contravariant=inverse,
        gradient_covariant=gradient,
        kinetic_coefficient=epsilon,
        potential=scalar_potential,
        gradient_norm_squared=norm,
        covariant=stress,
        mixed=mixed,
        trace=trace,
    )


def massless_scalar_flat_flrw(
    gravitational_coupling: Any,
    time: Any,
    reference_time: Any,
    reference_scale_factor: Any,
    reference_scalar: Any,
    *,
    branch: Any = 1,
) -> MasslessScalarFlatFLRW:
    r"""Return an exact nonvacuum solution of the declared ``V=0`` action.

    On ``t>0``, the spatially flat metric has
    ``a=a0*(t/t0)^(1/3)`` and the homogeneous scalar is
    ``phi=phi0+s*sqrt(2/(3*kappa))*log(t/t0)`` for ``s=+1`` or ``-1``.
    The function checks the two independent Einstein equations, scalar equation,
    and continuity equation before returning.  The solution is singular at
    ``t=0`` and makes no breather or refractive-index identification.
    """

    kappa = _positive_exact(gravitational_coupling, "gravitational_coupling")
    if not isinstance(time, sp.Symbol):
        raise ValueError("time must be a positive SymPy Symbol")
    current_time = _positive_exact(time, "time")
    time_zero = _positive_exact(reference_time, "reference_time")
    scale_zero = _positive_exact(reference_scale_factor, "reference_scale_factor")
    scalar_zero = _exact_real(reference_scalar, "reference_scalar")
    sign = _exact_real(branch, "branch")
    if sign not in {-1, 1}:
        raise ValueError("branch must be exactly -1 or 1")

    scale_factor = sp.simplify(scale_zero * (current_time / time_zero) ** sp.Rational(1, 3))
    scalar_field = sp.simplify(
        scalar_zero
        + sign * sp.sqrt(sp.Rational(2, 3) / kappa) * sp.log(current_time / time_zero)
    )
    hubble = sp.simplify(sp.diff(scale_factor, current_time) / scale_factor)
    scalar_velocity = sp.simplify(sp.diff(scalar_field, current_time))
    energy_density = sp.simplify(scalar_velocity**2 / 2)
    pressure = energy_density
    friedmann = sp.simplify(3 * hubble**2 - kappa * energy_density)
    spatial = sp.simplify(
        -(2 * sp.diff(hubble, current_time) + 3 * hubble**2) - kappa * pressure
    )
    scalar_residual = sp.simplify(
        sp.diff(scalar_field, current_time, 2) + 3 * hubble * scalar_velocity
    )
    continuity = sp.simplify(
        sp.diff(energy_density, current_time)
        + 3 * hubble * (energy_density + pressure)
    )
    if any(value != 0 for value in (friedmann, spatial, scalar_residual, continuity)):
        raise AssertionError("massless Einstein-scalar FLRW equations do not close")

    ricci = sp.simplify(6 * (sp.diff(hubble, current_time) + 2 * hubble**2))
    acceleration_ratio = sp.simplify(sp.diff(scale_factor, current_time, 2) / scale_factor)
    kretschmann = sp.simplify(12 * (acceleration_ratio**2 + hubble**4))
    return MasslessScalarFlatFLRW(
        gravitational_coupling=kappa,
        time=current_time,
        reference_time=time_zero,
        reference_scale_factor=scale_zero,
        reference_scalar=scalar_zero,
        branch=sign,
        scale_factor=scale_factor,
        scalar_field=scalar_field,
        hubble_rate=hubble,
        scalar_velocity=scalar_velocity,
        energy_density=energy_density,
        pressure=pressure,
        friedmann_residual=friedmann,
        spatial_einstein_residual=spatial,
        scalar_equation_residual=scalar_residual,
        continuity_residual=continuity,
        ricci_scalar=ricci,
        kretschmann_scalar=kretschmann,
    )
