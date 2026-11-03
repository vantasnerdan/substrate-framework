"""Conditional Bogomolny bounds for a declared sextic-plus-potential energy.

This module owns exact algebra and convention conversions only.  It does not
establish that the declared energy is a physical action, that a saturation
equation has a solution in any degree sector, or that map degree labels a
physical state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive_symbolic(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _nonnegative_symbolic(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return expression


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be a positive integer")
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _nonzero_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be a nonzero integer")
    if value == 0:
        raise ValueError(f"{name} must be a nonzero integer")
    return value


def _orientation(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value not in (-1, 1):
        raise ValueError("orientation must be +1 or -1")
    return value


@dataclass(frozen=True)
class BogomolnyDensityDecomposition:
    """Exact pointwise square decomposition for one orientation branch."""

    energy_density: sp.Expr
    square_density: sp.Expr
    cross_density: sp.Expr
    saturation_residual: sp.Expr
    orientation: int

    @property
    def identity_residual(self) -> sp.Expr:
        """Return the exact residual of energy minus square and cross terms."""

        return sp.simplify(
            self.energy_density - self.square_density - self.cross_density
        )


@dataclass(frozen=True)
class NearBpsDifference:
    """Ledger for a signed difference of two controlled mass expansions."""

    degree_balance: int
    bps_term: sp.Expr
    linear_coefficient: sp.Expr
    remainder: sp.Expr
    expression: sp.Expr


def target_three_sphere_volume() -> sp.Expr:
    """Return the unit round three-sphere volume ``2*pi**2``."""

    return 2 * sp.pi**2


def normalized_sqrt_potential_average(sqrt_potential_integral: Any) -> sp.Expr:
    """Normalize ``integral_{S^3} sqrt(V) Omega`` by the target volume."""

    integral = _nonnegative_symbolic(
        sqrt_potential_integral,
        "sqrt_potential_integral",
    )
    return sp.simplify(integral / target_three_sphere_volume())


def degree_weighted_target_pairing(degree: int, target_average: Any) -> sp.Expr:
    """Return ``B*<sqrt(V)>`` from the declared oriented degree theorem.

    The caller supplies a target average and is responsible for the regularity
    and compactification hypotheses of the pullback integration theorem.
    """

    signed_degree = _nonzero_integer(degree, "degree")
    average = _nonnegative_symbolic(target_average, "target_average")
    return sp.expand(signed_degree * average)


def bogomolny_density_decomposition(
    baryon_density: Any,
    potential_value: Any,
    coupling_lambda: Any,
    coupling_mu: Any,
    *,
    orientation: int,
) -> BogomolnyDensityDecomposition:
    """Return the exact square completion for the declared local energy.

    ``baryon_density`` is normalized so that its spatial integral is the
    signed integer degree.  The branch orientation is normally ``sign(B)``.
    """

    density = sp.sympify(baryon_density)
    potential = _nonnegative_symbolic(potential_value, "potential_value")
    lam = _positive_symbolic(coupling_lambda, "coupling_lambda")
    mu = _positive_symbolic(coupling_mu, "coupling_mu")
    sign = _orientation(orientation)
    weighted_density = lam * sp.pi**2 * density
    weighted_potential = mu * sp.sqrt(potential)
    residual = weighted_density - sign * weighted_potential
    energy = weighted_density**2 + weighted_potential**2
    square = residual**2
    cross = 2 * sign * weighted_density * weighted_potential
    return BogomolnyDensityDecomposition(
        energy_density=energy,
        square_density=square,
        cross_density=cross,
        saturation_residual=residual,
        orientation=sign,
    )


def bps_bound_per_absolute_degree(
    coupling_lambda: Any,
    coupling_mu: Any,
    target_average: Any,
) -> sp.Expr:
    """Return the conditional bound coefficient ``2*lambda*mu*pi**2*W``."""

    lam = _positive_symbolic(coupling_lambda, "coupling_lambda")
    mu = _positive_symbolic(coupling_mu, "coupling_mu")
    average = _nonnegative_symbolic(target_average, "target_average")
    return sp.expand(2 * lam * mu * sp.pi**2 * average)


def bps_topological_lower_bound(
    degree: int,
    coupling_lambda: Any,
    coupling_mu: Any,
    target_average: Any,
) -> sp.Expr:
    """Return the degree-sector lower bound, without assuming attainment."""

    signed_degree = _nonzero_integer(degree, "degree")
    coefficient = bps_bound_per_absolute_degree(
        coupling_lambda,
        coupling_mu,
        target_average,
    )
    return sp.expand(abs(signed_degree) * coefficient)


def conditional_attained_bps_sector_energy(
    degree: int,
    coupling_lambda: Any,
    coupling_mu: Any,
    target_average: Any,
) -> sp.Expr:
    """Return the sector energy conditional on actual bound attainment.

    This deliberately named API does not infer existence from the lower bound.
    """

    return bps_topological_lower_bound(
        degree,
        coupling_lambda,
        coupling_mu,
        target_average,
    )


def near_bps_mass_difference(
    base_degree: int,
    composite_degree: int,
    *,
    multiplicity: int,
    bps_energy_per_degree: Any,
    epsilon: Any,
    base_correction: Any,
    composite_correction: Any,
    base_remainder: Any = 0,
    composite_remainder: Any = 0,
) -> NearBpsDifference:
    """Expand ``n*M(A)-M(C)`` while keeping corrections and remainders visible.

    The declared positive-degree masses are
    ``M(D)=K*D+epsilon*Delta_D+r_D``.  The BPS term cancels only when
    ``composite_degree == multiplicity*base_degree``.  An ``O(epsilon)``
    conclusion additionally requires controlled remainder inputs.
    """

    base = _positive_integer(base_degree, "base_degree")
    composite = _positive_integer(composite_degree, "composite_degree")
    count = _positive_integer(multiplicity, "multiplicity")
    coefficient = _nonnegative_symbolic(
        bps_energy_per_degree,
        "bps_energy_per_degree",
    )
    small_parameter = sp.sympify(epsilon)
    correction_base = sp.sympify(base_correction)
    correction_composite = sp.sympify(composite_correction)
    remainder_base = sp.sympify(base_remainder)
    remainder_composite = sp.sympify(composite_remainder)
    degree_balance = count * base - composite
    bps_term = sp.expand(coefficient * degree_balance)
    linear_coefficient = sp.expand(
        count * correction_base - correction_composite
    )
    remainder = sp.expand(count * remainder_base - remainder_composite)
    expression = sp.expand(
        bps_term + small_parameter * linear_coefficient + remainder
    )
    return NearBpsDifference(
        degree_balance=degree_balance,
        bps_term=bps_term,
        linear_coefficient=linear_coefficient,
        remainder=remainder,
        expression=expression,
    )
