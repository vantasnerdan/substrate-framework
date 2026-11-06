"""Integrated moment identities for localized conserved stress tensors.

The APIs encode kinematics and conservation in flat three-dimensional space.
They do not decide whether a moment sources or radiates any gravitational
field.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _column_three(vector: Any, name: str) -> sp.Matrix:
    value = sp.Matrix(vector)
    if value.shape not in ((3, 1), (1, 3)):
        raise ValueError(f"{name} must have three components")
    return value if value.shape == (3, 1) else value.T


def symmetric_trace_free(matrix: Any) -> sp.Matrix:
    """Return ``(A+A.T)/2 - I*trace(A)/3`` exactly."""

    value = sp.Matrix(matrix)
    if value.shape != (3, 3):
        raise ValueError("matrix must be 3 by 3")
    symmetric = sp.simplify((value + value.T) / 2)
    return sp.simplify(symmetric - sp.eye(3) * sp.trace(symmetric) / 3)


@dataclass(frozen=True)
class DiscreteMassMoments:
    """Cartesian moments of finitely many declared point masses."""

    monopole: sp.Expr
    dipole: sp.Matrix
    second_moment: sp.Matrix
    trace_free_second_moment: sp.Matrix
    triple_normalized_quadrupole: sp.Matrix


@dataclass(frozen=True)
class RadialDensitySecondMoments:
    """Second moments of a radial density with an optional axisymmetric P2 factor."""

    scalar_radial_moment: sp.Expr
    second_moment: sp.Matrix
    trace_free_second_moment: sp.Matrix
    triple_normalized_quadrupole: sp.Matrix


def axisymmetric_p2_density_second_moments(
    scalar_radial_moment: Any,
    deformation_amplitude: Any = 0,
) -> RadialDensitySecondMoments:
    """Return moments for ``rho=f(r)*(1+a*P2(cos(theta)))`` exactly.

    ``scalar_radial_moment`` means ``integral rho*r^2*d^3x``.  Because the
    angular mean of ``P2`` vanishes, it is also ``4*pi*integral f(r)*r^4 dr``
    for every ``a``.  The symmetry axis is z.  Setting ``a=0`` gives an
    arbitrary spherical density; no field equation or positivity assumption
    is needed for this angular theorem.
    """

    scalar = sp.sympify(scalar_radial_moment)
    amplitude = sp.sympify(deformation_amplitude)
    second = sp.diag(
        scalar * (sp.Rational(1, 3) - amplitude / 15),
        scalar * (sp.Rational(1, 3) - amplitude / 15),
        scalar * (sp.Rational(1, 3) + 2 * amplitude / 15),
    )
    trace_free = symmetric_trace_free(second)
    return RadialDensitySecondMoments(
        scalar_radial_moment=scalar,
        second_moment=second,
        trace_free_second_moment=trace_free,
        triple_normalized_quadrupole=sp.simplify(3 * trace_free),
    )


def spherical_density_second_moments(
    scalar_radial_moment: Any,
) -> RadialDensitySecondMoments:
    """Return the isotropic second moment and exact STF null of a radial density."""

    return axisymmetric_p2_density_second_moments(scalar_radial_moment, 0)


def discrete_mass_moments(
    masses: Sequence[Any],
    positions: Sequence[Sequence[Any]],
) -> DiscreteMassMoments:
    """Return ``M``, ``D_i``, ``I_ij``, and two STF conventions."""

    if len(masses) == 0 or len(masses) != len(positions):
        raise ValueError("masses and positions must have equal nonzero length")
    mass_values = tuple(sp.sympify(mass) for mass in masses)
    position_values = tuple(
        _column_three(position, "each position") for position in positions
    )
    monopole = sp.simplify(sum(mass_values, sp.Integer(0)))
    dipole = sp.simplify(
        sum(
            (mass * position for mass, position in zip(mass_values, position_values)),
            sp.zeros(3, 1),
        )
    )
    second = sp.simplify(
        sum(
            (
                mass * position * position.T
                for mass, position in zip(mass_values, position_values)
            ),
            sp.zeros(3),
        )
    )
    trace_free = symmetric_trace_free(second)
    return DiscreteMassMoments(
        monopole=monopole,
        dipole=dipole,
        second_moment=second,
        trace_free_second_moment=trace_free,
        triple_normalized_quadrupole=sp.simplify(3 * trace_free),
    )


@dataclass(frozen=True)
class IsolatedMomentRates:
    """Moment derivatives implied by localized symmetric stress conservation."""

    monopole_rate: sp.Expr
    momentum_rate: sp.Matrix
    dipole_rate: sp.Matrix
    dipole_acceleration: sp.Matrix
    second_moment_acceleration: sp.Matrix
    trace_free_second_moment_acceleration: sp.Matrix
    triple_normalized_quadrupole_acceleration: sp.Matrix


def isolated_conserved_stress_moment_rates(
    total_momentum: Any,
    integrated_spatial_stress: Any,
) -> IsolatedMomentRates:
    """Return exact isolated-source rates through the second mass moment.

    Assumptions are ``partial_mu T^{mu nu}=0``, ``T^{mu nu}=T^{nu mu}``, and
    vanishing surface terms, including the coordinate-weighted terms needed
    for second moments.  The convention is
    ``I_ij=integral T00*x_i*x_j`` and therefore
    ``ddot(I)_ij=2*integral T_ij``.
    """

    momentum = _column_three(total_momentum, "total_momentum")
    stress = sp.Matrix(integrated_spatial_stress)
    if stress.shape != (3, 3):
        raise ValueError("integrated_spatial_stress must be 3 by 3")
    if sp.simplify(stress - stress.T) != sp.zeros(3):
        raise ValueError("integrated_spatial_stress must be symmetric")
    second_acceleration = sp.simplify(2 * stress)
    trace_free_acceleration = symmetric_trace_free(second_acceleration)
    return IsolatedMomentRates(
        monopole_rate=sp.Integer(0),
        momentum_rate=sp.zeros(3, 1),
        dipole_rate=momentum,
        dipole_acceleration=sp.zeros(3, 1),
        second_moment_acceleration=second_acceleration,
        trace_free_second_moment_acceleration=trace_free_acceleration,
        triple_normalized_quadrupole_acceleration=sp.simplify(
            3 * trace_free_acceleration
        ),
    )
