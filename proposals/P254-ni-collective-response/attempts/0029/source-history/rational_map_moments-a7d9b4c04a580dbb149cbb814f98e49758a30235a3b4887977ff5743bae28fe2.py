"""Conditional intrinsic moments of declared rational-map energy densities.

The APIs factor a separately declared rational-map local density into radial
and angular moments.  They do not turn the rational-map approximation into a
full three-dimensional field solution, a conserved physical stress tensor, a
particle or nucleus, or a gravitational radiation source.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
import sympy as sp

from .numerics import trapezoid_integral
from .rational_map_radial import RationalMapRadialProfileEvidence

FloatMatrix = NDArray[np.float64]


@dataclass(frozen=True)
class RationalMapAngularSTFMoments:
    """Exact sphere averages and unnormalized angular STF tensors."""

    degree: int
    conformal_jacobian: sp.Expr
    normalized_mean: sp.Expr
    normalized_squared_mean: sp.Expr
    linear_stf: sp.Matrix
    quadratic_stf: sp.Matrix


@dataclass(frozen=True)
class FactorizedRationalMapEnergyMoments:
    """Exact monopole and two accepted STF conventions after factorization."""

    monopole: sp.Expr
    normalized_stf: sp.Matrix
    triple_normalized_quadrupole: sp.Matrix


@dataclass(frozen=True)
class RationalMapProfileMomentEvidence:
    """Resolution-bounded intrinsic moments of one declared radial branch."""

    degree: int
    angular_integral: float
    monopole: float
    normalized_stf: FloatMatrix
    triple_normalized_quadrupole: FloatMatrix
    linear_second_radial_integral: float
    quadratic_second_radial_integral: float
    origin_monopole_estimate: float
    tail_monopole_estimate: float
    origin_second_moment_estimate: float
    tail_second_moment_estimate: float
    profile_energy_closure_relative_error: float
    coordinate_units: str = "declared dimensionless radial coordinate squared"

    @property
    def normalized_axial_ratio(self) -> float:
        """Return ``I_STF_zz/M`` in declared coordinate-squared units."""

        return float(self.normalized_stf[2, 2] / self.monopole)

    @property
    def triple_axial_ratio(self) -> float:
        """Return ``Q_zz/M=3*I_STF_zz/M`` in the accepted Q convention."""

        return float(self.triple_normalized_quadrupole[2, 2] / self.monopole)


def rational_map_local_energy_density(
    field: Any,
    radial_derivative: Any,
    radius: Any,
    conformal_jacobian: Any,
) -> sp.Expr:
    r"""Return the declared angular-resolved rational-map density.

    With ``J`` denoting the conformal Jacobian, the declaration is

    ``f'^2 + 2*sin(f)^2*(1+f'^2)*J/r^2 + sin(f)^4*J^2/r^4``.

    Its sphere average reproduces the accepted reduced radial density only
    when the separately established inputs ``<J>=B`` and ``<J^2>=I`` hold.
    """

    f = sp.sympify(field)
    fp = sp.sympify(radial_derivative)
    r = sp.sympify(radius)
    jacobian = sp.sympify(conformal_jacobian)
    if r.is_number and r.is_positive is not True:
        raise ValueError("radius must be positive")
    return sp.simplify(
        fp**2
        + 2 * sp.sin(f) ** 2 * (1 + fp**2) * jacobian / r**2
        + sp.sin(f) ** 4 * jacobian**2 / r**4
    )


def degree_one_rational_map_angular_stf_moments() -> RationalMapAngularSTFMoments:
    """Return the exact isotropic angular data for ``R(z)=z``."""

    zero = sp.zeros(3)
    return RationalMapAngularSTFMoments(
        degree=1,
        conformal_jacobian=sp.Integer(1),
        normalized_mean=sp.Integer(1),
        normalized_squared_mean=sp.Integer(1),
        linear_stf=zero,
        quadratic_stf=zero,
    )


def degree_two_axial_rational_map_angular_stf_moments() -> RationalMapAngularSTFMoments:
    r"""Return exact angular STF tensors for the axial map ``R(z)=z^2``.

    In ``u=cos(theta)``, its conformal Jacobian is
    ``J=4*(1-u^2)/(1+u^2)^2``.  The returned tensors are the full-sphere
    integrals of ``J^k*(n_i*n_j-delta_ij/3)`` for ``k=1,2``; they are not
    normalized sphere averages.
    """

    u = sp.Symbol("u", real=True)
    jacobian = 4 * (1 - u**2) / (1 + u**2) ** 2
    linear_zz = sp.Rational(8, 3) * sp.pi * (3 * sp.pi - 10)
    quadratic_zz = sp.Rational(8, 9) * sp.pi * (3 * sp.pi - 16)
    linear = sp.diag(-linear_zz / 2, -linear_zz / 2, linear_zz)
    quadratic = sp.diag(
        -quadratic_zz / 2,
        -quadratic_zz / 2,
        quadratic_zz,
    )
    return RationalMapAngularSTFMoments(
        degree=2,
        conformal_jacobian=jacobian,
        normalized_mean=sp.Integer(2),
        normalized_squared_mean=sp.pi + sp.Rational(8, 3),
        linear_stf=sp.simplify(linear),
        quadratic_stf=sp.simplify(quadratic),
    )


def factorized_rational_map_energy_moments(
    angular: RationalMapAngularSTFMoments,
    *,
    isotropic_monopole_radial: Any,
    linear_monopole_radial: Any,
    quadratic_monopole_radial: Any,
    linear_second_radial: Any,
    quadratic_second_radial: Any,
) -> FactorizedRationalMapEnergyMoments:
    r"""Combine radial integrals with exact angular data.

    The radial monopole inputs multiply ``1``, ``J``, and ``J^2`` after the
    volume factor is included.  The two second-moment inputs multiply the
    angular STF tensors of ``J`` and ``J^2``; the isotropic term vanishes.
    """

    m0 = sp.sympify(isotropic_monopole_radial)
    m1 = sp.sympify(linear_monopole_radial)
    m2 = sp.sympify(quadratic_monopole_radial)
    h1 = sp.sympify(linear_second_radial)
    h2 = sp.sympify(quadratic_second_radial)
    monopole = sp.simplify(
        4
        * sp.pi
        * (
            m0
            + angular.normalized_mean * m1
            + angular.normalized_squared_mean * m2
        )
    )
    normalized = sp.simplify(h1 * angular.linear_stf + h2 * angular.quadratic_stf)
    return FactorizedRationalMapEnergyMoments(
        monopole=monopole,
        normalized_stf=normalized,
        triple_normalized_quadrupole=sp.simplify(3 * normalized),
    )


def _profile_arrays(
    profile: RationalMapRadialProfileEvidence,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    if profile.degree != 2:
        raise ValueError("profile must have declared degree two")
    expected_angular = float(sp.pi + sp.Rational(8, 3))
    if not np.isclose(profile.angular_integral, expected_angular, rtol=0.0, atol=2e-13):
        raise ValueError("profile must use the exact axial degree-two angular integral")
    radius = np.asarray(profile.radius, dtype=np.float64)
    field = np.asarray(profile.field, dtype=np.float64)
    derivative = np.asarray(profile.radial_derivative, dtype=np.float64)
    if radius.ndim != 1 or radius.size < 101 or not (
        radius.shape == field.shape == derivative.shape
    ):
        raise ValueError("profile arrays must be equal one-dimensional sampled vectors")
    if (
        not np.all(np.isfinite(radius))
        or not np.all(np.isfinite(field))
        or not np.all(np.isfinite(derivative))
    ):
        raise ValueError("profile arrays must be finite")
    if np.any(radius <= 0.0) or np.any(np.diff(radius) <= 0.0):
        raise ValueError("profile radius must be positive and strictly increasing")
    return radius, field, derivative


def degree_two_profile_intrinsic_moments(
    profile: RationalMapRadialProfileEvidence,
) -> RationalMapProfileMomentEvidence:
    r"""Evaluate the factorized moment of one corrected degree-two branch.

    Sampled domain terms use the canonical trapezoidal helper.  Omitted origin
    and tail intervals use the leading powers already recorded by the profile;
    they are estimates, not exact nonlinear endpoint integrals.  The returned
    ratio has units of the declared dimensionless radial coordinate squared.
    """

    radius, field, derivative = _profile_arrays(profile)
    sine_squared = np.sin(field) ** 2
    m0_domain = trapezoid_integral(radius**2 * derivative**2, radius)
    m1_domain = trapezoid_integral(
        2.0 * sine_squared * (1.0 + derivative**2),
        radius,
    )
    m2_domain = trapezoid_integral(sine_squared**2 / radius**2, radius)
    h1_domain = trapezoid_integral(
        2.0 * radius**2 * sine_squared * (1.0 + derivative**2),
        radius,
    )
    h2_domain = trapezoid_integral(sine_squared**2, radius)

    inner = float(radius[0])
    outer = float(radius[-1])
    amplitude = float(profile.origin_amplitude)
    sigma = float(profile.origin_power)
    tail_power = float(profile.tail_power)
    tail_amplitude = float(field[-1] * outer**tail_power)

    m0_origin = amplitude**2 * sigma**2 * inner ** (2 * sigma + 1) / (2 * sigma + 1)
    m1_origin = (
        2 * amplitude**2 * inner ** (2 * sigma + 1) / (2 * sigma + 1)
        + 2
        * amplitude**4
        * sigma**2
        * inner ** (4 * sigma - 1)
        / (4 * sigma - 1)
    )
    m2_origin = amplitude**4 * inner ** (4 * sigma - 1) / (4 * sigma - 1)
    h1_origin = (
        2 * amplitude**2 * inner ** (2 * sigma + 3) / (2 * sigma + 3)
        + 2
        * amplitude**4
        * sigma**2
        * inner ** (4 * sigma + 1)
        / (4 * sigma + 1)
    )
    h2_origin = amplitude**4 * inner ** (4 * sigma + 1) / (4 * sigma + 1)

    m0_tail = (
        tail_power**2
        * tail_amplitude**2
        * outer ** (1 - 2 * tail_power)
        / (2 * tail_power - 1)
    )
    m1_tail = (
        2
        * tail_amplitude**2
        * outer ** (1 - 2 * tail_power)
        / (2 * tail_power - 1)
        + 2
        * tail_power**2
        * tail_amplitude**4
        * outer ** (-4 * tail_power - 1)
        / (4 * tail_power + 1)
    )
    m2_tail = (
        tail_amplitude**4
        * outer ** (-4 * tail_power - 1)
        / (4 * tail_power + 1)
    )
    if tail_power <= 1.5:
        raise ValueError("degree-two second moment requires tail power above three halves")
    h1_tail = (
        2
        * tail_amplitude**2
        * outer ** (3 - 2 * tail_power)
        / (2 * tail_power - 3)
        + 2
        * tail_power**2
        * tail_amplitude**4
        * outer ** (1 - 4 * tail_power)
        / (4 * tail_power - 1)
    )
    h2_tail = (
        tail_amplitude**4
        * outer ** (1 - 4 * tail_power)
        / (4 * tail_power - 1)
    )

    angular = degree_two_axial_rational_map_angular_stf_moments()
    factored = factorized_rational_map_energy_moments(
        angular,
        isotropic_monopole_radial=m0_domain + m0_origin + m0_tail,
        linear_monopole_radial=m1_domain + m1_origin + m1_tail,
        quadratic_monopole_radial=m2_domain + m2_origin + m2_tail,
        linear_second_radial=h1_domain + h1_origin + h1_tail,
        quadratic_second_radial=h2_domain + h2_origin + h2_tail,
    )
    monopole = float(sp.N(factored.monopole, 17))
    normalized = np.asarray(factored.normalized_stf, dtype=np.float64)
    triple = np.asarray(factored.triple_normalized_quadrupole, dtype=np.float64)
    profile_energy = float(profile.two_derivative_energy + profile.four_derivative_energy)
    closure = abs(monopole - profile_energy) / profile_energy
    origin_monopole = float(
        4
        * np.pi
        * (
            m0_origin
            + 2 * m1_origin
            + float(angular.normalized_squared_mean) * m2_origin
        )
    )
    tail_monopole = float(
        4
        * np.pi
        * (
            m0_tail
            + 2 * m1_tail
            + float(angular.normalized_squared_mean) * m2_tail
        )
    )
    linear_zz = float(angular.linear_stf[2, 2])
    quadratic_zz = float(angular.quadratic_stf[2, 2])
    return RationalMapProfileMomentEvidence(
        degree=2,
        angular_integral=float(angular.normalized_squared_mean),
        monopole=monopole,
        normalized_stf=normalized,
        triple_normalized_quadrupole=triple,
        linear_second_radial_integral=float(h1_domain + h1_origin + h1_tail),
        quadratic_second_radial_integral=float(h2_domain + h2_origin + h2_tail),
        origin_monopole_estimate=origin_monopole,
        tail_monopole_estimate=tail_monopole,
        origin_second_moment_estimate=float(
            linear_zz * h1_origin + quadratic_zz * h2_origin
        ),
        tail_second_moment_estimate=float(
            linear_zz * h1_tail + quadratic_zz * h2_tail
        ),
        profile_energy_closure_relative_error=float(closure),
    )
