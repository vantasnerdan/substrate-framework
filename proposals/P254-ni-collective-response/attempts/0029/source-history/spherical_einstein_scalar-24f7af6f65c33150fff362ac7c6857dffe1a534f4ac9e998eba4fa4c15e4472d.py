"""Exact phase-averaged spherical Einstein--sine-Gordon reduction.

The physical theory is the canonical scalar sector declared by
``einstein_scalar`` with

``phi = F*u`` and ``V(phi) = mu**2*F**2*(1-cos(u))``.

Using ``x=mu*r``, ``tau=mu*t``, ``alpha=kappa*F**2`` and the static areal
metric ``ds^2=-exp(2*Phi) d tau^2 + dx^2/f + x^2 d Omega^2`` with
``f=1-2*m/x``, this module derives the single-harmonic, phase-averaged reduced
equations for ``u=a(x)*cos(Omega*tau)``.  It deliberately also exposes the
discarded harmonics: a solution of this reduction is not thereby a pointwise
solution of the full time-dependent Einstein--scalar equations.
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


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    result = _exact_real(value, name)
    if result.is_nonnegative is not True:
        raise ValueError(f"{name} must be declared nonnegative")
    return result


@dataclass(frozen=True)
class SineGordonGravityScaling:
    """Natural-unit nondimensionalization of the canonical physical action."""

    gravitational_coupling: sp.Expr
    field_scale: sp.Expr
    mass_scale: sp.Expr
    dimensionless_field: sp.Expr
    physical_radius: sp.Expr
    physical_time: sp.Expr
    geometric_mass: sp.Expr
    physical_field: sp.Expr
    physical_potential: sp.Expr
    dimensionless_radius: sp.Expr
    dimensionless_time: sp.Expr
    dimensionless_mass: sp.Expr
    dimensionless_coupling: sp.Expr
    mass_dimensions: tuple[tuple[str, int], ...]


@dataclass(frozen=True)
class AveragedSphericalScalarStress:
    """Dimensionless orthonormal stress averaged over the scalar phase."""

    energy_density: sp.Expr
    radial_pressure: sp.Expr
    tangential_pressure: sp.Expr


@dataclass(frozen=True)
class StaticSphericalHarmonicReduction:
    """Exact equations and truncation defects for the averaged ansatz."""

    radius: sp.Symbol
    amplitude: sp.Expr
    mass: sp.Expr
    lapse_exponent: sp.Expr
    frequency: sp.Expr
    dimensionless_coupling: sp.Expr
    lapse: sp.Expr
    radial_metric_function: sp.Expr
    stress: AveragedSphericalScalarStress
    mass_constraint_residual: sp.Expr
    lapse_constraint_residual: sp.Expr
    scalar_equation_residual: sp.Expr
    conservation_residual: sp.Expr
    conservation_factor: sp.Expr
    conservation_identity_residual: sp.Expr
    discarded_scalar_third_harmonic: sp.Expr
    pointwise_energy_density_second_harmonic: sp.Expr


@dataclass(frozen=True)
class RegularOriginSphericalData:
    """Leading regular-origin coefficients for fixed central amplitude."""

    central_amplitude: sp.Expr
    central_lapse_exponent: sp.Expr
    central_energy_density: sp.Expr
    central_radial_pressure: sp.Expr
    amplitude_second_derivative: sp.Expr
    mass_cubic_coefficient: sp.Expr
    lapse_second_derivative: sp.Expr


def sine_gordon_gravity_scaling(
    gravitational_coupling: Any,
    field_scale: Any,
    mass_scale: Any,
    dimensionless_field: Any,
    physical_radius: Any,
    physical_time: Any,
    geometric_mass: Any,
) -> SineGordonGravityScaling:
    """Return the exact scale map used by the reduced equations.

    Natural mass dimensions in four spacetime dimensions are
    ``[kappa]=-2``, ``[F]=[mu]=[phi]=1``, ``[r]=[t]=[M_geo]=-1`` and
    ``[V]=4``.  Thus every returned dimensionless variable has dimension zero.
    ``geometric_mass`` is the length-valued areal mass appearing in
    ``1-2*M_geo/r``; it is not a particle mass or a fitted physical constant.
    """

    kappa = _positive_exact(gravitational_coupling, "gravitational_coupling")
    scale = _positive_exact(field_scale, "field_scale")
    mu = _positive_exact(mass_scale, "mass_scale")
    field = _exact_real(dimensionless_field, "dimensionless_field")
    radius = _exact_real(physical_radius, "physical_radius")
    time = _exact_real(physical_time, "physical_time")
    mass = _exact_real(geometric_mass, "geometric_mass")

    return SineGordonGravityScaling(
        gravitational_coupling=kappa,
        field_scale=scale,
        mass_scale=mu,
        dimensionless_field=field,
        physical_radius=radius,
        physical_time=time,
        geometric_mass=mass,
        physical_field=sp.simplify(scale * field),
        physical_potential=sp.simplify(mu**2 * scale**2 * (1 - sp.cos(field))),
        dimensionless_radius=sp.simplify(mu * radius),
        dimensionless_time=sp.simplify(mu * time),
        dimensionless_mass=sp.simplify(mu * mass),
        dimensionless_coupling=sp.simplify(kappa * scale**2),
        mass_dimensions=(
            ("gravitational_coupling", -2),
            ("field_scale", 1),
            ("mass_scale", 1),
            ("physical_field", 1),
            ("physical_potential", 4),
            ("physical_radius", -1),
            ("physical_time", -1),
            ("geometric_mass", -1),
            ("dimensionless_radius", 0),
            ("dimensionless_time", 0),
            ("dimensionless_mass", 0),
            ("dimensionless_coupling", 0),
        ),
    )


def static_spherical_sine_gordon_reduction(
    radius: sp.Symbol,
    amplitude: Any,
    mass: Any,
    lapse_exponent: Any,
    frequency: Any,
    dimensionless_coupling: Any,
) -> StaticSphericalHarmonicReduction:
    r"""Derive the exact phase-averaged single-harmonic reduced equations.

    The result is valid on a domain where ``radius>0``, ``exp(Phi)>0`` and
    ``f=1-2*m/radius>0``.  Constraint residuals vanish for a reduced solution.
    The conservation identity is derived from the same projected scalar
    equation and is returned as a consistency identity, not an independent
    oracle for the full Einstein--scalar PDE.
    """

    if not isinstance(radius, sp.Symbol) or radius.is_positive is not True:
        raise ValueError("radius must be a positive SymPy Symbol")
    x = radius
    a = _exact_real(amplitude, "amplitude")
    m = _exact_real(mass, "mass")
    phi = _exact_real(lapse_exponent, "lapse_exponent")
    omega = _positive_exact(frequency, "frequency")
    alpha = _nonnegative_exact(dimensionless_coupling, "dimensionless_coupling")

    lapse = sp.exp(phi)
    radial_metric = sp.simplify(1 - 2 * m / x)
    if radial_metric.is_nonpositive is True:
        raise ValueError("radial_metric_function must be positive on the domain")
    first = sp.diff(a, x)

    density = sp.simplify(
        omega**2 * a**2 / (4 * lapse**2)
        + radial_metric * first**2 / 4
        + 1
        - sp.besselj(0, a)
    )
    radial_pressure = sp.simplify(
        radial_metric * first**2 / 4
        + omega**2 * a**2 / (4 * lapse**2)
        - (1 - sp.besselj(0, a))
    )
    tangential_pressure = sp.simplify(
        omega**2 * a**2 / (4 * lapse**2)
        - radial_metric * first**2 / 4
        - (1 - sp.besselj(0, a))
    )
    stress = AveragedSphericalScalarStress(
        energy_density=density,
        radial_pressure=radial_pressure,
        tangential_pressure=tangential_pressure,
    )

    mass_constraint = sp.simplify(sp.diff(m, x) - alpha * x**2 * density / 2)
    lapse_constraint = sp.simplify(
        sp.diff(phi, x)
        - (m + alpha * x**3 * radial_pressure / 2) / (x * (x - 2 * m))
    )
    scalar_residual = sp.simplify(
        sp.diff(a, x, 2)
        + (
            sp.diff(phi, x)
            + sp.diff(radial_metric, x) / (2 * radial_metric)
            + 2 / x
        )
        * first
        + (omega**2 * a / lapse**2 - 2 * sp.besselj(1, a))
        / radial_metric
    )
    conservation = sp.simplify(
        sp.diff(radial_pressure, x)
        + (density + radial_pressure) * sp.diff(phi, x)
        + 2 * (radial_pressure - tangential_pressure) / x
    )
    factor = sp.simplify(radial_metric * first * scalar_residual / 2)
    identity = sp.simplify(sp.expand_func(conservation - factor))
    if identity != 0:
        raise AssertionError("averaged stress conservation does not factorize")

    return StaticSphericalHarmonicReduction(
        radius=x,
        amplitude=a,
        mass=m,
        lapse_exponent=phi,
        frequency=omega,
        dimensionless_coupling=alpha,
        lapse=lapse,
        radial_metric_function=radial_metric,
        stress=stress,
        mass_constraint_residual=mass_constraint,
        lapse_constraint_residual=lapse_constraint,
        scalar_equation_residual=scalar_residual,
        conservation_residual=conservation,
        conservation_factor=factor,
        conservation_identity_residual=identity,
        discarded_scalar_third_harmonic=sp.simplify(2 * sp.besselj(3, a)),
        pointwise_energy_density_second_harmonic=sp.simplify(
            -omega**2 * a**2 / (4 * lapse**2)
            + radial_metric * first**2 / 4
            + 2 * sp.besselj(2, a)
        ),
    )


def regular_origin_sine_gordon_data(
    central_amplitude: Any,
    central_lapse_exponent: Any,
    frequency: Any,
    dimensionless_coupling: Any,
) -> RegularOriginSphericalData:
    """Return the regular Taylor data at the areal origin.

    With the returned values, ``a=A+a_second*x**2/2+O(x**4)``,
    ``m=mass_cubic*x**3+O(x**5)`` and
    ``Phi=Phi0+lapse_second*x**2/2+O(x**4)``.
    """

    amplitude = _exact_real(central_amplitude, "central_amplitude")
    phi_zero = _exact_real(central_lapse_exponent, "central_lapse_exponent")
    omega = _positive_exact(frequency, "frequency")
    alpha = _nonnegative_exact(dimensionless_coupling, "dimensionless_coupling")
    lapse_zero = sp.exp(phi_zero)
    density = sp.simplify(
        omega**2 * amplitude**2 / (4 * lapse_zero**2)
        + 1
        - sp.besselj(0, amplitude)
    )
    radial_pressure = sp.simplify(
        omega**2 * amplitude**2 / (4 * lapse_zero**2)
        - (1 - sp.besselj(0, amplitude))
    )
    amplitude_second = sp.simplify(
        (2 * sp.besselj(1, amplitude) - omega**2 * amplitude / lapse_zero**2) / 3
    )
    mass_cubic = sp.simplify(alpha * density / 6)
    lapse_second = sp.simplify(alpha * (density / 6 + radial_pressure / 2))
    return RegularOriginSphericalData(
        central_amplitude=amplitude,
        central_lapse_exponent=phi_zero,
        central_energy_density=density,
        central_radial_pressure=radial_pressure,
        amplitude_second_derivative=amplitude_second,
        mass_cubic_coefficient=mass_cubic,
        lapse_second_derivative=lapse_second,
    )
