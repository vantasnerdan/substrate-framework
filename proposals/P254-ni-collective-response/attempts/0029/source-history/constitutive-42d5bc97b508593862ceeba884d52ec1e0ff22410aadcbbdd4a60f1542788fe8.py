"""Exact utilities for a declared co-scaled electromagnetic response ansatz."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


@dataclass(frozen=True)
class SIConstitutiveDimensionLedger:
    """Exact SI dimension columns in ``(M,L,T,I)`` order."""

    base_dimensions: tuple[str, str, str, str]
    permittivity: sp.ImmutableMatrix
    permeability: sp.ImmutableMatrix
    inverse_permeability: sp.ImmutableMatrix
    mass_density: sp.ImmutableMatrix
    stiffness: sp.ImmutableMatrix
    energy_density: sp.ImmutableMatrix
    speed: sp.ImmutableMatrix
    newton_constant: sp.ImmutableMatrix
    newton_over_speed_squared: sp.ImmutableMatrix
    mechanical_conversion: sp.ImmutableMatrix


@dataclass(frozen=True)
class MechanicalMediumConversion:
    """A declared SI electromagnetic-to-mechanical conversion ledger."""

    permittivity: sp.Expr
    inverse_permeability: sp.Expr
    inertia_conversion: sp.Expr
    stiffness_conversion: sp.Expr
    strain_amplitude: sp.Expr
    mass_density: sp.Expr
    stiffness: sp.Expr
    electromagnetic_speed_squared: sp.Expr
    mechanical_speed_squared: sp.Expr
    speed_squared_ratio: sp.Expr
    strain_energy_density: sp.Expr
    mass_equivalent_density: sp.Expr


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise TypeError(f"{name} must be exact")
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise TypeError(f"{name} must be exact")
    if expression.is_number and expression.is_real is not True:
        raise ValueError(f"{name} must be real")
    return expression


def si_constitutive_dimension_ledger() -> SIConstitutiveDimensionLedger:
    """Return SI dimensions needed for an electromagnetic medium dictionary.

    Columns use mass, length, time, and electric-current order.  In particular,
    SI permittivity is not a mass density and inverse permeability is not a
    stiffness or energy density.  Both mechanical identifications require the
    same dimensioned conversion column ``M**2*T**-4*I**-2`` if their ratio is
    to retain the electromagnetic wave speed.
    """

    permittivity = sp.ImmutableMatrix([-1, -3, 4, 2])
    permeability = sp.ImmutableMatrix([1, 1, -2, -2])
    inverse_permeability = -permeability
    mass_density = sp.ImmutableMatrix([1, -3, 0, 0])
    stiffness = sp.ImmutableMatrix([1, -1, -2, 0])
    speed = sp.ImmutableMatrix([0, 1, -1, 0])
    newton = sp.ImmutableMatrix([-1, 3, -2, 0])
    conversion = sp.ImmutableMatrix(mass_density - permittivity)
    return SIConstitutiveDimensionLedger(
        base_dimensions=("M", "L", "T", "I"),
        permittivity=permittivity,
        permeability=permeability,
        inverse_permeability=sp.ImmutableMatrix(inverse_permeability),
        mass_density=mass_density,
        stiffness=stiffness,
        energy_density=stiffness,
        speed=speed,
        newton_constant=newton,
        newton_over_speed_squared=sp.ImmutableMatrix(newton - 2 * speed),
        mechanical_conversion=conversion,
    )


def mechanical_medium_conversion(
    permittivity: Any,
    inverse_permeability: Any,
    inertia_conversion: Any,
    *,
    stiffness_conversion: Any | None = None,
    strain_amplitude: Any = 1,
) -> MechanicalMediumConversion:
    """Compose a declared SI electromagnetic-to-mechanical dictionary.

    The inputs declare ``rho=a*epsilon`` and ``K=b*mu_inverse``.  Their
    mechanical wave speed is ``(b/a)*mu_inverse/epsilon`` and equals the
    electromagnetic value only when the positive conversion factors agree.
    With dimensionless strain ``xi``, the mechanical energy density is
    ``K*xi**2/2`` and its mass equivalent is ``a*epsilon*xi**2/2``.  The
    function derives consequences of supplied conversions; it does not infer
    a material, calibration, strain, or field amplitude from SI constants.
    """

    epsilon = _exact_positive(permittivity, "permittivity")
    inverse_mu = _exact_positive(
        inverse_permeability, "inverse_permeability"
    )
    inertia_factor = _exact_positive(
        inertia_conversion, "inertia_conversion"
    )
    stiffness_factor = (
        inertia_factor
        if stiffness_conversion is None
        else _exact_positive(stiffness_conversion, "stiffness_conversion")
    )
    strain = _exact_real(strain_amplitude, "strain_amplitude")
    density = sp.simplify(inertia_factor * epsilon)
    stiffness = sp.simplify(stiffness_factor * inverse_mu)
    electromagnetic_speed_squared = sp.simplify(inverse_mu / epsilon)
    mechanical_speed_squared = sp.simplify(stiffness / density)
    energy_density = sp.simplify(stiffness * strain**2 / 2)
    return MechanicalMediumConversion(
        permittivity=epsilon,
        inverse_permeability=inverse_mu,
        inertia_conversion=inertia_factor,
        stiffness_conversion=stiffness_factor,
        strain_amplitude=strain,
        mass_density=density,
        stiffness=stiffness,
        electromagnetic_speed_squared=electromagnetic_speed_squared,
        mechanical_speed_squared=mechanical_speed_squared,
        speed_squared_ratio=sp.simplify(stiffness_factor / inertia_factor),
        strain_energy_density=energy_density,
        mass_equivalent_density=sp.simplify(
            energy_density / mechanical_speed_squared
        ),
    )


def co_scaled_permittivity(
    density: Any,
    thermal_scale: Any,
    reference_speed: Any,
) -> sp.Expr:
    """Return the declared response ``epsilon = rho*Theta/c**2``."""

    rho = _positive(density, "density")
    theta = _positive(thermal_scale, "thermal_scale")
    speed = _positive(reference_speed, "reference_speed")
    return rho * theta / speed**2


def co_scaled_inverse_permeability(
    density: Any,
    thermal_scale: Any,
) -> sp.Expr:
    """Return the declared response ``mu**-1 = rho*Theta``."""

    rho = _positive(density, "density")
    theta = _positive(thermal_scale, "thermal_scale")
    return rho * theta


def local_wave_speed(
    permittivity: Any,
    inverse_permeability: Any,
) -> sp.Expr:
    """Return ``1/sqrt(epsilon*mu) = sqrt(mu_inverse/epsilon)``."""

    epsilon = _positive(permittivity, "permittivity")
    inverse_mu = _positive(inverse_permeability, "inverse_permeability")
    return sp.sqrt(inverse_mu / epsilon)


def co_scaled_wave_speed(
    density: Any,
    thermal_scale: Any,
    reference_speed: Any,
) -> sp.Expr:
    """Return the wave speed of the declared co-scaled responses."""

    return sp.simplify(
        local_wave_speed(
            co_scaled_permittivity(
                density, thermal_scale, reference_speed
            ),
            co_scaled_inverse_permeability(density, thermal_scale),
        )
    )


def lattice_debye_energy(
    action_scale: Any,
    reference_speed: Any,
    length_scale: Any,
    speed_ratio: Any = 1,
) -> sp.Expr:
    """Return the declared lattice scale ``Theta=kappa*S*c/a``.

    This is a conditional Debye/zero-point premise. Dimensional analysis makes
    the monomial available but does not select it dynamically or determine the
    dimensionless speed ratio ``kappa``.
    """

    action = _positive(action_scale, "action_scale")
    speed = _positive(reference_speed, "reference_speed")
    length = _positive(length_scale, "length_scale")
    ratio = _positive(speed_ratio, "speed_ratio")
    return sp.simplify(ratio * action * speed / length)


def lattice_reduced_responses(
    action_scale: Any,
    reference_speed: Any,
    length_scale: Any,
    speed_ratio: Any = 1,
    mass_density_ratio: Any = sp.Rational(1, 2),
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Compose declared lattice and Debye premises with co-scaled responses.

    Returns ``(Theta, epsilon, mu_inverse, mass_density)`` conditional on
    ``n=a**-3``, ``Theta=kappa*S*c/a``, the co-scaled response laws, and
    ``mass_density=mass_density_ratio*epsilon``.
    """

    action = _positive(action_scale, "action_scale")
    speed = _positive(reference_speed, "reference_speed")
    length = _positive(length_scale, "length_scale")
    ratio = _positive(speed_ratio, "speed_ratio")
    density_ratio = _positive(mass_density_ratio, "mass_density_ratio")
    number_density = length**-3
    thermal_scale = lattice_debye_energy(action, speed, length, ratio)
    epsilon = co_scaled_permittivity(number_density, thermal_scale, speed)
    inverse_mu = co_scaled_inverse_permeability(number_density, thermal_scale)
    return (
        thermal_scale,
        sp.simplify(epsilon),
        sp.simplify(inverse_mu),
        sp.simplify(density_ratio * epsilon),
    )
