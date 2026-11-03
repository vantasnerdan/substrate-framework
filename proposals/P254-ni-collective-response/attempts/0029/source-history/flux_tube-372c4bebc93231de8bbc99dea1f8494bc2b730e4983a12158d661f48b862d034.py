"""Exact conditional algebra for a uniform fixed-area flux tube."""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def uniform_tube_field(flux: Any, cross_section: Any) -> sp.Expr:
    """Return the field ``flux/cross_section`` under uniform-cap Gauss data."""

    flux_value = _positive(flux, "flux")
    area_value = _positive(cross_section, "cross_section")
    return sp.simplify(flux_value / area_value)


def tube_energy_slope(flux: Any, cross_section: Any) -> sp.Expr:
    """Return field energy per length for density ``field**2/2``."""

    field = uniform_tube_field(flux, cross_section)
    area_value = _positive(cross_section, "cross_section")
    return sp.simplify(field**2 * area_value / 2)


def tube_field_energy(length: Any, flux: Any, cross_section: Any) -> sp.Expr:
    """Return stored field energy for a uniform fixed-area tube segment."""

    length_value = _positive(length, "length")
    return sp.simplify(length_value * tube_energy_slope(flux, cross_section))


def endpoint_force_slope(charge: Any, flux: Any, cross_section: Any) -> sp.Expr:
    """Return the constant endpoint-force magnitude ``charge*field``."""

    charge_value = _positive(charge, "charge")
    return sp.simplify(
        charge_value * uniform_tube_field(flux, cross_section)
    )


def endpoint_potential(
    length: Any, charge: Any, flux: Any, cross_section: Any
) -> sp.Expr:
    """Return endpoint work from the declared constant-force law."""

    length_value = _positive(length, "length")
    return sp.simplify(
        length_value * endpoint_force_slope(charge, flux, cross_section)
    )


def charge_for_slope_equality(flux: Any) -> sp.Expr:
    """Return the charge required to equate endpoint and field-energy slopes."""

    return sp.simplify(_positive(flux, "flux") / 2)


def spherical_field(flux: Any, radius: Any) -> sp.Expr:
    """Return the uniform radial field on a sphere of the declared radius."""

    flux_value = _positive(flux, "flux")
    radius_value = _positive(radius, "radius")
    return sp.simplify(flux_value / (4 * sp.pi * radius_value**2))
