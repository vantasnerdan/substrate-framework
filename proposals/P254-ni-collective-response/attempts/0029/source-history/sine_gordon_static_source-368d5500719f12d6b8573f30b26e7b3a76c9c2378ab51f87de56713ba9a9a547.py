"""Supplied-coupling Poisson enclosure for a sine-Gordon breather fiber.

The source is the accepted dimensional sine-Gordon breather embedded along a
line.  The Newton constant is an explicit positive input.  Positivity and the
exact breather second moment bound the difference between the extended-source
potential at transverse radius ``R`` and its line-profile monopole limit.
The returned compactness exposes, but does not assume or prove, the weak-field
regime.  Retardation, stationarity, and a physical three-dimensional source
map remain separate obligations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_breather_observables,
    dimensional_sine_gordon_scales,
)
from .exact_symbolic import positive_exact as _positive_exact
from .linearized_einstein import WeakFieldMonopole, weak_field_monopole
from .sine_gordon import (
    breather_energy_second_moment_extrema,
    breather_inverse_width,
)


@dataclass(frozen=True)
class SineGordonStaticSourceEnclosure:
    """Transverse line-source potential with a second-moment error bound."""

    frequency: sp.Expr
    inverse_width: sp.Expr
    observer_radius_in_profile_lengths: sp.Expr
    observer_radius: sp.Expr
    source_energy: sp.Expr
    source_mass: sp.Expr
    source_second_moment_maximum: sp.Expr
    newton_constant: sp.Expr
    signal_speed: sp.Expr
    monopole: WeakFieldMonopole
    monopole_compactness: sp.Expr
    maximum_relative_profile_error: sp.Expr
    exact_potential_lower_bound: sp.Expr
    exact_potential_upper_bound: sp.Expr


def sine_gordon_static_source_enclosure(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    newton_constant: Any,
    observer_radius_in_profile_lengths: Any,
) -> SineGordonStaticSourceEnclosure:
    r"""Return a supplied-``G`` line-profile Poisson enclosure.

    The exact conditional potential is
    ``-G/c**2*integral epsilon(x)/sqrt(R**2+x**2) dx``.  For positive energy
    density,
    ``0 <= 1-R/sqrt(R**2+x**2) <= x**2/(2*R**2)`` bounds its relative
    difference from ``-G*M/R`` by the exact cycle-maximum second moment.  A
    radius whose profile bound is not below one is rejected.  The returned
    ``monopole_compactness=G*M/(c**2*R)`` lets the caller impose a separately
    justified weak-field threshold.
    """

    newton = _positive_exact(newton_constant, "newton_constant")
    radius_multiple = _positive_exact(
        observer_radius_in_profile_lengths,
        "observer_radius_in_profile_lengths",
    )
    scales = dimensional_sine_gordon_scales(coefficients)
    observables = dimensional_breather_observables(frequency, coefficients)
    exact_frequency = sp.sympify(frequency)
    eta = breather_inverse_width(exact_frequency)
    radius = sp.simplify(radius_multiple * scales.length / eta)
    source_energy = observables.energy
    source_mass = sp.simplify(source_energy / scales.signal_speed**2)
    _, normalized_second_moment_maximum = breather_energy_second_moment_extrema(
        exact_frequency
    )
    physical_second_moment_maximum = sp.simplify(
        scales.energy * scales.length**2 * normalized_second_moment_maximum
    )
    relative_error = sp.simplify(
        physical_second_moment_maximum / (2 * source_energy * radius**2)
    )
    _positive_exact(1 - relative_error, "one minus profile error bound")
    monopole = weak_field_monopole(
        newton,
        scales.signal_speed,
        source_mass,
        radius,
    )
    lower = monopole.newtonian_potential
    upper = sp.simplify(lower * (1 - relative_error))
    compactness = sp.simplify(-monopole.potential_over_speed_squared)
    if sp.simplify(
        relative_error
        - eta * normalized_second_moment_maximum / (32 * radius_multiple**2)
    ) != 0:
        raise AssertionError("physical static-profile scaling did not reduce")

    return SineGordonStaticSourceEnclosure(
        frequency=exact_frequency,
        inverse_width=eta,
        observer_radius_in_profile_lengths=radius_multiple,
        observer_radius=radius,
        source_energy=source_energy,
        source_mass=source_mass,
        source_second_moment_maximum=physical_second_moment_maximum,
        newton_constant=newton,
        signal_speed=scales.signal_speed,
        monopole=monopole,
        monopole_compactness=compactness,
        maximum_relative_profile_error=relative_error,
        exact_potential_lower_bound=lower,
        exact_potential_upper_bound=upper,
    )
