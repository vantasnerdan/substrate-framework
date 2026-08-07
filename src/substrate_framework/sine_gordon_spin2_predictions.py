"""Static and radiative predictions using the derived wall-network coupling.

The source is the accepted dimensional sine-Gordon breather in both cases.
The Newton constant is never an input: it is obtained from
``sine_gordon_wall_network_scales``.  The static result is stated with a
rigorous finite-profile error bound, while the radiation result reuses the
accepted full-retardation source transform.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_breather_observables,
)
from .exact_symbolic import positive_exact as _positive_exact
from .emergent_spin2 import sine_gordon_spin_two_coupling_ledger
from .linearized_einstein import WeakFieldMonopole, weak_field_monopole
from .sine_gordon import (
    breather_energy_second_moment_extrema,
    breather_inverse_width,
)
from .sine_gordon_breather_radiation import (
    SineGordonBreatherRadiationEvidence,
    sine_gordon_breather_radiation_evidence,
)


@dataclass(frozen=True)
class ControlledStaticBreatherPrediction:
    """Transverse static monopole with an exact profile-error enclosure.

    The exact line-source potential is
    ``-G/c**2 * integral dx epsilon(x)/sqrt(R**2+x**2)``.  Positivity of the
    breather energy density and
    ``0 <= 1-R/sqrt(R**2+x**2) <= x**2/(2*R**2)`` enclose it using the exact
    cycle-maximum second moment.
    """

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
    maximum_relative_profile_error: sp.Expr
    exact_potential_lower_bound: sp.Expr
    exact_potential_upper_bound: sp.Expr


def controlled_static_breather_prediction(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    observer_radius_in_profile_lengths: Any,
) -> ControlledStaticBreatherPrediction:
    """Return a derived-coupling static prediction with error below one.

    ``observer_radius_in_profile_lengths`` is the transverse radius divided by
    ``ell/eta``.  The function rejects a radius for which the moment enclosure
    does not establish a controlled (sub-unit) relative error.
    """

    coupling = sine_gordon_spin_two_coupling_ledger(coefficients)
    wall = coupling.wall_scales
    radius_multiple = _positive_exact(
        observer_radius_in_profile_lengths,
        "observer_radius_in_profile_lengths",
    )
    observables = dimensional_breather_observables(frequency, coefficients)
    exact_frequency = sp.sympify(frequency)
    eta = breather_inverse_width(exact_frequency)
    radius = sp.simplify(radius_multiple * wall.profile_length / eta)
    source_energy = observables.energy
    source_mass = sp.simplify(source_energy / wall.signal_speed**2)
    _, normalized_second_moment_maximum = breather_energy_second_moment_extrema(
        exact_frequency
    )
    physical_second_moment_maximum = sp.simplify(
        wall.energy_scale
        * wall.profile_length**2
        * normalized_second_moment_maximum
    )
    relative_error = sp.simplify(
        physical_second_moment_maximum / (2 * source_energy * radius**2)
    )
    _positive_exact(1 - relative_error, "one minus profile error bound")
    monopole = weak_field_monopole(
        coupling.newton_constant,
        wall.signal_speed,
        source_mass,
        radius,
    )
    lower = monopole.newtonian_potential
    upper = sp.simplify(lower * (1 - relative_error))
    if sp.simplify(
        relative_error
        - eta * normalized_second_moment_maximum / (32 * radius_multiple**2)
    ) != 0:
        raise AssertionError("physical static-profile scaling did not reduce")

    return ControlledStaticBreatherPrediction(
        frequency=exact_frequency,
        inverse_width=eta,
        observer_radius_in_profile_lengths=radius_multiple,
        observer_radius=radius,
        source_energy=source_energy,
        source_mass=source_mass,
        source_second_moment_maximum=physical_second_moment_maximum,
        newton_constant=coupling.newton_constant,
        signal_speed=wall.signal_speed,
        monopole=monopole,
        maximum_relative_profile_error=relative_error,
        exact_potential_lower_bound=lower,
        exact_potential_upper_bound=upper,
    )


def sine_gordon_spin2_radiation_evidence(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    *,
    spatial_extent_in_inverse_widths: float = 15.0,
    spatial_points: int = 4097,
    time_points: int = 257,
    direction_points: int = 2001,
    harmonic_count: int = 6,
) -> SineGordonBreatherRadiationEvidence:
    """Compose the accepted retarded source transform with the derived ``G``."""

    coupling = sine_gordon_spin_two_coupling_ledger(coefficients)
    wall = coupling.wall_scales
    result = sine_gordon_breather_radiation_evidence(
        frequency,
        coefficients,
        coupling.newton_constant,
        spatial_extent_in_inverse_widths=spatial_extent_in_inverse_widths,
        spatial_points=spatial_points,
        time_points=time_points,
        direction_points=direction_points,
        harmonic_count=harmonic_count,
    )
    expected_dimensionless_coupling = sp.Rational(5, 1024) / sp.pi
    if sp.simplify(
        result.newton_constant
        * coefficients.onsite
        / result.signal_speed**4
        - expected_dimensionless_coupling
    ) != 0:
        raise AssertionError("derived radiation coupling retained a hidden scale")
    return result
