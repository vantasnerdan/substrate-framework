"""Exact kinematics and conditional TT observables for a circular point pair.

The point paths are declared inputs.  These APIs do not provide a binding
stress, an orbital equation, a breather embedding, or a physical gravitational
field law.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .conserved_moments import DiscreteMassMoments, discrete_mass_moments
from .tt_angular import (
    TTPolarizationBasis,
    conditional_tt_power,
    frobenius_inner_product,
    tt_polarization_basis,
    tt_project_symmetric,
)


def _time_symbol(time: Any) -> sp.Symbol:
    if not isinstance(time, sp.Symbol):
        raise ValueError("time must be a SymPy Symbol")
    return time


def equal_mass_circular_pair_moments(
    mass: Any,
    orbital_radius: Any,
    angular_frequency: Any,
    time: Any,
) -> DiscreteMassMoments:
    """Return exact moments for masses at opposite points of a circle.

    ``orbital_radius`` is each mass's distance from the origin, so the pair
    separation is twice this value.
    """

    time_symbol = _time_symbol(time)
    mass_value = sp.sympify(mass)
    radius = sp.sympify(orbital_radius)
    frequency = sp.sympify(angular_frequency)
    position = sp.Matrix(
        [
            radius * sp.cos(frequency * time_symbol),
            radius * sp.sin(frequency * time_symbol),
            0,
        ]
    )
    return discrete_mass_moments(
        [mass_value, mass_value],
        [position, -position],
    )


@dataclass(frozen=True)
class ConditionalCircularPairWaveform:
    """Circular-pair moments and coefficients under a declared TT waveform."""

    moments: DiscreteMassMoments
    normalized_stf_second_derivative: sp.Matrix
    normalized_stf_third_derivative: sp.Matrix
    polarization_basis: TTPolarizationBasis
    tt_tensor: sp.Matrix
    normalized_plus_coordinate: sp.Expr
    normalized_cross_coordinate: sp.Expr
    conventional_plus: sp.Expr
    conventional_cross: sp.Expr


def conditional_equal_mass_circular_waveform(
    mass: Any,
    orbital_radius: Any,
    angular_frequency: Any,
    time: Any,
    inclination: Any,
    waveform_prefactor: Any,
    distance: Any,
) -> ConditionalCircularPairWaveform:
    """Return exact conditional TT coefficients for a circular point pair.

    The orbit lies in the x-y plane.  ``inclination`` is measured from the
    positive z axis, with line of sight ``(sin i, 0, cos i)`` and transverse
    frame ``p=(cos i,0,-sin i)``, ``v=(0,1,0)``.  The declared field premise is
    ``h_TT=(waveform_prefactor/distance)*TT(I_STF_ddot)``.

    Normalized basis coordinates use C-GW-002.  The conventional matrix
    read-offs are ``h_plus=(h_pp-h_vv)/2`` and ``h_cross=h_pv`` and therefore
    equal the normalized coordinates divided by ``sqrt(2)``.
    """

    time_symbol = _time_symbol(time)
    angle = sp.sympify(inclination)
    wave = sp.sympify(waveform_prefactor)
    radius_to_observer = sp.sympify(distance)
    if radius_to_observer == 0:
        raise ValueError("distance must be nonzero")
    moments = equal_mass_circular_pair_moments(
        mass,
        orbital_radius,
        angular_frequency,
        time_symbol,
    )
    second = sp.simplify(
        moments.trace_free_second_moment.diff(time_symbol, 2)
    )
    third = sp.simplify(
        moments.trace_free_second_moment.diff(time_symbol, 3)
    )
    line_of_sight = sp.Matrix([sp.sin(angle), 0, sp.cos(angle)])
    transverse_reference = sp.Matrix([sp.cos(angle), 0, -sp.sin(angle)])
    basis = tt_polarization_basis(line_of_sight, transverse_reference)
    tt_tensor = sp.simplify(
        wave
        / radius_to_observer
        * tt_project_symmetric(second, line_of_sight)
    )
    plus = sp.simplify(frobenius_inner_product(tt_tensor, basis.plus))
    cross = sp.simplify(frobenius_inner_product(tt_tensor, basis.cross))
    return ConditionalCircularPairWaveform(
        moments=moments,
        normalized_stf_second_derivative=second,
        normalized_stf_third_derivative=third,
        polarization_basis=basis,
        tt_tensor=tt_tensor,
        normalized_plus_coordinate=plus,
        normalized_cross_coordinate=cross,
        conventional_plus=sp.simplify(plus / sp.sqrt(2)),
        conventional_cross=sp.simplify(cross / sp.sqrt(2)),
    )


def conditional_equal_mass_circular_power(
    mass: Any,
    orbital_radius: Any,
    angular_frequency: Any,
    waveform_prefactor: Any,
    flux_prefactor: Any,
) -> sp.Expr:
    """Return the constant angular power under declared waveform/flux inputs."""

    time = sp.Symbol("_circular_pair_time", real=True)
    moments = equal_mass_circular_pair_moments(
        mass,
        orbital_radius,
        angular_frequency,
        time,
    )
    third = sp.simplify(moments.trace_free_second_moment.diff(time, 3))
    return sp.trigsimp(
        conditional_tt_power(
            third,
            waveform_prefactor=waveform_prefactor,
            flux_prefactor=flux_prefactor,
        )
    )
