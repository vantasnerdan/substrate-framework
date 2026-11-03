"""Exact passive scattering at a declared right-half-line boundary.

For a massless scalar on ``x >= 0``, the boundary condition
``phi_t - zeta*phi_x = 0`` with positive speed ``zeta`` removes bulk energy.
This module derives its harmonic reflection ledger.  It does not derive that
boundary law from an action or identify its channels with chirality, charge,
detector outcomes, weak interactions, or material absorption.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .branching import two_channel_allocation


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return sp.simplify(expression)


@dataclass(frozen=True)
class PassiveHalfLineScatteringLedger:
    """Exact amplitude, power, energy-rate, and reciprocity data.

    ``reference_contrast`` additionally assumes a separate channel with unit
    reflected power.  It is the difference of the two normalized shares from
    :func:`substrate_framework.branching.two_channel_allocation`, not an
    independently derived physical observable.
    """

    wave_speed: sp.Expr
    boundary_speed: sp.Expr
    normalized_impedance: sp.Expr
    amplitude_reflection: sp.Expr
    reflected_power_fraction: sp.Expr
    absorbed_power_fraction: sp.Expr
    bulk_energy_rate_per_spatial_trace_squared: sp.Expr
    reciprocal_boundary_speed: sp.Expr
    reciprocal_amplitude_reflection: sp.Expr
    reciprocal_reflected_power_fraction: sp.Expr
    reciprocal_absorbed_power_fraction: sp.Expr
    reference_contrast: sp.Expr
    reciprocal_reference_contrast: sp.Expr
    contrast_as_absorbed_transform_residual: sp.Expr


def passive_half_line_scattering_ledger(
    wave_speed: Any,
    boundary_speed: Any,
) -> PassiveHalfLineScatteringLedger:
    """Return the exact passive harmonic ledger on ``x >= 0``.

    The declared bulk equation is ``phi_tt-c**2*phi_xx=0``.  An incoming wave
    is proportional to ``exp(-i*omega*(t+x/c))`` and the reflected wave to
    ``exp(-i*omega*(t-x/c))``.  Substitution into
    ``phi_t-zeta*phi_x=0`` gives ``r=(z-1)/(z+1)`` for ``z=zeta/c``.

    With the canonical scalar energy, the right-half-line bulk rate contributed
    by the boundary is ``-c**2*zeta*phi_x(0,t)**2``.  Therefore positive
    ``zeta`` is passive.  The returned absorbed fraction assumes a steady
    harmonic experiment with no additional boundary storage or flux channel.
    """

    speed = _positive_exact(wave_speed, "wave_speed")
    boundary = _positive_exact(boundary_speed, "boundary_speed")
    impedance = sp.simplify(boundary / speed)
    reflection = sp.factor((impedance - 1) / (impedance + 1))
    reflected_power = sp.factor(reflection**2)
    absorbed_power = sp.factor(1 - reflected_power)

    reciprocal_boundary = sp.simplify(speed**2 / boundary)
    reciprocal_impedance = sp.simplify(reciprocal_boundary / speed)
    reciprocal_reflection = sp.factor(
        (reciprocal_impedance - 1) / (reciprocal_impedance + 1)
    )
    reciprocal_reflected_power = sp.factor(reciprocal_reflection**2)
    reciprocal_absorbed_power = sp.factor(1 - reciprocal_reflected_power)

    allocation = two_channel_allocation(sp.S.One, reflected_power)
    contrast = sp.factor(
        allocation.first_fraction - allocation.second_fraction
    )
    reciprocal_allocation = two_channel_allocation(
        sp.S.One,
        reciprocal_reflected_power,
    )
    reciprocal_contrast = sp.factor(
        reciprocal_allocation.first_fraction
        - reciprocal_allocation.second_fraction
    )

    return PassiveHalfLineScatteringLedger(
        wave_speed=speed,
        boundary_speed=boundary,
        normalized_impedance=impedance,
        amplitude_reflection=reflection,
        reflected_power_fraction=reflected_power,
        absorbed_power_fraction=absorbed_power,
        bulk_energy_rate_per_spatial_trace_squared=sp.factor(
            -speed**2 * boundary
        ),
        reciprocal_boundary_speed=reciprocal_boundary,
        reciprocal_amplitude_reflection=reciprocal_reflection,
        reciprocal_reflected_power_fraction=reciprocal_reflected_power,
        reciprocal_absorbed_power_fraction=reciprocal_absorbed_power,
        reference_contrast=contrast,
        reciprocal_reference_contrast=reciprocal_contrast,
        contrast_as_absorbed_transform_residual=sp.simplify(
            contrast - absorbed_power / (2 - absorbed_power)
        ),
    )
