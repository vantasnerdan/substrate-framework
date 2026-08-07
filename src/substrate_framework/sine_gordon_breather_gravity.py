"""Static weak field of one sine-Gordon breather with derived gravity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalBreatherObservables,
    DimensionalSineGordonCoefficients,
    dimensional_breather_observables,
)
from .linearized_einstein import WeakFieldMonopole, weak_field_monopole
from .sine_gordon import breather_inverse_width
from .sine_gordon_teleparallel import (
    SineGordonTeleparallelGravity,
    sine_gordon_teleparallel_gravity,
)


@dataclass(frozen=True)
class TeleparallelBreatherWeakFieldPrediction:
    """Exact breather energy and its derived-coupling static Einstein field."""

    frequency: sp.Expr
    inverse_width: sp.Expr
    observables: DimensionalBreatherObservables
    gravity: SineGordonTeleparallelGravity
    source_mass: sp.Expr
    profile_length: sp.Expr
    schwarzschild_radius: sp.Expr
    profile_compactness: sp.Expr
    field: WeakFieldMonopole


def teleparallel_breather_weak_field_prediction(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    radius: Any,
) -> TeleparallelBreatherWeakFieldPrediction:
    r"""Return the leading mostly-plus field of the exact rest breather.

    The source mass is not independent: ``M=E_breather/c**2``.  The Newton
    constant is not supplied: it is ``c**4/(8*pi*mu)`` from the SG
    teleparallel ledger.  Thus the exact compactness at the profile length is
    ``4*eta**2/pi`` for inverse width ``eta=sqrt(1-omega**2)``.
    """

    exact_frequency = sp.sympify(frequency)
    observables = dimensional_breather_observables(exact_frequency, coefficients)
    inverse_width = breather_inverse_width(exact_frequency)
    gravity = sine_gordon_teleparallel_gravity(coefficients)
    mass = sp.simplify(observables.energy / gravity.scales.signal_speed**2)
    field = weak_field_monopole(
        gravity.newton_constant,
        gravity.scales.signal_speed,
        mass,
        radius,
    )
    schwarzschild_radius = sp.simplify(
        2
        * gravity.newton_constant
        * mass
        / gravity.scales.signal_speed**2
    )
    compactness = sp.simplify(schwarzschild_radius / observables.profile_length)
    expected_compactness = sp.simplify(4 * inverse_width**2 / sp.pi)
    if sp.simplify(compactness - expected_compactness) != 0:
        raise AssertionError("breather compactness did not reduce upstream")
    return TeleparallelBreatherWeakFieldPrediction(
        frequency=exact_frequency,
        inverse_width=inverse_width,
        observables=observables,
        gravity=gravity,
        source_mass=mass,
        profile_length=observables.profile_length,
        schwarzschild_radius=schwarzschild_radius,
        profile_compactness=compactness,
        field=field,
    )
