"""Static Einstein weak field of one exact sine-Gordon breather energy."""

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
from .sine_gordon_weave import SineGordonWeaveGravity, sine_gordon_weave_gravity


@dataclass(frozen=True)
class BreatherWeakFieldPrediction:
    """Exact upstream source data and its leading static Einstein field."""

    frequency: sp.Expr
    inverse_width: sp.Expr
    observables: DimensionalBreatherObservables
    weave_gravity: SineGordonWeaveGravity
    mass: sp.Expr
    profile_length: sp.Expr
    schwarzschild_radius: sp.Expr
    profile_compactness: sp.Expr
    field: WeakFieldMonopole


def breather_weak_field_prediction(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    radius: Any,
) -> BreatherWeakFieldPrediction:
    """Return the leading static field of the exact rest-breather energy."""

    exact_frequency = sp.sympify(frequency)
    observables = dimensional_breather_observables(exact_frequency, coefficients)
    inverse_width = breather_inverse_width(exact_frequency)
    weave = sine_gordon_weave_gravity(coefficients)
    mass = sp.simplify(observables.energy / weave.scales.signal_speed**2)
    profile_length = observables.profile_length
    field = weak_field_monopole(
        weave.horizon_equilibrium.newton_constant,
        weave.scales.signal_speed,
        mass,
        radius,
    )
    schwarzschild_radius = sp.simplify(
        2 * weave.horizon_equilibrium.newton_constant
        * mass
        / weave.scales.signal_speed**2
    )
    return BreatherWeakFieldPrediction(
        frequency=exact_frequency,
        inverse_width=inverse_width,
        observables=observables,
        weave_gravity=weave,
        mass=mass,
        profile_length=profile_length,
        schwarzschild_radius=schwarzschild_radius,
        profile_compactness=sp.simplify(schwarzschild_radius / profile_length),
        field=field,
    )
