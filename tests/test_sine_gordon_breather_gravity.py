from __future__ import annotations

import sympy as sp

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import dimensional_sine_gordon_coefficients
from substrate_framework.sine_gordon_breather_gravity import (
    BreatherWeakFieldPrediction,
    breather_weak_field_prediction,
)


def test_public_package_exports_breather_gravity_api() -> None:
    assert framework.BreatherWeakFieldPrediction is BreatherWeakFieldPrediction
    assert framework.breather_weak_field_prediction is breather_weak_field_prediction


def test_same_small_breather_has_a_controlled_static_weak_field() -> None:
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    inverse_width = sp.Rational(3, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    profile_length = 1 / inverse_width
    prediction = breather_weak_field_prediction(
        frequency,
        coefficients,
        profile_length,
    )
    assert prediction.observables.energy == 16 * inverse_width
    assert prediction.mass == 16 * inverse_width
    assert prediction.profile_length == profile_length
    assert prediction.schwarzschild_radius == 8 * inverse_width / sp.log(2)
    assert prediction.profile_compactness == 8 * inverse_width**2 / sp.log(2)
    assert prediction.field.potential_over_speed_squared == (
        -4 * inverse_width**2 / sp.log(2)
    )
    assert float(prediction.profile_compactness) < 0.011
