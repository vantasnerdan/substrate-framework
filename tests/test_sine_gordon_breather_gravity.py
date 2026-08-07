from __future__ import annotations

import sympy as sp

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import (
    dimensional_sine_gordon_coefficients,
    rescale_dimensional_sine_gordon_coefficients,
)
from substrate_framework.sine_gordon_breather_gravity import (
    TeleparallelBreatherWeakFieldPrediction,
    teleparallel_breather_weak_field_prediction,
)


def test_public_package_exports_derived_breather_static_api() -> None:
    assert (
        framework.TeleparallelBreatherWeakFieldPrediction
        is TeleparallelBreatherWeakFieldPrediction
    )
    assert (
        framework.teleparallel_breather_weak_field_prediction
        is teleparallel_breather_weak_field_prediction
    )


def test_exact_breather_energy_and_derived_give_a_controlled_weak_field() -> None:
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    inverse_width = sp.Rational(3, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    profile_length = 1 / inverse_width
    prediction = teleparallel_breather_weak_field_prediction(
        frequency,
        coefficients,
        profile_length,
    )
    assert prediction.observables.energy == 16 * inverse_width
    assert prediction.source_mass == 16 * inverse_width
    assert prediction.profile_length == profile_length
    assert prediction.gravity.newton_constant == 1 / (8 * sp.pi)
    assert prediction.schwarzschild_radius == inverse_width * 4 / sp.pi
    assert prediction.profile_compactness == 4 * inverse_width**2 / sp.pi
    assert prediction.field.potential_over_speed_squared == (
        -2 * inverse_width**2 / sp.pi
    )
    assert float(prediction.profile_compactness) < 0.0012


def test_common_sg_multiplier_changes_source_and_g_reciprocally() -> None:
    inverse_width = sp.Rational(1, 20)
    frequency = sp.sqrt(1 - inverse_width**2)
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    baseline = teleparallel_breather_weak_field_prediction(
        frequency,
        coefficients,
        1 / inverse_width,
    )
    scaled = teleparallel_breather_weak_field_prediction(
        frequency,
        rescale_dimensional_sine_gordon_coefficients(coefficients, 7),
        1 / inverse_width,
    )
    assert scaled.source_mass == 7 * baseline.source_mass
    assert scaled.gravity.newton_constant == baseline.gravity.newton_constant / 7
    assert scaled.schwarzschild_radius == baseline.schwarzschild_radius
    assert scaled.field.newtonian_potential == baseline.field.newtonian_potential
    assert scaled.field.metric_00 == baseline.field.metric_00
    assert scaled.field.spatial_metric_factor == baseline.field.spatial_metric_factor
