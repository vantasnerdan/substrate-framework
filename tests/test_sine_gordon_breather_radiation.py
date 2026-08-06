from __future__ import annotations

import sympy as sp
import pytest

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import dimensional_sine_gordon_coefficients
from substrate_framework.sine_gordon_breather_radiation import (
    BreatherRadiationEvidence,
    breather_radiation_evidence,
    leading_small_amplitude_breather_power_ratio,
)


def test_public_package_exports_breather_radiation_api() -> None:
    assert framework.BreatherRadiationEvidence is BreatherRadiationEvidence
    assert framework.breather_radiation_evidence is breather_radiation_evidence
    assert (
        framework.leading_small_amplitude_breather_power_ratio
        is leading_small_amplitude_breather_power_ratio
    )


def test_exact_breather_stress_has_converged_nonzero_full_retardation_power() -> None:
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    frequency = sp.sqrt(9991) / 100
    coarse = breather_radiation_evidence(
        frequency,
        coefficients,
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    refined = breather_radiation_evidence(
        frequency,
        coefficients,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    assert abs(
        refined.power_over_onsite_speed - coarse.power_over_onsite_speed
    ) < 3.0e-10
    assert abs(
        refined.power_over_onsite_speed - 0.002606646892787329
    ) < 2.0e-15
    assert refined.power > 0.0
    assert refined.source_size_frequency_ratio > 30.0
    assert refined.relative_leading_difference < 7.3e-4
    assert refined.harmonic_pressure_maxima[1] < 3.3e-6
    assert refined.harmonic_pressure_maxima[2] < 1.1e-9
    assert refined.field_cycle_energy_loss_fraction < 0.035


def test_common_microscopic_scale_changes_absolute_power_not_dimensionless_shape() -> None:
    frequency = sp.sqrt(9991) / 100
    baseline = breather_radiation_evidence(
        frequency,
        dimensional_sine_gordon_coefficients(1, 1, 1),
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    scaled = breather_radiation_evidence(
        frequency,
        dimensional_sine_gordon_coefficients(5, 5, 5),
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    assert scaled.newton_constant == baseline.newton_constant / 5
    assert scaled.power_scale == 5 * baseline.power_scale
    assert scaled.power_over_onsite_speed == baseline.power_over_onsite_speed
    assert scaled.power == 5 * baseline.power


@pytest.mark.parametrize(
    ("keyword", "value", "message"),
    [
        ("spatial_extent_in_inverse_widths", -1.0, "spatial_extent"),
        ("spatial_points", 512, "spatial_points"),
        ("time_points", 32, "time_points"),
        ("direction_points", 128, "direction_points"),
        ("harmonic_count", 0, "harmonic_count"),
        ("harmonic_count", 9, "four samples"),
    ],
)
def test_radiation_resolution_inputs_fail_with_context(
    keyword: str,
    value: object,
    message: str,
) -> None:
    arguments = {
        "spatial_points": 257,
        "time_points": 33,
        "direction_points": 129,
        "harmonic_count": 2,
        keyword: value,
    }
    with pytest.raises(ValueError, match=message):
        breather_radiation_evidence(
            sp.Rational(4, 5),
            dimensional_sine_gordon_coefficients(1, 1, 1),
            **arguments,
        )


def test_radiation_requires_exact_frequency_and_numeric_coefficients() -> None:
    with pytest.raises(ValueError, match="frequency"):
        breather_radiation_evidence(
            0.8,
            dimensional_sine_gordon_coefficients(1, 1, 1),
        )
    scale = sp.symbols("scale", positive=True)
    with pytest.raises(ValueError, match="coefficients.onsite"):
        breather_radiation_evidence(
            sp.Rational(4, 5),
            dimensional_sine_gordon_coefficients(scale, scale, scale),
            spatial_points=257,
            time_points=33,
            direction_points=129,
            harmonic_count=2,
        )
