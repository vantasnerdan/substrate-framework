from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import (
    dimensional_breather_field,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_stress,
)
from substrate_framework.sine_gordon_breather_radiation import (
    TeleparallelBreatherRadiationEvidence,
    _normalized_breather_longitudinal_stress,
    leading_small_amplitude_teleparallel_breather_power_ratio,
    teleparallel_breather_radiation_evidence,
)


def test_public_package_exports_full_retardation_breather_api() -> None:
    assert (
        framework.TeleparallelBreatherRadiationEvidence
        is TeleparallelBreatherRadiationEvidence
    )
    assert (
        framework.leading_small_amplitude_teleparallel_breather_power_ratio
        is leading_small_amplitude_teleparallel_breather_power_ratio
    )
    assert (
        framework.teleparallel_breather_radiation_evidence
        is teleparallel_breather_radiation_evidence
    )


def test_numeric_pressure_kernel_is_the_exact_conserved_sg_stress_component() -> None:
    x, t = sp.symbols("x t", real=True)
    frequency = sp.Rational(4, 5)
    inverse_width = sp.Rational(3, 5)
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    field = dimensional_breather_field(x, t, frequency, coefficients)
    symbolic_pressure = dimensional_sine_gordon_stress(
        field,
        x,
        t,
        coefficients,
    ).longitudinal_pressure
    evaluator = sp.lambdify((x, t), symbolic_pressure, "numpy")
    coordinates = np.array([-1.25, -0.2, 0.0, 0.7])
    phases = np.array([0.0, 0.3, 0.9, np.pi])
    direct = evaluator(coordinates[:, None], phases[None, :] / float(frequency))
    kernel = _normalized_breather_longitudinal_stress(
        float(inverse_width) * coordinates,
        phases,
        float(frequency),
        float(inverse_width),
    )
    np.testing.assert_allclose(kernel, direct, rtol=3.0e-14, atol=3.0e-14)


def test_exact_stress_has_converged_nonzero_full_retardation_power() -> None:
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    frequency = sp.sqrt(9991) / 100
    coarse = teleparallel_breather_radiation_evidence(
        frequency,
        coefficients,
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    refined = teleparallel_breather_radiation_evidence(
        frequency,
        coefficients,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    assert abs(
        refined.power_over_onsite_speed - coarse.power_over_onsite_speed
    ) < 4.0e-11
    assert abs(
        refined.power_over_onsite_speed - 0.0002875595507880884
    ) < 3.0e-17
    assert refined.power > 0.0
    assert refined.source_size_frequency_ratio > 30.0
    assert refined.relative_leading_difference < 7.3e-4
    assert refined.harmonic_pressure_maxima[1] < 3.3e-6
    assert refined.harmonic_pressure_maxima[2] < 1.1e-9
    assert refined.field_cycle_energy_loss_fraction < 0.004


def test_common_microscopic_scale_changes_absolute_power_not_shape() -> None:
    frequency = sp.sqrt(9991) / 100
    baseline = teleparallel_breather_radiation_evidence(
        frequency,
        dimensional_sine_gordon_coefficients(1, 1, 1),
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    scaled = teleparallel_breather_radiation_evidence(
        frequency,
        dimensional_sine_gordon_coefficients(5, 5, 5),
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    assert scaled.gravity.newton_constant == baseline.gravity.newton_constant / 5
    assert scaled.power_scale == 5 * baseline.power_scale
    assert scaled.power_over_onsite_speed == pytest.approx(
        baseline.power_over_onsite_speed,
        rel=3.0e-16,
    )
    assert scaled.power == pytest.approx(5 * baseline.power, rel=3.0e-16)


def test_small_width_transform_has_the_derived_analytic_asymptote() -> None:
    inverse_width = sp.Rational(1, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    exact_transform = leading_small_amplitude_teleparallel_breather_power_ratio(
        frequency,
        direction_points=4001,
    )
    asymptotic = 32.0 * float(inverse_width) ** 3 * float(frequency) / 3.0
    assert abs(exact_transform / asymptotic - 1.0) < 4.1e-5


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
        teleparallel_breather_radiation_evidence(
            sp.Rational(4, 5),
            dimensional_sine_gordon_coefficients(1, 1, 1),
            **arguments,
        )


def test_radiation_requires_exact_frequency_and_numeric_coefficients() -> None:
    with pytest.raises(ValueError, match="frequency"):
        teleparallel_breather_radiation_evidence(
            0.8,
            dimensional_sine_gordon_coefficients(1, 1, 1),
        )
    scale = sp.symbols("scale", positive=True)
    with pytest.raises(ValueError, match="coefficients.onsite"):
        teleparallel_breather_radiation_evidence(
            sp.Rational(4, 5),
            dimensional_sine_gordon_coefficients(scale, scale, scale),
            spatial_points=257,
            time_points=33,
            direction_points=129,
            harmonic_count=2,
        )
