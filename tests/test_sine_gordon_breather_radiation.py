from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import simpson
import sympy as sp

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import (
    dimensional_breather_field,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_stress,
)
from substrate_framework.sine_gordon_breather_radiation import (
    SineGordonBreatherRadiationEvidence,
    SineGordonBreatherRetardedStressEvidence,
    _normalized_breather_longitudinal_stress,
    leading_small_amplitude_breather_retarded_integral,
    sine_gordon_breather_radiation_evidence,
    sine_gordon_breather_retarded_stress_evidence,
)


def test_public_package_exports_full_retardation_breather_api() -> None:
    assert (
        framework.SineGordonBreatherRetardedStressEvidence
        is SineGordonBreatherRetardedStressEvidence
    )
    assert (
        framework.SineGordonBreatherRadiationEvidence
        is SineGordonBreatherRadiationEvidence
    )
    assert (
        framework.leading_small_amplitude_breather_retarded_integral
        is leading_small_amplitude_breather_retarded_integral
    )
    assert (
        framework.sine_gordon_breather_retarded_stress_evidence
        is sine_gordon_breather_retarded_stress_evidence
    )
    assert (
        framework.sine_gordon_breather_radiation_evidence
        is sine_gordon_breather_radiation_evidence
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


def test_source_transform_converges_without_a_gravitational_coupling() -> None:
    frequency = sp.sqrt(9991) / 100
    coarse = sine_gordon_breather_retarded_stress_evidence(
        frequency,
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    refined = sine_gordon_breather_retarded_stress_evidence(
        frequency,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    expected = 16.0 * np.pi * 0.0002875595507880884
    assert abs(
        refined.dimensionless_retarded_integral
        - coarse.dimensionless_retarded_integral
    ) < 16.0 * np.pi * 4.0e-11
    assert refined.dimensionless_retarded_integral == pytest.approx(
        expected,
        abs=2.0e-15,
    )
    assert refined.radiation_frequency_size_ratio > 60.0
    assert refined.relative_leading_difference < 7.3e-4
    assert refined.harmonic_pressure_maxima[1] < 3.3e-6
    assert refined.harmonic_pressure_maxima[2] < 1.1e-9


def test_independent_simpson_transform_matches_the_package_integral() -> None:
    inverse_width = 0.03
    frequency = np.sqrt(1.0 - inverse_width**2)
    scaled_coordinate = np.linspace(-15.0, 15.0, 1025)
    coordinate = scaled_coordinate / inverse_width
    phase = np.linspace(0.0, np.pi, 129)
    sech = 1.0 / np.cosh(scaled_coordinate)[:, None]
    sine = np.sin(phase)[None, :]
    cosine = np.cos(phase)[None, :]
    tangent = np.tanh(scaled_coordinate)[:, None]
    argument = (inverse_width / frequency) * sech * sine
    denominator = 1.0 + argument**2
    field_time = 4.0 * inverse_width * sech * cosine / denominator
    field_space = -4.0 * inverse_width * argument * tangent / denominator
    potential = 8.0 * argument**2 / denominator**2
    pressure = 0.5 * (field_time**2 + field_space**2) - potential

    modes = np.arange(1, 5, dtype=np.float64)
    pressure_modes = np.asarray(
        [
            2.0
            / np.pi
            * simpson(
                pressure * np.cos(2.0 * mode * phase)[None, :],
                x=phase,
                axis=1,
            )
            for mode in modes
        ]
    )
    directions = np.linspace(0.0, 1.0, 501)
    mean_derivative_squared = np.zeros_like(directions)
    for mode, pressure_mode in zip(modes, pressure_modes):
        wavenumber = 2.0 * mode * frequency
        amplitudes = np.asarray(
            [
                simpson(
                    pressure_mode
                    * np.cos(wavenumber * direction * coordinate),
                    x=coordinate,
                )
                for direction in directions
            ]
        )
        mean_derivative_squared += 0.5 * (wavenumber * amplitudes) ** 2
    independent_integral = 2.0 * simpson(
        (1.0 - directions**2) ** 2 * mean_derivative_squared,
        x=directions,
    )
    package = sine_gordon_breather_retarded_stress_evidence(
        sp.sqrt(9991) / 100,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    assert independent_integral == pytest.approx(
        package.dimensionless_retarded_integral,
        abs=1.2e-10,
    )


def test_supplied_newton_constant_is_visible_and_load_bearing() -> None:
    frequency = sp.sqrt(9991) / 100
    coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    baseline = sine_gordon_breather_radiation_evidence(
        frequency,
        coefficients,
        1 / (8 * sp.pi),
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    doubled = sine_gordon_breather_radiation_evidence(
        frequency,
        coefficients,
        1 / (4 * sp.pi),
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    assert baseline.power_over_onsite_speed == pytest.approx(
        0.00028755957070134653,
        abs=3.0e-17,
    )
    assert doubled.newton_constant == 2 * baseline.newton_constant
    assert doubled.source == baseline.source
    assert doubled.power == pytest.approx(2.0 * baseline.power, rel=3.0e-16)
    assert doubled.field_cycle_energy_loss_fraction == pytest.approx(
        2.0 * baseline.field_cycle_energy_loss_fraction,
        rel=3.0e-16,
    )


def test_small_width_transform_has_the_analytic_asymptote() -> None:
    inverse_width = sp.Rational(1, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    exact_transform = leading_small_amplitude_breather_retarded_integral(
        frequency,
        direction_points=4001,
    )
    asymptotic = (
        512.0
        * np.pi
        * float(inverse_width) ** 3
        * float(frequency)
        / 3.0
    )
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
        sine_gordon_breather_retarded_stress_evidence(
            sp.Rational(4, 5),
            **arguments,
        )


def test_radiation_requires_exact_frequency_and_explicit_positive_g() -> None:
    with pytest.raises(ValueError, match="frequency"):
        sine_gordon_breather_retarded_stress_evidence(0.8)
    with pytest.raises(ValueError, match="newton_constant"):
        sine_gordon_breather_radiation_evidence(
            sp.Rational(4, 5),
            dimensional_sine_gordon_coefficients(1, 1, 1),
            -1,
            spatial_points=257,
            time_points=33,
            direction_points=129,
            harmonic_count=2,
        )
