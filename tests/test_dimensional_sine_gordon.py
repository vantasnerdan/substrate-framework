from __future__ import annotations

import pytest
import sympy as sp
import substrate_framework as framework

from substrate_framework.dimensional_sine_gordon import (
    DimensionalBreatherObservables,
    DimensionalSineGordonCoefficients,
    DimensionalSineGordonLinearSpectrum,
    DimensionalSineGordonScales,
    DimensionalSineGordonStress,
    DimensionalSineGordonTimeHarmonicLedger,
    dimensional_breather_field,
    dimensional_breather_observables,
    dimensional_sine_gordon_coefficient_dimension_matrix,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_coefficients_from_speed_gap,
    dimensional_sine_gordon_hamiltonian_density,
    dimensional_sine_gordon_lagrangian_density,
    dimensional_sine_gordon_log_ratio_jacobian,
    dimensional_sine_gordon_linear_spectrum,
    dimensional_sine_gordon_linearized_residual,
    dimensional_sine_gordon_normalized_coordinates,
    dimensional_sine_gordon_physical_coordinates,
    dimensional_sine_gordon_residual,
    dimensional_sine_gordon_scales,
    dimensional_sine_gordon_stress,
    dimensional_sine_gordon_tail_coefficient,
    dimensional_sine_gordon_time_harmonic_ledger,
    evanescent_half_line_matching_matrix,
    linear_wave_energy_density,
    linear_wave_residual,
    linear_wave_traveling_field,
    rescale_dimensional_sine_gordon_coefficients,
)


def test_public_package_exports_dimensional_sine_gordon_api() -> None:
    assert framework.DimensionalSineGordonCoefficients is DimensionalSineGordonCoefficients
    assert framework.DimensionalSineGordonScales is DimensionalSineGordonScales
    assert framework.DimensionalSineGordonStress is DimensionalSineGordonStress
    assert framework.DimensionalSineGordonLinearSpectrum is DimensionalSineGordonLinearSpectrum
    assert (
        framework.DimensionalSineGordonTimeHarmonicLedger
        is DimensionalSineGordonTimeHarmonicLedger
    )
    assert framework.DimensionalBreatherObservables is DimensionalBreatherObservables
    assert framework.dimensional_breather_field is dimensional_breather_field
    assert framework.dimensional_breather_observables is dimensional_breather_observables
    assert framework.dimensional_sine_gordon_stress is dimensional_sine_gordon_stress


def test_coefficient_ratios_and_inverse_family_retain_the_common_scale() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    scales = dimensional_sine_gordon_scales(coefficients)
    assert scales.signal_speed == 2
    assert scales.gap_frequency == 3
    assert scales.length == sp.Rational(2, 3)
    assert scales.energy == 12
    assert scales.action == 4
    assert scales.length == scales.signal_speed / scales.gap_frequency
    assert scales.energy == scales.gap_frequency * scales.action

    inverse = dimensional_sine_gordon_coefficients_from_speed_gap(2, 2, 3)
    assert inverse == coefficients


def test_common_multiplier_preserves_dynamics_but_rescales_energy_and_action() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    scaled = rescale_dimensional_sine_gordon_coefficients(coefficients, 5)
    original_scales = dimensional_sine_gordon_scales(coefficients)
    scaled_scales = dimensional_sine_gordon_scales(scaled)
    assert scaled_scales.signal_speed == original_scales.signal_speed
    assert scaled_scales.gap_frequency == original_scales.gap_frequency
    assert scaled_scales.length == original_scales.length
    assert scaled_scales.energy == 5 * original_scales.energy
    assert scaled_scales.action == 5 * original_scales.action


def test_ratio_jacobian_exposes_one_common_coefficient_direction() -> None:
    jacobian = dimensional_sine_gordon_log_ratio_jacobian()
    assert jacobian.rank() == 2
    assert jacobian * sp.ones(3, 1) == sp.zeros(3, 1)
    assert jacobian.T.nullspace() == [sp.Matrix([-1, 1, 1])]


def test_coefficient_dimension_matrix_matches_energy_density_convention() -> None:
    matrix = dimensional_sine_gordon_coefficient_dimension_matrix()
    assert matrix == sp.Matrix([[1, 1, 1], [-1, 1, -1], [2, 0, 0]])
    assert matrix.det() == -4
    assert matrix.rank() == 3
    assert matrix.nullspace() == []


def test_coordinate_maps_are_exact_inverses() -> None:
    x, t, X, tau = sp.symbols("x t X tau", real=True)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    normalized = dimensional_sine_gordon_normalized_coordinates(x, t, coefficients)
    assert normalized == (sp.Rational(3, 2) * x, 3 * t)
    physical = dimensional_sine_gordon_physical_coordinates(*normalized, coefficients)
    assert physical == (x, t)
    inverse_physical = dimensional_sine_gordon_physical_coordinates(X, tau, coefficients)
    assert dimensional_sine_gordon_normalized_coordinates(
        *inverse_physical,
        coefficients,
    ) == (X, tau)


def test_declared_lagrangian_and_hamiltonian_densities_keep_all_coefficients() -> None:
    x, t = sp.symbols("x t", real=True)
    field = sp.Function("u")(x, t)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    kinetic = sp.diff(field, t) ** 2
    gradient = sp.diff(field, x) ** 2
    potential = 1 - sp.cos(field)
    assert dimensional_sine_gordon_lagrangian_density(
        field,
        x,
        t,
        coefficients,
    ) == kinetic - 4 * gradient - 18 * potential
    assert dimensional_sine_gordon_hamiltonian_density(
        field,
        x,
        t,
        coefficients,
    ) == kinetic + 4 * gradient + 18 * potential


def test_physical_stress_has_the_exact_off_shell_divergence() -> None:
    x, t = sp.symbols("x t", real=True)
    field = sp.Function("u")(x, t)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    field_t = sp.diff(field, t)
    field_x = sp.diff(field, x)
    potential = 1 - sp.cos(field)
    residual = 2 * sp.diff(field, t, 2) - 8 * sp.diff(field, x, 2) + 18 * sp.sin(field)
    stress = dimensional_sine_gordon_stress(field, x, t, coefficients)
    assert stress.energy_density == field_t**2 + 4 * field_x**2 + 18 * potential
    assert stress.energy_flux == -8 * field_t * field_x
    assert stress.longitudinal_pressure == field_t**2 + 4 * field_x**2 - 18 * potential
    assert stress.contravariant == sp.Matrix(
        [
            [stress.energy_density, -4 * field_t * field_x],
            [-4 * field_t * field_x, stress.longitudinal_pressure],
        ]
    )
    assert stress.field_equation_residual == residual
    assert stress.divergence == sp.Matrix([field_t * residual / 2, -field_x * residual])


def test_pulled_back_breather_solves_the_dimensional_equation() -> None:
    x, t = sp.symbols("x t", real=True)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    field = dimensional_breather_field(x, t, sp.Rational(1, 2), coefficients)
    assert sp.simplify(
        dimensional_sine_gordon_residual(field, x, t, coefficients)
    ) == 0


def test_linearized_residual_and_plane_wave_dispersion_are_exact() -> None:
    x, t = sp.symbols("x t", real=True)
    k, angular = sp.symbols("k Omega", real=True)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    field = sp.exp(sp.I * (k * x - angular * t))
    characteristic = sp.simplify(
        dimensional_sine_gordon_linearized_residual(
            field,
            x,
            t,
            coefficients,
        )
        / field
    )
    assert characteristic == -2 * angular**2 + 8 * k**2 + 18
    assert sp.simplify(characteristic.subs(angular**2, 9 + 4 * k**2)) == 0


def test_linear_spectrum_tracks_gap_group_and_phase_velocities() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    zero = dimensional_sine_gordon_linear_spectrum(0, coefficients)
    assert zero.angular_frequency == 3
    assert zero.group_velocity == 0
    assert zero.phase_velocity is None

    spectrum = dimensional_sine_gordon_linear_spectrum(4, coefficients)
    assert spectrum.angular_frequency_squared == 73
    assert spectrum.angular_frequency == sp.sqrt(73)
    assert spectrum.group_velocity == 16 / sp.sqrt(73)
    assert spectrum.phase_velocity == sp.sqrt(73) / 4
    assert sp.simplify(spectrum.group_velocity * spectrum.phase_velocity - 4) == 0


def test_time_harmonic_ledger_separates_half_line_tails_from_global_l2_modes() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    subgap = dimensional_sine_gordon_time_harmonic_ledger(1, coefficients)
    assert subgap.behavior == "evanescent"
    assert subgap.spatial_coefficient == 2
    assert subgap.spatial_rate == sp.sqrt(2)
    assert subgap.right_half_line_l2_dimension == 1
    assert subgap.left_half_line_l2_dimension == 1
    assert subgap.whole_line_l2_dimension == 0

    threshold = dimensional_sine_gordon_time_harmonic_ledger(3, coefficients)
    oscillatory = dimensional_sine_gordon_time_harmonic_ledger(5, coefficients)
    assert threshold.behavior == "threshold"
    assert threshold.spatial_rate == 0
    assert threshold.whole_line_l2_dimension == 0
    assert oscillatory.behavior == "oscillatory"
    assert oscillatory.spatial_rate == 2
    assert oscillatory.whole_line_l2_dimension == 0


def test_evanescent_half_line_branches_have_only_the_zero_smooth_global_match() -> None:
    kappa = sp.symbols("kappa", positive=True, real=True)
    matching = evanescent_half_line_matching_matrix(kappa)
    assert matching == sp.Matrix([[1, -1], [-kappa, -kappa]])
    assert matching.det() == -2 * kappa
    assert matching.rank() == 2
    assert matching.nullspace() == []


def test_tail_coefficient_matches_exact_breather_inverse_width() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    omega = sp.Rational(1, 2)
    physical_angular = 3 * omega
    coefficient = dimensional_sine_gordon_tail_coefficient(
        physical_angular,
        coefficients,
    )
    observables = dimensional_breather_observables(omega, coefficients)
    assert sp.simplify(coefficient - observables.inverse_width**2) == 0


def test_gapless_traveling_profile_is_exact_and_has_translated_energy_density() -> None:
    x, t, z = sp.symbols("x t z", real=True)
    profile = sp.Function("F")
    field = linear_wave_traveling_field(profile, x, t, 2)
    assert field == profile(x - 2 * t)
    assert sp.simplify(linear_wave_residual(field, x, t, 2)) == 0
    density = linear_wave_energy_density(field, x, t, 2, 8)
    expected = 8 * sp.Subs(sp.Derivative(profile(z), z), z, x - 2 * t) ** 2
    assert sp.simplify(density - expected) == 0


def test_dimensional_breather_observables_include_physical_scales() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    observables = dimensional_breather_observables(sp.Rational(1, 2), coefficients)
    assert observables.angular_frequency == sp.Rational(3, 2)
    assert observables.period == 4 * sp.pi / 3
    assert observables.inverse_width == 3 * sp.sqrt(3) / 4
    assert sp.simplify(observables.profile_length - 4 / (3 * sp.sqrt(3))) == 0
    assert observables.energy == 96 * sp.sqrt(3)
    assert observables.action == 64 * sp.pi / 3


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: dimensional_sine_gordon_coefficients(0, 1, 1), "inertia"),
        (lambda: dimensional_sine_gordon_coefficients(1, -1, 1), "gradient"),
        (lambda: dimensional_sine_gordon_coefficients(1, 1, 1.0), "exact"),
        (
            lambda: dimensional_sine_gordon_coefficients_from_speed_gap(1, 0, 1),
            "signal_speed",
        ),
        (
            lambda: dimensional_breather_observables(
                1,
                dimensional_sine_gordon_coefficients(1, 1, 1),
            ),
            "frequency",
        ),
        (
            lambda: rescale_dimensional_sine_gordon_coefficients(
                dimensional_sine_gordon_coefficients(1, 1, 1),
                0,
            ),
            "multiplier",
        ),
        (
            lambda: dimensional_sine_gordon_residual(
                sp.Function("u")(sp.Symbol("x"), sp.Symbol("t")),
                sp.Symbol("x"),
                sp.Symbol("t"),
                DimensionalSineGordonCoefficients(1, -1, 1),
            ),
            "gradient",
        ),
        (
            lambda: dimensional_sine_gordon_time_harmonic_ledger(
                sp.Symbol("Omega", nonnegative=True),
                dimensional_sine_gordon_coefficients(1, 4, 9),
            ),
            "decidable relation",
        ),
        (lambda: evanescent_half_line_matching_matrix(0), "rate"),
        (
            lambda: linear_wave_traveling_field(
                sp.Function("F"),
                sp.Symbol("x"),
                sp.Symbol("t"),
                1,
                direction=0,
            ),
            "direction",
        ),
    ],
)
def test_exact_positive_domains_are_enforced(call, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()
