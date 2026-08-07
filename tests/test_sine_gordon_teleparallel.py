from __future__ import annotations

import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import (
    dimensional_sine_gordon_coefficients,
    rescale_dimensional_sine_gordon_coefficients,
)
from substrate_framework.linearized_einstein import harmonic_gauge_einstein_ledger
from substrate_framework.sine_gordon_teleparallel import (
    CompactTeleparallelActionLedger,
    CompactTorsionChannelLedger,
    SineGordonLinkCoframe,
    SineGordonTeleparallelGravity,
    TeleparallelCoframeLedger,
    collective_sine_gordon_link_coframe,
    compact_teleparallel_action_ledger,
    compact_torsion_channel_ledger,
    sine_gordon_teleparallel_gravity,
    teleparallel_constitutive_matrix_mostly_plus,
    teleparallel_constitutive_spectral_basis,
    teleparallel_coframe_ledger,
)


def test_public_package_exports_teleparallel_construction() -> None:
    assert framework.CompactTeleparallelActionLedger is CompactTeleparallelActionLedger
    assert framework.CompactTorsionChannelLedger is CompactTorsionChannelLedger
    assert framework.SineGordonLinkCoframe is SineGordonLinkCoframe
    assert framework.SineGordonTeleparallelGravity is SineGordonTeleparallelGravity
    assert framework.TeleparallelCoframeLedger is TeleparallelCoframeLedger
    assert (
        framework.collective_sine_gordon_link_coframe
        is collective_sine_gordon_link_coframe
    )
    assert (
        framework.compact_teleparallel_action_ledger
        is compact_teleparallel_action_ledger
    )
    assert framework.compact_torsion_channel_ledger is compact_torsion_channel_ledger
    assert framework.sine_gordon_teleparallel_gravity is sine_gordon_teleparallel_gravity
    assert (
        framework.teleparallel_constitutive_matrix_mostly_plus
        is teleparallel_constitutive_matrix_mostly_plus
    )
    assert (
        framework.teleparallel_constitutive_spectral_basis
        is teleparallel_constitutive_spectral_basis
    )
    assert framework.teleparallel_coframe_ledger is teleparallel_coframe_ledger


def test_microscopic_link_covectors_construct_the_collective_metric_and_measure() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    length = sp.Rational(2, 3)
    time_scale, x_scale, y_scale, z_scale = sp.symbols(
        "a_t a_x a_y a_z",
        positive=True,
    )
    links = length * sp.diag(time_scale, x_scale, y_scale, z_scale)
    geometry = collective_sine_gordon_link_coframe(links, coefficients)
    assert geometry.microscopic_link_covectors == links
    assert geometry.cell_length == length
    assert geometry.collective_coframe == sp.diag(
        time_scale,
        x_scale,
        y_scale,
        z_scale,
    )
    assert geometry.metric_covariant == sp.diag(
        -time_scale**2,
        x_scale**2,
        y_scale**2,
        z_scale**2,
    )
    assert geometry.volume_density == time_scale * x_scale * y_scale * z_scale


def test_cosine_cell_has_an_exact_chord_torsion_quadratic_form() -> None:
    phase = sp.symbols("q", real=True)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    channel = compact_torsion_channel_ledger(
        phase,
        coefficients,
        constitutive_weight=sp.Rational(-1, 2),
    )
    assert channel.cell_length == sp.Rational(2, 3)
    assert channel.transverse_line_density == sp.Rational(9, 4)
    assert channel.chord_torsion == 3 * sp.sin(phase / 2)
    assert channel.geometric_small_phase_torsion == 3 * phase / 2
    assert channel.microscopic_lagrangian_density == (
        sp.Rational(81, 4) * (1 - sp.cos(phase))
    )
    assert channel.quadratic_collective_density == (
        sp.Rational(81, 2) * sp.sin(phase / 2) ** 2
    )
    assert channel.cosine_chord_identity_residual == 0
    relative_chord = sp.simplify(
        channel.chord_torsion / channel.geometric_small_phase_torsion
    )
    assert sp.series(relative_chord, phase, 0, 5) == (
        1 - phase**2 / 24 + phase**4 / 1920 + sp.Order(phase**5)
    )


def test_all_compact_channels_have_an_explicit_periodic_map_to_tegr() -> None:
    phase = sp.symbols("q", real=True)
    phases = [0] * 24
    phases[0] = phase
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    action = compact_teleparallel_action_ledger(phases, coefficients)
    constitutive = teleparallel_constitutive_matrix_mostly_plus()
    basis, spectral = teleparallel_constitutive_spectral_basis()
    assert action.constitutive_matrix == constitutive
    assert action.spectral_basis == basis
    assert action.spectral_weights == spectral
    assert constitutive.shape == (24, 24)
    assert constitutive == constitutive.T
    assert constitutive.rank() == 24
    eigenvalues = list(spectral.diagonal())
    assert eigenvalues.count(-2) == 3
    assert eigenvalues.count(-1) == 8
    assert eigenvalues.count(sp.Rational(-1, 2)) == 1
    assert eigenvalues.count(sp.Rational(1, 2)) == 3
    assert eigenvalues.count(1) == 8
    assert eigenvalues.count(2) == 1
    assert basis.T * basis == sp.eye(24)
    assert sp.simplify(basis.T * constitutive * basis - spectral) == sp.zeros(24)
    assert action.microscopic_lagrangian_density == 81 * (1 - sp.cos(phase))
    assert action.collective_teleparallel_density == 162 * sp.sin(phase / 2) ** 2
    assert action.compact_collective_identity_residual == 0

    shifted = compact_teleparallel_action_ledger(
        [phase + 2 * sp.pi, *([0] * 23)],
        coefficients,
    )
    assert sp.trigsimp(
        shifted.microscopic_lagrangian_density
        - action.microscopic_lagrangian_density
    ) == 0
    assert sp.trigsimp(
        shifted.collective_teleparallel_density
        - action.collective_teleparallel_density
    ) == 0


def test_sg_cell_scale_closes_the_remaining_teleparallel_multiplier() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    gravity = sine_gordon_teleparallel_gravity(coefficients)
    assert gravity.cell_length == sp.Rational(2, 3)
    assert gravity.torsion_invariant_weights == (
        sp.Rational(1, 4),
        sp.Rational(1, 2),
        -1,
    )
    assert gravity.quadratic_coframe_parameter_count == 3
    assert gravity.local_lorentz_ratio_dimension == 1
    assert gravity.teleparallel_action_coefficient == -9
    assert gravity.einstein_hilbert_coefficient == 9
    assert gravity.einstein_stress_coupling == sp.Rational(1, 18)
    assert gravity.newton_constant == 1 / (9 * sp.pi)
    assert gravity.radiative_tensor_polarizations == 2
    assert gravity.independent_gravitational_parameters == 0
    response = harmonic_gauge_einstein_ledger(
        gravity.newton_constant,
        gravity.scales.signal_speed,
    )
    assert response.einstein_stress_coupling == gravity.einstein_stress_coupling


def test_common_microscopic_multiplier_rescales_gravity_without_a_new_parameter() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    baseline = sine_gordon_teleparallel_gravity(coefficients)
    scaled = sine_gordon_teleparallel_gravity(
        rescale_dimensional_sine_gordon_coefficients(coefficients, 5)
    )
    assert scaled.scales.signal_speed == baseline.scales.signal_speed
    assert scaled.cell_length == baseline.cell_length
    assert scaled.einstein_hilbert_coefficient == 5 * baseline.einstein_hilbert_coefficient
    assert scaled.einstein_stress_coupling == baseline.einstein_stress_coupling / 5
    assert scaled.newton_constant == baseline.newton_constant / 5


def test_nonlinear_flrw_coframe_closes_the_mostly_plus_tegr_identity() -> None:
    x0, x, y, z = sp.symbols("x0 x y z", real=True)
    hubble = sp.symbols("H", real=True)
    scale = sp.exp(hubble * x0)
    ledger = teleparallel_coframe_ledger(
        sp.diag(1, scale, scale, scale),
        (x0, x, y, z),
    )
    assert ledger.metric_covariant == sp.diag(-1, scale**2, scale**2, scale**2)
    assert ledger.volume_density == scale**3
    assert ledger.torsion_invariant_one == -6 * hubble**2
    assert ledger.torsion_invariant_two == -3 * hubble**2
    assert ledger.torsion_vector_norm_squared == -9 * hubble**2
    assert ledger.torsion_scalar == 6 * hubble**2
    assert ledger.constitutive_quadratic_residual == 0
    assert ledger.boundary_divergence == 18 * hubble**2
    assert ledger.levi_civita_ricci_scalar == 12 * hubble**2
    assert ledger.einstein_teleparallel_identity_residual == 0


def test_local_frame_cancellations_and_invariant_weight_mutations_are_visible() -> None:
    x0, x, y, z = sp.symbols("x0 x y z", real=True)
    slope = sp.symbols("b", real=True)
    rapidity = slope * x0
    boost = sp.Matrix(
        [
            [sp.cosh(rapidity), sp.sinh(rapidity), 0, 0],
            [sp.sinh(rapidity), sp.cosh(rapidity), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    ledger = teleparallel_coframe_ledger(boost, (x0, x, y, z))
    assert ledger.metric_covariant == sp.diag(-1, 1, 1, 1)
    assert ledger.torsion_invariant_one == 2 * slope**2
    assert ledger.torsion_invariant_two == slope**2
    assert ledger.torsion_vector_norm_squared == slope**2
    assert ledger.torsion_scalar == 0
    assert ledger.einstein_teleparallel_identity_residual == 0

    mutated = teleparallel_coframe_ledger(
        boost,
        (x0, x, y, z),
        torsion_invariant_weights=(sp.Rational(1, 3), sp.Rational(1, 2), -1),
    )
    assert mutated.einstein_teleparallel_identity_residual == slope**2 / 6


@pytest.mark.parametrize(
    ("operation", "message"),
    [
        (
            lambda: compact_torsion_channel_ledger(
                0.2,
                dimensional_sine_gordon_coefficients(1, 1, 1),
            ),
            "phase",
        ),
        (
            lambda: teleparallel_coframe_ledger(sp.eye(3), sp.symbols("x:4")),
            "4 by 4",
        ),
        (
            lambda: teleparallel_coframe_ledger(sp.zeros(4), sp.symbols("x:4")),
            "invertible",
        ),
        (
            lambda: teleparallel_coframe_ledger(
                sp.eye(4),
                sp.symbols("x:4"),
                torsion_invariant_weights=(1, 2),
            ),
            "exactly three",
        ),
    ],
)
def test_teleparallel_inputs_fail_with_context(operation, message: str) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        operation()
