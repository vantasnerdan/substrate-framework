from __future__ import annotations

import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.linearized_einstein import (
    HarmonicGaugeEinsteinLedger,
    LongitudinalLineRadiationLedger,
    WeakFieldMonopole,
    harmonic_gauge_einstein_ledger,
    longitudinal_line_radiation_ledger,
    weak_field_monopole,
)


def test_public_package_exports_linearized_einstein_api() -> None:
    assert framework.HarmonicGaugeEinsteinLedger is HarmonicGaugeEinsteinLedger
    assert framework.LongitudinalLineRadiationLedger is LongitudinalLineRadiationLedger
    assert framework.WeakFieldMonopole is WeakFieldMonopole
    assert framework.harmonic_gauge_einstein_ledger is harmonic_gauge_einstein_ledger
    assert (
        framework.longitudinal_line_radiation_ledger
        is longitudinal_line_radiation_ledger
    )
    assert framework.weak_field_monopole is weak_field_monopole


def test_harmonic_gauge_coefficients_share_one_newton_normalization() -> None:
    newton, speed = sp.symbols("G c", positive=True)
    ledger = harmonic_gauge_einstein_ledger(newton, speed)
    assert ledger.einstein_stress_coupling == 8 * sp.pi * newton / speed**4
    assert ledger.trace_reversed_wave_source_coefficient == (
        -16 * sp.pi * newton / speed**4
    )
    assert ledger.retarded_far_zone_coefficient == 4 * newton / speed**4
    assert ledger.isaacson_flux_coefficient == speed**3 / (32 * sp.pi * newton)


def test_static_monopole_follows_from_the_linearized_source_equation() -> None:
    newton, speed, mass, radius = sp.symbols("G c M r", positive=True)
    field = weak_field_monopole(newton, speed, mass, radius)
    assert field.newtonian_potential == -newton * mass / radius
    assert field.potential_over_speed_squared == -newton * mass / (radius * speed**2)
    assert field.trace_reversed_metric_00 == -4 * newton * mass / (radius * speed**2)
    assert field.metric_00 == 1 - 2 * newton * mass / (radius * speed**2)
    assert field.spatial_metric_factor == 1 + 2 * newton * mass / (radius * speed**2)


def test_full_retardation_line_source_has_the_exact_tt_angular_weight() -> None:
    newton, speed = sp.symbols("G c", positive=True)
    cosine = sp.symbols("zeta", real=True)
    ledger = longitudinal_line_radiation_ledger(newton, speed, cosine)
    assert ledger.projected_norm_squared == (1 - cosine**2) ** 2 / 2
    assert ledger.differential_solid_angle_power_coefficient == (
        newton * (1 - cosine**2) ** 2 / (4 * sp.pi * speed**5)
    )
    assert ledger.azimuth_power_coefficient == (
        newton * (1 - cosine**2) ** 2 / (2 * speed**5)
    )
    assert (
        longitudinal_line_radiation_ledger(newton, speed, 1).projected_norm_squared
        == 0
    )
    assert longitudinal_line_radiation_ledger(
        newton, speed, 0
    ).projected_norm_squared == sp.Rational(1, 2)


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: harmonic_gauge_einstein_ledger(0, 1), "newton_constant"),
        (lambda: weak_field_monopole(1, 1, 1, 0), "radius"),
        (lambda: longitudinal_line_radiation_ledger(1, 1, 2), "direction_cosine"),
    ],
)
def test_linearized_einstein_inputs_are_guarded(call, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()
