from __future__ import annotations

import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.linearized_einstein import (
    HarmonicGaugeEinsteinLedger,
    RetardedLineSourceRadiationLedger,
    WeakFieldMonopole,
    harmonic_gauge_einstein_ledger,
    retarded_line_source_radiation_ledger,
    weak_field_monopole,
)


def test_public_package_exports_linearized_einstein_api() -> None:
    assert framework.HarmonicGaugeEinsteinLedger is HarmonicGaugeEinsteinLedger
    assert (
        framework.RetardedLineSourceRadiationLedger
        is RetardedLineSourceRadiationLedger
    )
    assert framework.WeakFieldMonopole is WeakFieldMonopole
    assert framework.harmonic_gauge_einstein_ledger is harmonic_gauge_einstein_ledger
    assert (
        framework.retarded_line_source_radiation_ledger
        is retarded_line_source_radiation_ledger
    )
    assert framework.weak_field_monopole is weak_field_monopole


def test_harmonic_gauge_ledger_freezes_mostly_plus_signs() -> None:
    newton, speed = sp.symbols("G c", positive=True)
    ledger = harmonic_gauge_einstein_ledger(newton, speed)
    assert ledger.minkowski_metric_covariant == sp.diag(-1, 1, 1, 1)
    assert ledger.minkowski_metric_covariant != sp.diag(1, -1, -1, -1)
    assert ledger.d_alembertian_time_coefficient == -1 / speed**2
    assert ledger.d_alembertian_spatial_coefficient == 1
    assert ledger.einstein_stress_coupling == 8 * sp.pi * newton / speed**4
    assert ledger.trace_reversed_wave_source_coefficient == (
        -16 * sp.pi * newton / speed**4
    )
    assert ledger.retarded_far_zone_coefficient == 4 * newton / speed**4
    assert sp.simplify(
        ledger.retarded_far_zone_coefficient
        + ledger.trace_reversed_wave_source_coefficient / (4 * sp.pi)
    ) == 0
    assert ledger.isaacson_flux_coefficient == speed**3 / (32 * sp.pi * newton)


def test_static_monopole_closes_the_mostly_plus_poisson_equation() -> None:
    newton, speed, mass, radius = sp.symbols("G c M r", positive=True)
    field = weak_field_monopole(newton, speed, mass, radius)
    assert field.newtonian_potential == -newton * mass / radius
    assert field.potential_over_speed_squared == -newton * mass / (
        radius * speed**2
    )
    assert field.trace_reversed_poisson_source_coefficient == (
        -16 * sp.pi * newton * mass / speed**2
    )
    assert sp.simplify(
        -4 * sp.pi * radius * field.trace_reversed_metric_00
        - field.trace_reversed_poisson_source_coefficient
    ) == 0
    assert field.trace_reversed_metric_00 == 4 * newton * mass / (
        radius * speed**2
    )
    assert field.metric_perturbation_00 == 2 * newton * mass / (
        radius * speed**2
    )
    assert field.metric_00 == -1 + 2 * newton * mass / (radius * speed**2)
    assert field.spatial_metric_factor == 1 + 2 * newton * mass / (
        radius * speed**2
    )
    assert sp.limit(field.metric_00, radius, sp.oo) == -1
    assert sp.limit(field.spatial_metric_factor, radius, sp.oo) == 1


def test_mostly_minus_static_sign_mutations_break_the_source_closure() -> None:
    newton, speed, mass, radius = sp.symbols("G c M r", positive=True)
    field = weak_field_monopole(newton, speed, mass, radius)
    reduced = field.potential_over_speed_squared
    wrong_trace_reversed_metric_00 = 4 * reduced
    wrong_metric_00 = 1 + 2 * reduced
    assert sp.simplify(
        -4 * sp.pi * radius * wrong_trace_reversed_metric_00
        - field.trace_reversed_poisson_source_coefficient
    ) != 0
    assert sp.simplify(field.metric_00 - wrong_metric_00) != 0


def test_retarded_line_source_has_the_exact_tt_angular_weight() -> None:
    newton, speed = sp.symbols("G c", positive=True)
    cosine = sp.symbols("zeta", real=True)
    ledger = retarded_line_source_radiation_ledger(newton, speed, cosine)
    assert ledger.projected_norm_squared == (1 - cosine**2) ** 2 / 2
    assert ledger.differential_solid_angle_power_coefficient == (
        newton * (1 - cosine**2) ** 2 / (4 * sp.pi * speed**5)
    )
    assert ledger.azimuth_power_coefficient == (
        newton * (1 - cosine**2) ** 2 / (2 * speed**5)
    )
    assert (
        retarded_line_source_radiation_ledger(
            newton, speed, 1
        ).projected_norm_squared
        == 0
    )
    assert retarded_line_source_radiation_ledger(
        newton, speed, 0
    ).projected_norm_squared == sp.Rational(1, 2)


def test_retarded_source_sign_mutation_breaks_the_green_function_relation() -> None:
    newton, speed = sp.symbols("G c", positive=True)
    ledger = harmonic_gauge_einstein_ledger(newton, speed)
    wrong_retarded_coefficient = -ledger.retarded_far_zone_coefficient
    assert sp.simplify(
        wrong_retarded_coefficient
        + ledger.trace_reversed_wave_source_coefficient / (4 * sp.pi)
    ) != 0


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: harmonic_gauge_einstein_ledger(0, 1), "newton_constant"),
        (lambda: harmonic_gauge_einstein_ledger(1, 0), "signal_speed"),
        (lambda: weak_field_monopole(1, 1, 0, 1), "mass"),
        (lambda: weak_field_monopole(1, 1, 1, 0), "radius"),
        (
            lambda: retarded_line_source_radiation_ledger(1, 1, 2),
            r"direction_cosine must lie in \[-1, 1\]; got 2",
        ),
    ],
)
def test_linearized_einstein_inputs_are_guarded(call, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()
