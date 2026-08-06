#!/usr/bin/env python3
"""Primary verifier for the proposed P226 sine-Gordon gravity chain."""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
from scipy.integrate import quad
import sympy as sp
import yaml

from substrate_framework.dimensional_sine_gordon import (
    dimensional_sine_gordon_coefficients,
)
from substrate_framework.local_horizon_gravity import local_horizon_einstein_ledger
from substrate_framework.linearized_einstein import (
    harmonic_gauge_einstein_ledger,
    longitudinal_line_radiation_ledger,
)
from substrate_framework.sine_gordon_breather_gravity import (
    breather_weak_field_prediction,
)
from substrate_framework.sine_gordon_breather_radiation import (
    breather_radiation_evidence,
    leading_small_amplitude_breather_power_ratio,
)
from substrate_framework.sine_gordon_fiber_source import (
    canonical_sine_gordon_fiber_stress,
)
from substrate_framework.sine_gordon_weave import (
    collective_weave_metric,
    sine_gordon_weave_gravity,
)
from substrate_framework.source_audit import audit_numpy_trapezoid_compatibility
from substrate_framework.verification import CheckLedger


ROOT = Path(__file__).resolve().parents[2]
PROPOSAL = ROOT / "proposals/P226-sg-emergent-gravity"


def _independent_leading_power(frequency: float, inverse_width: float) -> float:
    scale = np.pi * frequency / inverse_width

    def integrand(direction: float) -> float:
        if direction == 0.0:
            ratio = 1.0 / scale**2
        else:
            ratio = (direction / np.sinh(scale * direction)) ** 2
        return ratio * (1.0 - direction**2) ** 2

    integral, _ = quad(integrand, 0.0, 1.0, epsabs=1.0e-14, epsrel=1.0e-13)
    return 128.0 * np.pi**2 * frequency**4 * integral / np.log(2.0)


def main() -> int:
    checks = CheckLedger("P226-C-EGR-001-004")
    manifest = yaml.safe_load((PROPOSAL / "proposal.yaml").read_text())
    checks.check(
        "proposal reserves exactly the four emergent-gravity claims",
        manifest["claims_proposed"]
        == ["C-EGR-001", "C-EGR-002", "C-EGR-003", "C-EGR-004"],
    )
    comparison = yaml.safe_load(
        (PROPOSAL / "evidence/candidate-comparison.yaml").read_text()
    )
    checks.check(
        "candidate selection preceded every empirical comparator",
        comparison["decisions"]["A"]["decision"] == "selected"
        and not comparison["comparator_selected_candidate"]
        and comparison["empirical_comparators_opened"] == [],
    )

    lam, tension, onsite = sp.symbols("lambda T mu", positive=True)
    coefficients = dimensional_sine_gordon_coefficients(lam, tension, onsite)
    weave = sine_gordon_weave_gravity(coefficients)
    expected_newton = tension**2 / (
        4 * lam**2 * onsite * sp.log(2)
    )
    checks.check(
        "one microscopic coefficient triple fixes cell speed action and total G",
        weave.geometry.cell_length == sp.sqrt(tension / onsite)
        and weave.scales.signal_speed == sp.sqrt(tension / lam)
        and weave.scales.action == sp.sqrt(lam * tension)
        and sp.simplify(
            weave.horizon_equilibrium.newton_constant - expected_newton
        )
        == 0,
    )
    checks.check(
        "the Einstein source coefficient closes without an additive baseline",
        sp.simplify(
            weave.horizon_equilibrium.einstein_stress_coupling
            - 2 * sp.pi / (onsite * sp.log(2))
        )
        == 0
        and tuple(inspect.signature(sine_gordon_weave_gravity).parameters)
        == ("coefficients",),
    )
    checks.mutation_sensitive(
        "topological sector degeneracy fixes the entropy coefficient",
        lambda sector_count: sp.simplify(
            local_horizon_einstein_ledger(
                sp.log(sector_count) / weave.scales.length**2,
                weave.scales.action,
                weave.scales.signal_speed,
            ).newton_constant
            - weave.horizon_equilibrium.newton_constant
        )
        == 0,
        2,
        [3, 4],
    )
    checks.mutation_sensitive(
        "the unique sine-Gordon action is the Unruh action",
        lambda action: sp.simplify(
            local_horizon_einstein_ledger(
                weave.entropy_area_density,
                action,
                weave.scales.signal_speed,
            ).newton_constant
            - weave.horizon_equilibrium.newton_constant
        )
        == 0,
        weave.scales.action,
        [2 * weave.scales.action, weave.scales.action / 2],
    )

    time_scale, x_scale, y_scale, z_scale = sp.symbols(
        "a_t a_x a_y a_z", positive=True
    )
    collective_metric = collective_weave_metric(
        sp.diag(time_scale, x_scale, y_scale, z_scale)
    )
    checks.check(
        "clock and tetrahedral link moments reconstruct a 3+1 metric and measure",
        collective_metric.covariant_metric
        == sp.diag(time_scale**2, -x_scale**2, -y_scale**2, -z_scale**2)
        and collective_metric.link_second_moment
        == sp.diag(0, x_scale**2 / 3, y_scale**2 / 3, z_scale**2 / 3)
        and collective_metric.volume_density
        == time_scale * x_scale * y_scale * z_scale,
    )
    checks.check(
        "local equilibrium yields nonlinear Einstein dynamics and two tensor modes",
        weave.horizon_equilibrium.cosmological_constant == 0
        and weave.horizon_equilibrium.radiative_tensor_polarizations == 2
        and sp.simplify(
            weave.horizon_equilibrium.einstein_stress_coupling
            - 8
            * sp.pi
            * weave.horizon_equilibrium.newton_constant
            / weave.scales.signal_speed**4
        )
        == 0,
    )

    x, t, y, z = sp.symbols("x t y z", real=True)
    field = sp.Function("u")(x, t)
    source = canonical_sine_gordon_fiber_stress(
        field,
        x,
        t,
        y,
        z,
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    checks.check(
        "the transverse-delta fiber source is isolated and conserved on shell",
        source.stress_energy[2:, :] == sp.zeros(2, 4)
        and source.stress_energy[:, 2:] == sp.zeros(4, 2)
        and source.divergence.subs(
            source.sine_gordon_stress.field_equation_residual,
            0,
        )
        == sp.zeros(4, 1),
    )

    benchmark_coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    inverse_width = sp.Rational(3, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    profile_length = 1 / inverse_width
    static = breather_weak_field_prediction(
        frequency,
        benchmark_coefficients,
        profile_length,
    )
    checks.check(
        "the exact breather energy yields a controlled static weak field",
        static.observables.energy == 16 * inverse_width
        and static.field.potential_over_speed_squared
        == -4 * inverse_width**2 / sp.log(2)
        and static.profile_compactness == 8 * inverse_width**2 / sp.log(2)
        and float(static.profile_compactness) < 0.011,
    )
    linearized = harmonic_gauge_einstein_ledger(
        static.weave_gravity.horizon_equilibrium.newton_constant,
        static.weave_gravity.scales.signal_speed,
    )
    checks.check(
        "static and retarded limits use the same derived Einstein normalization",
        linearized.einstein_stress_coupling
        == static.weave_gravity.horizon_equilibrium.einstein_stress_coupling
        and longitudinal_line_radiation_ledger(
            linearized.newton_constant,
            linearized.signal_speed,
            0,
        ).retarded_far_zone_coefficient
        == linearized.retarded_far_zone_coefficient,
    )

    coarse = breather_radiation_evidence(
        frequency,
        benchmark_coefficients,
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    refined = breather_radiation_evidence(
        frequency,
        benchmark_coefficients,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    checks.check(
        "exact-stress full-retardation power is positive and refinement stable",
        refined.power_over_onsite_speed > 0.0
        and abs(
            refined.power_over_onsite_speed - coarse.power_over_onsite_speed
        )
        < 3.0e-10
        and abs(refined.power_over_onsite_speed - 0.002606646892787329)
        < 2.0e-15,
    )
    independent_leading = _independent_leading_power(
        float(frequency),
        float(inverse_width),
    )
    package_leading = leading_small_amplitude_breather_power_ratio(
        frequency,
        direction_points=4001,
    )
    checks.check(
        "independent adaptive integration reproduces the analytic leading power",
        abs(independent_leading - package_leading) < 2.0e-15
        and refined.relative_leading_difference < 7.3e-4,
    )
    checks.check(
        "full retardation is load bearing and one-way backreaction stays controlled",
        refined.source_size_frequency_ratio > 30.0
        and refined.field_cycle_energy_loss_fraction < 0.035
        and refined.harmonic_pressure_maxima[1] < 3.3e-6,
    )

    radiation_source = (
        ROOT / "src/substrate_framework/sine_gordon_breather_radiation.py"
    ).read_text()
    compatibility = audit_numpy_trapezoid_compatibility(
        radiation_source,
        filename="sine_gordon_breather_radiation.py",
    )
    checks.check(
        "radiation code reuses the numerical integration owner",
        compatibility.legacy_references == 0
        and compatibility.current_references == 0
        and compatibility.eager_legacy_default_fallbacks == 0,
    )
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
