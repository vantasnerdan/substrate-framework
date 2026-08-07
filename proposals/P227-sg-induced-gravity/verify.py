#!/usr/bin/env python3
"""Primary verifier for proposed P227 claims C-EGR-005 through C-EGR-008."""

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
from substrate_framework.linearized_einstein import harmonic_gauge_einstein_ledger
from substrate_framework.sine_gordon_breather_gravity import (
    teleparallel_breather_weak_field_prediction,
)
from substrate_framework.sine_gordon_breather_radiation import (
    leading_small_amplitude_teleparallel_breather_power_ratio,
    teleparallel_breather_radiation_evidence,
)
from substrate_framework.sine_gordon_fiber_source import (
    canonical_sine_gordon_fiber_stress,
)
from substrate_framework.sine_gordon_teleparallel import (
    collective_sine_gordon_link_coframe,
    compact_teleparallel_action_ledger,
    compact_torsion_channel_ledger,
    sine_gordon_teleparallel_gravity,
    teleparallel_coframe_ledger,
    teleparallel_constitutive_matrix_mostly_plus,
    teleparallel_constitutive_spectral_basis,
)
from substrate_framework.source_audit import audit_numpy_trapezoid_compatibility
from substrate_framework.verification import CheckLedger


ROOT = Path(__file__).resolve().parents[2]
PROPOSAL = ROOT / "proposals/P227-sg-induced-gravity"


def _independent_leading_power(frequency: float, inverse_width: float) -> float:
    scale = np.pi * frequency / inverse_width

    def integrand(direction: float) -> float:
        if direction == 0.0:
            ratio = 1.0 / scale**2
        else:
            ratio = (direction / np.sinh(scale * direction)) ** 2
        return ratio * (1.0 - direction**2) ** 2

    half_integral, _ = quad(
        integrand,
        0.0,
        1.0,
        epsabs=1.0e-14,
        epsrel=1.0e-13,
    )
    return 64.0 * np.pi * frequency**4 * half_integral


def main() -> int:
    checks = CheckLedger("P227-C-EGR-005-008")
    manifest = yaml.safe_load((PROPOSAL / "proposal.yaml").read_text())
    checks.check(
        "proposal reserves exactly the four new emergent-gravity claims",
        manifest["claims_proposed"]
        == ["C-EGR-005", "C-EGR-006", "C-EGR-007", "C-EGR-008"],
    )
    comparison = yaml.safe_load(
        (PROPOSAL / "evidence/candidate-comparison.yaml").read_text()
    )
    checks.check(
        "all candidates were structurally decided without empirical selection",
        comparison["decisions"]["C"]["decision"] == "selected_for_implementation"
        and all(
            comparison["decisions"][candidate]["decision"].startswith("rejected")
            for candidate in ("A", "B", "D", "E")
        )
        and comparison["comparators_opened"] == [],
    )

    lam, tension, onsite = sp.symbols("lambda T mu", positive=True)
    coefficients = dimensional_sine_gordon_coefficients(lam, tension, onsite)
    gravity = sine_gordon_teleparallel_gravity(coefficients)
    checks.check(
        "the SG length density and onsite scale fix total G with no new parameter",
        gravity.cell_length == sp.sqrt(tension / onsite)
        and gravity.transverse_line_density == onsite / tension
        and gravity.einstein_hilbert_coefficient == onsite / 2
        and gravity.einstein_stress_coupling == 1 / onsite
        and sp.simplify(
            gravity.newton_constant
            - tension**2 / (8 * sp.pi * lam**2 * onsite)
        )
        == 0
        and gravity.independent_gravitational_parameters == 0
        and tuple(inspect.signature(sine_gordon_teleparallel_gravity).parameters)
        == ("coefficients",),
    )
    response = harmonic_gauge_einstein_ledger(
        gravity.newton_constant,
        gravity.scales.signal_speed,
    )
    checks.check(
        "the derived nonlinear and linearized Einstein coefficients agree",
        response.einstein_stress_coupling == gravity.einstein_stress_coupling
        and response.newton_constant == gravity.newton_constant,
    )
    time_scale, x_scale, y_scale, z_scale = sp.symbols(
        "a_t a_x a_y a_z",
        positive=True,
    )
    link_geometry = collective_sine_gordon_link_coframe(
        gravity.cell_length
        * sp.diag(time_scale, x_scale, y_scale, z_scale),
        coefficients,
    )
    checks.check(
        "microscopic cell links explicitly construct the coframe metric and measure",
        link_geometry.collective_coframe
        == sp.diag(time_scale, x_scale, y_scale, z_scale)
        and link_geometry.metric_covariant
        == sp.diag(-time_scale**2, x_scale**2, y_scale**2, z_scale**2)
        and link_geometry.volume_density
        == time_scale * x_scale * y_scale * z_scale,
    )

    constitutive = teleparallel_constitutive_matrix_mostly_plus()
    basis, spectral = teleparallel_constitutive_spectral_basis()
    checks.check(
        "the Lorentz-selected 24-channel constitutive map is explicit and invertible",
        constitutive.shape == (24, 24)
        and constitutive == constitutive.T
        and constitutive.rank() == 24
        and basis.T * basis == sp.eye(24)
        and sp.simplify(basis.T * constitutive * basis - spectral) == sp.zeros(24)
        and sorted(spectral.diagonal())
        == [
            -2,
            -2,
            -2,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            sp.Rational(-1, 2),
            sp.Rational(1, 2),
            sp.Rational(1, 2),
            sp.Rational(1, 2),
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
        ],
    )
    phase = sp.symbols("q", real=True)
    compact = compact_teleparallel_action_ledger(
        [phase, *([0] * 23)],
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    shifted = compact_teleparallel_action_ledger(
        [phase + 2 * sp.pi, *([0] * 23)],
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    checks.check(
        "compact phases map exactly and periodically to the TEGR quadratic density",
        compact.compact_collective_identity_residual == 0
        and compact.microscopic_lagrangian_density == 81 * (1 - sp.cos(phase))
        and sp.trigsimp(
            shifted.microscopic_lagrangian_density
            - compact.microscopic_lagrangian_density
        )
        == 0
        and sp.trigsimp(
            shifted.collective_teleparallel_density
            - compact.collective_teleparallel_density
        )
        == 0,
    )
    channel = compact_torsion_channel_ledger(
        phase,
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    checks.mutation_sensitive(
        "the sine-Gordon cell length is load bearing in the chord map",
        lambda length: sp.simplify(
            channel.chord_torsion - 2 * sp.sin(phase / 2) / length
        )
        == 0,
        sp.Rational(2, 3),
        [sp.Rational(1, 3), sp.Rational(4, 3)],
    )

    x0, x, y, z = sp.symbols("x0 x y z", real=True)
    hubble = sp.symbols("H", real=True)
    scale = sp.exp(hubble * x0)
    flrw = teleparallel_coframe_ledger(
        sp.diag(1, scale, scale, scale),
        (x0, x, y, z),
    )
    checks.check(
        "an independently constructed nonlinear coframe closes R plus T minus B",
        flrw.torsion_scalar == 6 * hubble**2
        and flrw.boundary_divergence == 18 * hubble**2
        and flrw.levi_civita_ricci_scalar == 12 * hubble**2
        and flrw.constitutive_quadratic_residual == 0
        and flrw.einstein_teleparallel_identity_residual == 0,
    )
    mutated = teleparallel_coframe_ledger(
        sp.diag(1, scale, scale, scale),
        (x0, x, y, z),
        torsion_invariant_weights=(sp.Rational(1, 3), sp.Rational(1, 2), -1),
    )
    checks.check(
        "a non-TEGR invariant ratio fails the nonlinear identity",
        mutated.einstein_teleparallel_identity_residual == -hubble**2 / 2,
    )

    source_field = sp.Function("u")(x, x0)
    source = canonical_sine_gordon_fiber_stress(
        source_field,
        x,
        x0,
        y,
        z,
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    checks.check(
        "the embedded microscopic source remains conserved on shell",
        source.divergence.subs(
            source.sine_gordon_stress.field_equation_residual,
            0,
        )
        == sp.zeros(4, 1),
    )

    benchmark_coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    inverse_width = sp.Rational(3, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    profile_length = 1 / inverse_width
    static = teleparallel_breather_weak_field_prediction(
        frequency,
        benchmark_coefficients,
        profile_length,
    )
    checks.check(
        "the same breather energy and derived G give a controlled static field",
        static.observables.energy == 16 * inverse_width
        and static.source_mass == 16 * inverse_width
        and static.field.potential_over_speed_squared
        == -2 * inverse_width**2 / sp.pi
        and static.profile_compactness == 4 * inverse_width**2 / sp.pi
        and float(static.profile_compactness) < 0.0012,
    )

    coarse = teleparallel_breather_radiation_evidence(
        frequency,
        benchmark_coefficients,
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    refined = teleparallel_breather_radiation_evidence(
        frequency,
        benchmark_coefficients,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    checks.check(
        "the exact-stress full-retardation power is positive and refinement stable",
        refined.power_over_onsite_speed > 0.0
        and abs(refined.power_over_onsite_speed - coarse.power_over_onsite_speed)
        < 4.0e-11
        and abs(refined.power_over_onsite_speed - 0.0002875595507880884)
        < 3.0e-17,
    )
    independent_leading = _independent_leading_power(
        float(frequency),
        float(inverse_width),
    )
    package_leading = leading_small_amplitude_teleparallel_breather_power_ratio(
        frequency,
        direction_points=4001,
    )
    checks.check(
        "adaptive quadrature independently reproduces the analytic leading transform",
        abs(independent_leading - package_leading) < 3.0e-16
        and refined.relative_leading_difference < 7.3e-4,
    )
    checks.check(
        "full retardation and one-way backreaction gates remain explicit",
        refined.source_size_frequency_ratio > 30.0
        and refined.field_cycle_energy_loss_fraction < 0.004
        and refined.harmonic_pressure_maxima[1] < 3.3e-6,
    )
    checks.check(
        "static and radiative consumers retain identical gravity provenance",
        refined.gravity == static.gravity
        and refined.source_energy == static.observables.energy,
    )

    radiation_source = (
        ROOT / "src/substrate_framework/sine_gordon_breather_radiation.py"
    ).read_text()
    compatibility = audit_numpy_trapezoid_compatibility(
        radiation_source,
        filename="sine_gordon_breather_radiation.py",
    )
    checks.check(
        "radiation imports the shared integration owner without legacy aliases",
        compatibility.legacy_references == 0
        and compatibility.current_references == 0
        and compatibility.eager_legacy_default_fallbacks == 0,
    )
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
