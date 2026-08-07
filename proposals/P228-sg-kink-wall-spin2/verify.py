#!/usr/bin/env python3
"""Primary verifier for proposed P228 claims C-EGR-009 through C-EGR-012."""

from __future__ import annotations

import inspect
from pathlib import Path

import sympy as sp
import yaml

from substrate_framework.dimensional_sine_gordon import (
    dimensional_sine_gordon_coefficients,
)
from substrate_framework.emergent_spin2 import (
    collective_tensor_constraint_matrix,
    fierz_pauli_gauge_ledger,
    fierz_pauli_gauge_residual,
    sine_gordon_spin_two_coupling_ledger,
    spin_two_mode_ledger,
    tensor_mode_count_from_constraints,
)
from substrate_framework.sine_gordon_spin2_predictions import (
    controlled_static_breather_prediction,
    sine_gordon_spin2_radiation_evidence,
)
from substrate_framework.sine_gordon_wall_network import (
    SineGordonWallNetworkScales,
    icosahedral_wall_directions,
    icosahedral_wall_frame_ledger,
    local_sine_gordon_wall_network,
    sine_gordon_wall_network_scales,
    spatial_tensor_from_wall_channels,
    wall_channels_from_spatial_tensor,
)
from substrate_framework.verification import CheckLedger


ROOT = Path(__file__).resolve().parents[2]
PROPOSAL = ROOT / "proposals/P228-sg-kink-wall-spin2"


def _direction_projection(directions: tuple[sp.ImmutableMatrix, ...]) -> sp.Matrix:
    return sp.Matrix(
        [
            [
                direction[0] ** 2,
                direction[1] ** 2,
                direction[2] ** 2,
                2 * direction[0] * direction[1],
                2 * direction[0] * direction[2],
                2 * direction[1] * direction[2],
            ]
            for direction in directions
        ]
    )


def main() -> int:
    checks = CheckLedger("P228-C-EGR-009-012")
    manifest = yaml.safe_load((PROPOSAL / "proposal.yaml").read_text())
    checks.check(
        "proposal reserves exactly the four new emergent-gravity claims",
        manifest["claims_proposed"]
        == ["C-EGR-009", "C-EGR-010", "C-EGR-011", "C-EGR-012"],
    )
    comparison = yaml.safe_load(
        (PROPOSAL / "evidence/candidate-comparison.yaml").read_text()
    )
    checks.check(
        "candidate A was selected structurally before empirical comparison",
        comparison["decisions"]["A"]["decision"] == "selected_for_implementation"
        and all(
            comparison["decisions"][candidate]["decision"].startswith("rejected")
            for candidate in ("B", "C", "D")
        )
        and comparison["comparators_opened"] == [],
    )

    time, x, y, z = sp.symbols("t x y z", real=True)
    microscopic_fields = tuple(
        sp.Function(f"u_{index}")(time, x, y, z) for index in range(6)
    )
    microscopic = local_sine_gordon_wall_network(
        microscopic_fields,
        (x, y, z),
        time,
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    checks.check(
        "the microscopic model has explicit local transverse dynamics",
        microscopic.spatial_gradient_hessian == -18 * sp.eye(18)
        and microscopic.spatial_gradient_hessian.rank() == 18
        and all(
            microscopic.lagrangian_density.has(sp.diff(field, coordinate) ** 2)
            for field in microscopic_fields
            for coordinate in (x, y, z)
        ),
    )

    frame = icosahedral_wall_frame_ledger()
    checks.check(
        "the six wall directions are an invertible isotropic tight frame",
        frame.projection_rank == 6
        and frame.direction_second_moment == 2 * sp.eye(3)
        and frame.tight_frame_identity_residual == 0
        and frame.traceless_quadratic_weight == sp.Rational(4, 5),
    )
    tensor = sp.ImmutableMatrix([[1, 2, 3], [2, 5, 7], [3, 7, 11]])
    checks.check(
        "microscopic wall channels reconstruct a unique collective tensor",
        spatial_tensor_from_wall_channels(
            wall_channels_from_spatial_tensor(tensor)
        )
        == tensor,
    )
    directions = icosahedral_wall_directions()
    repeated_axis_projection = _direction_projection(
        directions[:-1] + (directions[0],)
    )
    checks.check(
        "repeating one wall family destroys the tensor reconstruction",
        repeated_axis_projection.rank() == 5,
    )

    lam, tension, onsite = sp.symbols("lambda T mu", positive=True)
    coefficients = dimensional_sine_gordon_coefficients(lam, tension, onsite)
    wall = sine_gordon_wall_network_scales(coefficients)
    checks.check(
        "the exact kink and one-cell lift fix the tensor normalization",
        wall.signal_speed == sp.sqrt(tension / lam)
        and wall.profile_length == sp.sqrt(tension / onsite)
        and wall.kink_gradient_integral == 8 * sp.sqrt(onsite / tension)
        and wall.wall_inertia_density == 8 * lam * onsite / tension
        and wall.wall_modulus == 8 * onsite
        and wall.spin2_time_inertia == 32 * lam * onsite / (5 * tension)
        and wall.spin2_spacetime_normalization == 32 * onsite / 5,
    )
    checks.check(
        "the wall-scale API contains no gravitational coupling shortcut",
        tuple(inspect.signature(sine_gordon_wall_network_scales).parameters)
        == ("coefficients",)
        and "einstein_coupling" not in SineGordonWallNetworkScales.__dataclass_fields__
        and "newton_constant" not in SineGordonWallNetworkScales.__dataclass_fields__,
    )

    gauge = fierz_pauli_gauge_ledger()
    checks.check(
        "the general quadratic gauge solve selects one Fierz-Pauli ray",
        gauge.sampled_constraint_rank == 3
        and gauge.allowed_coefficient_dimension == 1
        and gauge.normalized_coefficient_ray == (1, -2, 2, -1)
        and gauge.symbolic_gauge_residual == sp.zeros(10, 4),
    )
    checks.check(
        "a single coefficient mutation breaks the gauge kernel",
        fierz_pauli_gauge_residual((1, -2, 2, -1), (1, 2, 3, 5))
        == sp.zeros(10, 4)
        and fierz_pauli_gauge_residual(
            (1, -2, 2, sp.Rational(-9, 10)),
            (1, 2, 3, 5),
        ).rank()
        > 0,
    )

    mode = spin_two_mode_ledger((1, 2, 3), wall.signal_speed, wall.profile_length)
    checks.check(
        "volume and force balance derive two positive physical modes",
        mode.constraint_rank == 4
        and mode.physical_mode_count == 2
        and mode.projected_kinetic_rank == 2
        and mode.projected_frobenius_metric.is_positive_definite is True,
    )
    checks.check(
        "removing volume balance exposes an additional scalar mode",
        tensor_mode_count_from_constraints(
            collective_tensor_constraint_matrix((1, 2, 3))
        )
        == 2
        and tensor_mode_count_from_constraints(
            collective_tensor_constraint_matrix(
                (1, 2, 3),
                include_volume_constraint=False,
            )
        )
        == 3,
    )
    spacing = sp.Symbol("ell", positive=True)
    lattice_mode = spin_two_mode_ledger((3, 0, 0), 2, spacing)
    checks.check(
        "the local registry has a shared relativistic cone and visible cutoff term",
        lattice_mode.continuum_angular_frequency_squared == 36
        and sp.limit(lattice_mode.lattice_angular_frequency_squared, spacing, 0)
        == 36
        and lattice_mode.leading_lattice_correction == -27 * spacing**2
        and lattice_mode.relative_lattice_correction
        == -sp.Rational(3, 4) * spacing**2,
    )

    gravity = sine_gordon_spin_two_coupling_ledger(coefficients)
    checks.check(
        "gauge completion fixes kappa and G without a new parameter",
        gravity.normalized_fierz_pauli_coefficients
        == (
            -16 * onsite / 5,
            32 * onsite / 5,
            -32 * onsite / 5,
            16 * onsite / 5,
        )
        and gravity.einstein_coupling == 5 / (128 * onsite)
        and gravity.trace_reversed_source_coefficient == -5 / (64 * onsite)
        and sp.simplify(
            gravity.newton_constant
            - 5 * tension**2 / (1024 * sp.pi * lam**2 * onsite)
        )
        == 0
        and tuple(
            inspect.signature(sine_gordon_spin_two_coupling_ledger).parameters
        )
        == ("coefficients",),
    )
    scaled_gravity = sine_gordon_spin_two_coupling_ledger(
        dimensional_sine_gordon_coefficients(7 * lam, 7 * tension, 7 * onsite)
    )
    checks.check(
        "the common sine-Gordon scale is load bearing rather than hidden",
        scaled_gravity.wall_scales.signal_speed == gravity.wall_scales.signal_speed
        and scaled_gravity.wall_scales.profile_length
        == gravity.wall_scales.profile_length
        and scaled_gravity.newton_constant == gravity.newton_constant / 7,
    )

    benchmark_coefficients = dimensional_sine_gordon_coefficients(1, 1, 1)
    inverse_width = sp.Rational(3, 100)
    frequency = sp.sqrt(1 - inverse_width**2)
    static = controlled_static_breather_prediction(
        frequency,
        benchmark_coefficients,
        10,
    )
    checks.check(
        "the exact source moment controls the exterior static prediction",
        static.source_energy == sp.Rational(12, 25)
        and static.source_mass == sp.Rational(12, 25)
        and static.observer_radius == sp.Rational(1000, 3)
        and static.maximum_relative_profile_error < sp.Rational(1, 200)
        and static.exact_potential_lower_bound
        < static.exact_potential_upper_bound
        < 0,
    )
    checks.check(
        "doubling observer radius reduces the static profile bound fourfold",
        controlled_static_breather_prediction(
            frequency,
            benchmark_coefficients,
            20,
        ).maximum_relative_profile_error
        == static.maximum_relative_profile_error / 4,
    )

    coarse = sine_gordon_spin2_radiation_evidence(
        frequency,
        benchmark_coefficients,
        spatial_points=513,
        time_points=65,
        direction_points=257,
        harmonic_count=3,
    )
    refined = sine_gordon_spin2_radiation_evidence(
        frequency,
        benchmark_coefficients,
        spatial_points=1025,
        time_points=129,
        direction_points=501,
        harmonic_count=4,
    )
    checks.check(
        "full-retardation power is positive and refinement stable",
        refined.power_over_onsite_speed > 0
        and abs(
            refined.power_over_onsite_speed - coarse.power_over_onsite_speed
        )
        < 7.9e-13
        and abs(refined.power_over_onsite_speed - 1.1232794952659703e-5)
        < 3.0e-18,
    )
    checks.check(
        "radiation frequency and one-cycle backreaction gates remain explicit",
        refined.source.radiation_frequency_size_ratio > 60
        and refined.field_cycle_energy_loss_fraction < 1.5e-4
        and refined.source.relative_leading_difference < 7.3e-4,
    )
    checks.check(
        "static and radiative consumers use the identical derived coupling",
        refined.newton_constant == static.newton_constant
        == sp.Rational(5, 1024) / sp.pi
        and refined.signal_speed == static.signal_speed == 1,
    )
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
