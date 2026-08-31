"""Tests for the canonical P249 full-spatial M5 clock primitives."""

from __future__ import annotations

import sympy as sp

from substrate_framework.m5_exterior_clock import (
    aligned_axis_lock_potential,
    canonical_velocity_energy_density,
    charged_continuum_edge_squared,
    clock_inertia_density,
    clock_plane_basis,
    clock_tensor_norm_squared,
    exterior_channel_pencil,
    fixed_charge_derrick_residual,
    fixed_charge_legendre_data,
    full_clock_binding_witness,
    full_clock_potential,
    m5_radial_lower_bound,
    noether_charge_density,
    phase_lock_potential,
    phase_lock_target,
    projected_m5_static_potential,
    relative_equilibrium_velocity,
    scalar_amplitude_potential,
    tangent_clock_tensor,
    unsplit_charge_ratio,
    unsplit_full_lower_bound,
)


def test_clock_plane_basis_and_phase_target_are_equivariant() -> None:
    angle, real, imag = sp.symbols("angle real imag", real=True)
    generator, diagonal, partner = clock_plane_basis()
    rotation = sp.eye(3)
    rotation[1, 1] = sp.cos(angle)
    rotation[1, 2] = -sp.sin(angle)
    rotation[2, 1] = sp.sin(angle)
    rotation[2, 2] = sp.cos(angle)

    assert generator**2 == -sp.diag(0, 1, 1)
    assert clock_tensor_norm_squared(diagonal) == 1
    assert clock_tensor_norm_squared(partner) == 1
    assert sp.trace(diagonal * partner) == 0

    rotated_real = real * sp.cos(angle) - imag * sp.sin(angle)
    rotated_imag = real * sp.sin(angle) + imag * sp.cos(angle)
    target_after = sp.Matrix(phase_lock_target(rotated_real, rotated_imag))
    target_before = sp.Matrix(phase_lock_target(real, imag))
    residual = (target_after - rotation * target_before * rotation.T).applyfunc(
        sp.trigsimp
    )
    assert residual == sp.zeros(3)


def test_static_potential_and_axis_lock_have_exact_vacuum_and_split_values() -> None:
    split = sp.symbols("split", real=True)
    vacuum = sp.diag(1, 0, 0)
    clock_core = sp.diag(1, split, -split)

    assert projected_m5_static_potential(vacuum) == 0
    assert aligned_axis_lock_potential(vacuum) == 0
    assert sp.expand(
        projected_m5_static_potential(clock_core) - (3 * split**2 + 4 * split**4)
    ) == 0
    assert aligned_axis_lock_potential(clock_core) == 2 * split**2
    assert tangent_clock_tensor(clock_core) == sp.diag(0, split, -split)


def test_phase_lock_and_scalar_potential_are_positive_exact_squares() -> None:
    amplitude, split = sp.symbols("amplitude split", nonnegative=True, real=True)
    spatial = sp.diag(1, split, -split)

    assert scalar_amplitude_potential(amplitude) == sp.factor(
        amplitude**2 * (2 * (amplitude - 1) ** 2 + 1)
    )
    assert sp.expand(
        phase_lock_potential(spatial, amplitude, 0)
        - 6 * (split - amplitude**2) ** 2
    ) == 0
    assert clock_inertia_density(spatial, amplitude, 0) == amplitude**2 + 4 * split**2


def test_full_binding_witness_is_computed_independently_from_action() -> None:
    witness = full_clock_binding_witness()
    spatial = sp.diag(1, witness.clock_split, -witness.clock_split)
    potential = full_clock_potential(spatial, witness.scalar_amplitude, 0)
    inertia = clock_inertia_density(spatial, witness.scalar_amplitude, 0)

    assert potential == sp.Rational(45, 64)
    assert inertia == sp.Rational(1, 2)
    assert sp.factor(2 * potential / inertia) == sp.Rational(45, 16)
    assert witness.charged_edge_squared - 2 * potential / inertia == sp.Rational(
        19, 16
    )
    assert unsplit_charge_ratio(sp.Rational(1, 4)) == 5
    assert unsplit_full_lower_bound() == 4
    assert charged_continuum_edge_squared() - 2 * potential / inertia == sp.Rational(
        19, 16
    )


def test_exterior_pencil_has_all_six_positive_physical_channels() -> None:
    stiffness, kinetic = exterior_channel_pencil()
    generalized = kinetic.inv() * stiffness

    assert kinetic.det() == sp.Rational(1, 2)
    assert generalized == sp.diag(10, 10, 22, 4, 4, 22)
    assert charged_continuum_edge_squared() == 4


def test_full_noether_inertia_includes_weight_one_spatial_shears() -> None:
    u, v, p, q = sp.symbols("u v p q", real=True)
    spatial = sp.Matrix([[1, u, v], [u, p, q], [v, q, -p]])
    generator = sp.Matrix(clock_plane_basis()[0])
    tangent = generator * spatial - spatial * generator

    assert tangent == sp.Matrix(
        [[0, -v, u], [-v, -2 * q, 2 * p], [u, 2 * p, 2 * q]]
    )
    assert clock_inertia_density(spatial, 0, 0) == (
        u**2 + v**2 + 4 * p**2 + 4 * q**2
    )


def test_canonical_phase_space_reduces_exactly_on_relative_equilibrium() -> None:
    u, v, p, q, real, imag, omega = sp.symbols(
        "u v p q real imag omega", real=True
    )
    spatial = sp.Matrix([[1, u, v], [u, p, q], [v, q, -p]])
    tensor_velocity, scalar_velocity_real, scalar_velocity_imag = (
        relative_equilibrium_velocity(spatial, real, imag, omega)
    )
    inertia = clock_inertia_density(spatial, real, imag)

    assert noether_charge_density(
        spatial,
        tensor_velocity,
        real,
        imag,
        scalar_velocity_real,
        scalar_velocity_imag,
    ) == omega * inertia
    assert canonical_velocity_energy_density(
        tensor_velocity, scalar_velocity_real, scalar_velocity_imag
    ) == omega**2 * inertia / 2


def test_legendre_and_derrick_data_include_complete_fixed_charge_terms() -> None:
    charge, inertia = sp.symbols("charge inertia", positive=True, real=True)
    rotational, frequency, first, second, envelope = fixed_charge_legendre_data(
        charge, inertia
    )

    assert rotational == charge**2 / (2 * inertia)
    assert frequency == charge / inertia
    assert first == -frequency**2 / 2
    assert second == frequency**2 / inertia
    assert envelope == frequency
    assert fixed_charge_derrick_residual(2, 3, 4) == -1


def test_m5_lower_bound_and_input_validation() -> None:
    radius = sp.symbols("radius", nonnegative=True, real=True)
    assert m5_radial_lower_bound(radius) == sp.factor(
        (radius - 1) ** 2 * (2 * radius**2 + 2 * radius + 1) / 2
    )

    try:
        tangent_clock_tensor(sp.eye(2))
    except ValueError as error:
        assert "3x3" in str(error)
    else:
        raise AssertionError("a non-3x3 spatial tensor must be rejected")

    try:
        projected_m5_static_potential(sp.Matrix([[1, 1, 0], [0, 0, 0], [0, 0, 0]]))
    except ValueError as error:
        assert "symmetric" in str(error)
    else:
        raise AssertionError("a nonsymmetric spatial tensor must be rejected")
