"""Exact primitives for the P249 exterior M5 clock candidate.

The physical field space is a symmetric three-dimensional spatial M5 tensor
``S`` in a constrained auxiliary Lorentz frame, together with a complex scalar
``psi``.  The clock axis is the first adapted spatial basis vector; the last
two basis vectors span its oriented tangent plane.  All functions are pure and
return exact SymPy objects.

This module implements the conditional action used by campaign P249.  It does
not attach a particle interpretation, a two-clock force, or gravity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import sympy as sp


AXIS_LOCK_STRENGTH = sp.Integer(1)
PHASE_LOCK_STRENGTH = sp.Integer(6)


def _symmetric_spatial_matrix(
    value: sp.MatrixBase | Sequence[Sequence[object]], name: str
) -> sp.Matrix:
    matrix = sp.Matrix(value)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix")
    if matrix != matrix.T:
        raise ValueError(f"{name} must be symmetric")
    return matrix


def clock_plane_basis() -> tuple[sp.ImmutableMatrix, ...]:
    """Return the normalized tangent-plane ``(A,H,K)`` clock basis.

    ``A`` generates rotations in the last two adapted spatial directions,
    ``H`` is the diagonal traceless tensor, and ``K=(A H-H A)/2``.
    """

    generator = sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    diagonal = sp.diag(0, 1, -1)
    partner = (generator * diagonal - diagonal * generator) / 2
    return tuple(
        sp.ImmutableMatrix(matrix.applyfunc(sp.factor))
        for matrix in (generator, diagonal, partner)
    )


def tangent_clock_tensor(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
) -> sp.ImmutableMatrix:
    """Return the traceless part of ``S`` on the clock tangent plane."""

    spatial = _symmetric_spatial_matrix(spatial_tensor, "spatial_tensor")
    tangent = sp.diag(0, 1, 1)
    projected = tangent * spatial * tangent
    tangent_trace = sp.trace(projected)
    clock = projected - tangent_trace * tangent / 2
    return sp.ImmutableMatrix(clock.applyfunc(sp.factor))


def phase_lock_target(psi_real: object, psi_imag: object) -> sp.ImmutableMatrix:
    """Return the weight-two tensor ``Q(psi)`` for ``psi=real+i*imag``."""

    real = sp.sympify(psi_real)
    imag = sp.sympify(psi_imag)
    _, diagonal, partner = clock_plane_basis()
    target = (real**2 - imag**2) * diagonal + 2 * real * imag * partner
    return sp.ImmutableMatrix(target.applyfunc(sp.factor))


def clock_tensor_norm_squared(
    tensor: sp.MatrixBase | Sequence[Sequence[object]],
) -> sp.Expr:
    """Return ``Tr(B^2)/2`` on the normalized clock-tensor doublet."""

    matrix = _symmetric_spatial_matrix(tensor, "tensor")
    return sp.factor(sp.trace(matrix**2) / 2)


def projected_m5_static_potential(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
) -> sp.Expr:
    """Return the exact beta=c=1 projected M5 Landau--de Gennes potential."""

    spatial = _symmetric_spatial_matrix(spatial_tensor, "spatial_tensor")
    trace_two = sp.trace(spatial**2)
    return sp.factor(
        -trace_two / 2
        - sp.trace(spatial**3)
        + trace_two**2
        + sp.Rational(1, 2)
    )


def aligned_axis_lock_potential(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
    strength: object = AXIS_LOCK_STRENGTH,
) -> sp.Expr:
    """Return the positive axis lock in the adapted ``N=e_1`` frame."""

    spatial = _symmetric_spatial_matrix(spatial_tensor, "spatial_tensor")
    coefficient = sp.sympify(strength)
    return sp.factor(
        coefficient * (sp.trace(spatial**2) - spatial[0, 0] ** 2)
    )


def scalar_amplitude_potential(amplitude: object) -> sp.Expr:
    """Return ``W(f)=3 f^2-4 f^3+2 f^4`` for a nonnegative amplitude."""

    value = sp.sympify(amplitude)
    return sp.factor(3 * value**2 - 4 * value**3 + 2 * value**4)


def phase_lock_potential(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
    psi_real: object,
    psi_imag: object,
    strength: object = PHASE_LOCK_STRENGTH,
) -> sp.Expr:
    """Return ``strength*||B-Q(psi)||^2`` exactly."""

    clock = sp.Matrix(tangent_clock_tensor(spatial_tensor))
    target = sp.Matrix(phase_lock_target(psi_real, psi_imag))
    return sp.factor(
        sp.sympify(strength) * clock_tensor_norm_squared(clock - target)
    )


def full_clock_potential(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
    psi_real: object,
    psi_imag: object,
) -> sp.Expr:
    """Compose the complete nonnegative P249 static potential."""

    real = sp.sympify(psi_real)
    imag = sp.sympify(psi_imag)
    amplitude = sp.sqrt(real**2 + imag**2)
    return sp.factor(
        projected_m5_static_potential(spatial_tensor)
        + aligned_axis_lock_potential(spatial_tensor)
        + phase_lock_potential(spatial_tensor, real, imag)
        + scalar_amplitude_potential(amplitude)
    )


def clock_inertia_density(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
    psi_real: object,
    psi_imag: object,
) -> sp.Expr:
    """Return the full canonical diagonal-S1 inertia density.

    The spatial tensor generator contains a weight-one axis--tangent shear
    doublet in addition to the weight-two tangent-traceless clock tensor.
    Writing the tensor contribution as ``Tr([A,S]^T[A,S])/2`` keeps every
    charged spatial channel in the Noether norm.
    """

    real = sp.sympify(psi_real)
    imag = sp.sympify(psi_imag)
    spatial = _symmetric_spatial_matrix(spatial_tensor, "spatial_tensor")
    generator, _, _ = clock_plane_basis()
    tensor_generator = sp.Matrix(generator) * spatial - spatial * sp.Matrix(
        generator
    )
    return sp.factor(
        real**2
        + imag**2
        + sp.trace(tensor_generator.T * tensor_generator) / 2
    )


def relative_equilibrium_velocity(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
    psi_real: object,
    psi_imag: object,
    frequency: object,
) -> tuple[sp.ImmutableMatrix, sp.Expr, sp.Expr]:
    """Return ``(dot S, dot Re(psi), dot Im(psi))`` on the clock orbit."""

    spatial = _symmetric_spatial_matrix(spatial_tensor, "spatial_tensor")
    real = sp.sympify(psi_real)
    imag = sp.sympify(psi_imag)
    omega = sp.sympify(frequency)
    generator, _, _ = clock_plane_basis()
    tensor_velocity = omega * (
        sp.Matrix(generator) * spatial - spatial * sp.Matrix(generator)
    )
    return (
        sp.ImmutableMatrix(tensor_velocity.applyfunc(sp.factor)),
        sp.factor(-omega * imag),
        sp.factor(omega * real),
    )


def canonical_velocity_energy_density(
    spatial_velocity: sp.MatrixBase | Sequence[Sequence[object]],
    psi_velocity_real: object,
    psi_velocity_imag: object,
) -> sp.Expr:
    """Return the canonical positive velocity contribution to the energy."""

    tensor_velocity = _symmetric_spatial_matrix(
        spatial_velocity, "spatial_velocity"
    )
    scalar_real = sp.sympify(psi_velocity_real)
    scalar_imag = sp.sympify(psi_velocity_imag)
    return sp.factor(
        sp.trace(tensor_velocity.T * tensor_velocity) / 4
        + (scalar_real**2 + scalar_imag**2) / 2
    )


def noether_charge_density(
    spatial_tensor: sp.MatrixBase | Sequence[Sequence[object]],
    spatial_velocity: sp.MatrixBase | Sequence[Sequence[object]],
    psi_real: object,
    psi_imag: object,
    psi_velocity_real: object,
    psi_velocity_imag: object,
) -> sp.Expr:
    """Return the full diagonal-clock canonical Noether charge density."""

    spatial = _symmetric_spatial_matrix(spatial_tensor, "spatial_tensor")
    tensor_velocity = _symmetric_spatial_matrix(
        spatial_velocity, "spatial_velocity"
    )
    real = sp.sympify(psi_real)
    imag = sp.sympify(psi_imag)
    scalar_velocity_real = sp.sympify(psi_velocity_real)
    scalar_velocity_imag = sp.sympify(psi_velocity_imag)
    generator, _, _ = clock_plane_basis()
    tensor_generator = sp.Matrix(generator) * spatial - spatial * sp.Matrix(
        generator
    )
    return sp.factor(
        real * scalar_velocity_imag
        - imag * scalar_velocity_real
        + sp.trace(tensor_velocity.T * tensor_generator) / 2
    )


def m5_radial_lower_bound(norm: object) -> sp.Expr:
    """Return the exact eigenvalue-norm lower bound for the M5 potential."""

    radius = sp.sympify(norm)
    return sp.factor(
        (radius - 1) ** 2 * (2 * radius**2 + 2 * radius + 1) / 2
    )


def unsplit_charge_ratio(amplitude: object) -> sp.Expr:
    """Return the exact scalar part of the ``B=0`` bulk ``J/K`` ratio."""

    value = sp.sympify(amplitude)
    return sp.factor(16 * (value - sp.Rational(1, 4)) ** 2 + 5)


def unsplit_full_lower_bound() -> sp.Integer:
    """Return the full ``B=0`` bulk ``J/K`` lower bound.

    The weight-one shear contribution has ratio four, while the scalar
    contribution is at least five; all remaining static terms are
    nonnegative.  Their inertia-weighted combination is therefore at least
    four.
    """

    return sp.Integer(4)


def fixed_charge_legendre_data(
    charge: object, inertia: object
) -> tuple[sp.Expr, ...]:
    """Return ``(E_rot, omega, d_I E, d_I^2 E, d_Q E)`` exactly."""

    charge_value = sp.sympify(charge)
    inertia_value = sp.sympify(inertia)
    rotational = charge_value**2 / (2 * inertia_value)
    frequency = charge_value / inertia_value
    return tuple(
        map(
            sp.factor,
            (
                rotational,
                frequency,
                -charge_value**2 / (2 * inertia_value**2),
                charge_value**2 / inertia_value**3,
                frequency,
            ),
        )
    )


def exterior_channel_pencil() -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]:
    """Return the exact spatial ``(K,G)`` pencil in ``(a,t,p,u,v,q)`` order."""

    stiffness = sp.diag(5, 10, 22, 4, 4, 22)
    kinetic = sp.diag(sp.Rational(1, 2), 1, 1, 1, 1, 1)
    return sp.ImmutableMatrix(stiffness), sp.ImmutableMatrix(kinetic)


def charged_continuum_edge_squared() -> sp.Rational:
    """Return ``min(6, 4, 22/2^2)=4`` over all charged channels."""

    return sp.Integer(4)


@dataclass(frozen=True)
class FullClockBindingWitness:
    """Exact constant-core witness for the P249 hylomorphy inequality."""

    scalar_amplitude: sp.Rational
    clock_split: sp.Rational
    potential_density: sp.Rational
    inertia_density: sp.Rational
    charge_ratio: sp.Rational
    charged_edge_squared: sp.Rational
    binding_margin: sp.Rational
    unsplit_binding_margin: sp.Rational


def full_clock_binding_witness() -> FullClockBindingWitness:
    """Return the exact ``f=1/2, B=diag(1/4,-1/4)`` plateau data."""

    return FullClockBindingWitness(
        scalar_amplitude=sp.Rational(1, 2),
        clock_split=sp.Rational(1, 4),
        potential_density=sp.Rational(45, 64),
        inertia_density=sp.Rational(1, 2),
        charge_ratio=sp.Rational(45, 16),
        charged_edge_squared=sp.Integer(4),
        binding_margin=sp.Rational(19, 16),
        unsplit_binding_margin=sp.Rational(19, 16),
    )


def fixed_charge_derrick_residual(
    gradient_energy: object, potential_energy: object, rotational_energy: object
) -> sp.Expr:
    """Return the three-dimensional fixed-charge scaling residual."""

    return sp.factor(
        sp.sympify(gradient_energy)
        + 3 * sp.sympify(potential_energy)
        - 3 * sp.sympify(rotational_energy)
    )
