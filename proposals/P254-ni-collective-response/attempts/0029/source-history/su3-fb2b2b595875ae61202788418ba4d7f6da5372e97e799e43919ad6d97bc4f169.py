"""Exact invariants of the standard fundamental SU(3) representation."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from typing import Any

import sympy as sp


@dataclass(frozen=True)
class SU3Invariants:
    """Representation invariants in the ``T_a=lambda_a/2`` convention."""

    dynkin_index: sp.Expr
    fundamental_casimir: sp.Expr
    adjoint_casimir: sp.Expr


@dataclass(frozen=True)
class SU3SymmetricTensorEvidence:
    """Exact symmetric tensor and anticommutator data in the standard basis."""

    tensor: sp.ImmutableDenseNDimArray
    fully_symmetric: bool
    anticommutator_residuals: tuple[sp.ImmutableMatrix, ...]
    embedded_su2_all_zero: bool
    outside_nonzero_witness: tuple[int, int, int, sp.Expr]


def fundamental_generators() -> tuple[sp.Matrix, ...]:
    """Return the eight Hermitian generators ``T_a=lambda_a/2``."""

    imaginary = sp.I
    root_three = sp.sqrt(3)
    gell_mann = (
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -imaginary, 0], [imaginary, 0, 0], [0, 0, 0]]),
        sp.diag(1, -1, 0),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, -imaginary], [0, 0, 0], [imaginary, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -imaginary], [0, imaginary, 0]]),
        sp.diag(1 / root_three, 1 / root_three, -2 / root_three),
    )
    return tuple(matrix / 2 for matrix in gell_mann)


@cache
def fundamental_commutant_basis() -> tuple[sp.Matrix, ...]:
    """Derive a basis for matrices commuting with all fundamental generators.

    The nullspace is computed from the explicit generators rather than assuming
    Schur's lemma or preselecting scalar matrices.
    """

    entries = sp.symbols("m0:9")
    candidate = sp.Matrix(3, 3, entries)
    equations = tuple(
        entry
        for generator in fundamental_generators()
        for entry in candidate * generator - generator * candidate
    )
    coefficient_matrix, _ = sp.linear_eq_to_matrix(equations, entries)
    return tuple(
        sp.Matrix(3, 3, list(vector))
        for vector in coefficient_matrix.nullspace()
    )


def center_element(power: int = 1) -> sp.Matrix:
    """Return ``omega**power I_3`` with ``omega=exp(2*pi*i/3)`` exactly."""

    if not isinstance(power, int):
        raise TypeError("power must be an integer")
    omega = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
    phases = (sp.Integer(1), omega, sp.conjugate(omega))
    return phases[power % 3] * sp.eye(3)


def center_elements() -> tuple[sp.Matrix, ...]:
    """Return the three elements of the fundamental SU(3) center."""

    return tuple(center_element(power) for power in range(3))


def triality_phase(triality: int, center_power: int = 1) -> sp.Expr:
    """Return the center phase for an abstract integer triality sector."""

    if not isinstance(triality, int):
        raise TypeError("triality must be an integer")
    if not isinstance(center_power, int):
        raise TypeError("center_power must be an integer")
    return sp.simplify(center_element(center_power)[0, 0] ** (triality % 3))


def center_conjugation(matrix: Any, power: int = 1) -> sp.Matrix:
    """Conjugate a 3-by-3 matrix by a fundamental center element."""

    value = sp.Matrix(matrix)
    if value.shape != (3, 3):
        raise ValueError("matrix must be 3 by 3")
    center = center_element(power)
    return sp.simplify(center * value * center.inv())


def _structure_constant(
    generators: tuple[sp.Matrix, ...], a: int, b: int, c: int
) -> sp.Expr:
    commutator = generators[a] * generators[b] - generators[b] * generators[a]
    return sp.simplify(-2 * sp.I * sp.trace(commutator * generators[c]))


def _validate_generator_indices(a: int, b: int, c: int) -> None:
    if any(
        isinstance(index, bool)
        or not isinstance(index, (int, sp.Integer))
        or index < 0
        or index >= 8
        for index in (a, b, c)
    ):
        raise IndexError("SU(3) generator index must lie in 0..7")


def structure_constant(a: int, b: int, c: int) -> sp.Expr:
    """Return ``f_abc`` for ``[T_a,T_b]=i*f_abc*T_c`` (zero-based)."""

    _validate_generator_indices(a, b, c)
    return _structure_constant(fundamental_generators(), a, b, c)


def symmetric_structure_constant(a: int, b: int, c: int) -> sp.Expr:
    r"""Return ``d_abc=2*Tr({T_a,T_b}*T_c)`` in the standard basis."""

    _validate_generator_indices(a, b, c)
    generators = fundamental_generators()
    anticommutator = generators[a] * generators[b] + generators[b] * generators[a]
    return sp.simplify(2 * sp.trace(anticommutator * generators[c]))


@cache
def symmetric_tensor_evidence() -> SU3SymmetricTensorEvidence:
    r"""Derive the full fundamental ``d_abc`` tensor and its exact identities.

    The anticommutator convention is
    ``{T_a,T_b}=delta_ab*I_3/3+sum_c d_abc*T_c``.  The standard generators
    ``T_1,T_2,T_3`` form an embedded SU(2) on which every restricted ``d``
    component vanishes, while components outside that restriction are nonzero.
    This is representation algebra only and supplies no physical color sector.
    """

    generators = fundamental_generators()
    tensor = sp.ImmutableDenseNDimArray(
        [
            symmetric_structure_constant(a, b, c)
            for a in range(8)
            for b in range(8)
            for c in range(8)
        ],
        (8, 8, 8),
    )
    fully_symmetric = all(
        sp.simplify(tensor[a, b, c] - tensor[b, a, c]) == 0
        and sp.simplify(tensor[a, b, c] - tensor[a, c, b]) == 0
        for a in range(8)
        for b in range(8)
        for c in range(8)
    )
    identity = sp.eye(3)
    residuals = tuple(
        sp.ImmutableMatrix(
            (
                generators[a] * generators[b]
                + generators[b] * generators[a]
                - sp.Rational(1, 3) * (1 if a == b else 0) * identity
                - sum(
                    (tensor[a, b, c] * generators[c] for c in range(8)),
                    sp.zeros(3),
                )
            ).applyfunc(sp.simplify)
        )
        for a in range(8)
        for b in range(8)
    )
    embedded_zero = all(
        tensor[a, b, c] == 0
        for a in range(3)
        for b in range(3)
        for c in range(3)
    )
    witness = (0, 0, 7, sp.simplify(tensor[0, 0, 7]))
    return SU3SymmetricTensorEvidence(
        tensor=tensor,
        fully_symmetric=fully_symmetric,
        anticommutator_residuals=residuals,
        embedded_su2_all_zero=embedded_zero,
        outside_nonzero_witness=witness,
    )


@cache
def invariants() -> SU3Invariants:
    """Compute ``T_F``, ``C_F``, and ``C_A`` from explicit matrices."""

    generators = fundamental_generators()
    dynkin = sp.trace(generators[0] * generators[0])
    fundamental = sp.simplify(
        sum((generator * generator for generator in generators), sp.zeros(3))
    )
    adjoint_generators = tuple(
        sp.Matrix(
            8,
            8,
            lambda b, c: -sp.I * _structure_constant(generators, a, b, c),
        )
        for a in range(8)
    )
    adjoint = sp.simplify(
        sum((generator * generator for generator in adjoint_generators), sp.zeros(8))
    )
    if fundamental != fundamental[0, 0] * sp.eye(3):
        raise ValueError("fundamental Casimir is not scalar")
    if adjoint != adjoint[0, 0] * sp.eye(8):
        raise ValueError("adjoint Casimir is not scalar")
    return SU3Invariants(dynkin, fundamental[0, 0], adjoint[0, 0])


def conditional_one_loop_coefficient(
    flavor_count: Any,
    gauge_loop_weight: Any,
    matter_loop_weight: Any,
) -> sp.Expr:
    """Compose derived SU(3) invariants with declared one-loop weights."""

    flavors = sp.sympify(flavor_count)
    if flavors.is_number and (
        flavors.is_integer is not True or flavors.is_nonnegative is not True
    ):
        raise ValueError("flavor_count must be a nonnegative integer")
    gauge = sp.sympify(gauge_loop_weight)
    matter = sp.sympify(matter_loop_weight)
    values = invariants()
    return sp.simplify(
        gauge * values.adjoint_casimir
        - matter * values.dynkin_index * flavors
    )
