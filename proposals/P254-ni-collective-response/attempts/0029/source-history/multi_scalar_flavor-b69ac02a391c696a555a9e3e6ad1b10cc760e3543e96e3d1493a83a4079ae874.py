"""Exact mass-basis ledgers for a finite family of scalar coupling matrices.

The inputs in this module are declared finite matrices and scalar weights.  A
physical Yukawa interaction, scalar field content, vacuum expectation values,
generation map, decay rate, or experimental flavor bound must be supplied by
separate accepted claims.  Small matrix entries are not interpreted as
natural flavor conservation here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp


def _exact_expression(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    return expression


def _exact_square_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(value)
    if matrix.rows == 0 or matrix.rows != matrix.cols:
        raise ValueError(f"{name} must be non-empty and square")
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError(f"{name} must be exact rather than floating")
    return matrix


def _is_zero_matrix(matrix: sp.MatrixBase) -> bool:
    return all(sp.simplify(entry) == 0 for entry in matrix)


def _exact_unitary(value: Any, name: str, dimension: int) -> sp.ImmutableMatrix:
    matrix = _exact_square_matrix(value, name)
    if matrix.rows != dimension:
        raise ValueError(f"{name} has the wrong dimension")
    identity = sp.eye(dimension)
    if not _is_zero_matrix(matrix.adjoint() * matrix - identity) or not _is_zero_matrix(
        matrix * matrix.adjoint() - identity
    ):
        raise ValueError(f"{name} must be exactly unitary")
    return matrix


def off_diagonal_part(matrix: Any) -> sp.ImmutableMatrix:
    """Return a square matrix with its diagonal entries set to zero exactly."""

    value = _exact_square_matrix(matrix, "matrix")
    return sp.ImmutableMatrix(
        value.rows,
        value.cols,
        lambda row, column: (
            sp.Integer(0) if row == column else sp.simplify(value[row, column])
        ),
    )


@dataclass(frozen=True)
class MultiScalarMassBasisLedger:
    """Exact reconstruction of several couplings in one biunitary basis."""

    dimension: int
    scalar_weights: tuple[sp.Expr, ...]
    coupling_matrices: tuple[sp.ImmutableMatrix, ...]
    mass_matrix: sp.ImmutableMatrix
    left_basis: sp.ImmutableMatrix
    right_basis: sp.ImmutableMatrix
    diagonal_mass_matrix: sp.ImmutableMatrix
    mass_basis_couplings: tuple[sp.ImmutableMatrix, ...]
    off_diagonal_couplings: tuple[sp.ImmutableMatrix, ...]
    reconstructed_diagonal_mass_matrix: sp.ImmutableMatrix
    reconstruction_residual: sp.ImmutableMatrix
    all_couplings_diagonal: bool


def multi_scalar_mass_basis_ledger(
    coupling_matrices: Iterable[Any],
    scalar_weights: Iterable[Any],
    left_basis: Any,
    right_basis: Any,
) -> MultiScalarMassBasisLedger:
    r"""Transform a declared matrix family into one biunitary mass basis.

    Given square matrices ``Y_a`` and weights ``v_a``, this function constructs
    ``M=sum_a v_a*Y_a``.  The supplied unitary bases must make
    ``D=U_L.H*M*U_R`` exactly diagonal.  It then returns every
    ``Gamma_a=U_L.H*Y_a*U_R`` and verifies
    ``D=sum_a v_a*Gamma_a`` exactly.

    A diagonal sum does not make its summands diagonal: off-diagonal entries
    may cancel between scalar couplings.  Conditional on a separately declared
    multi-scalar interaction, flavor-diagonal neutral couplings in this fixed
    mass basis therefore require every returned ``Gamma_a`` to be diagonal.
    Repeated mass singular values can leave additional admissible basis choices,
    so callers must retain that degeneracy audit.
    """

    matrices = tuple(
        _exact_square_matrix(matrix, f"coupling_matrices[{index}]")
        for index, matrix in enumerate(coupling_matrices)
    )
    if not matrices:
        raise ValueError("at least one coupling matrix is required")
    dimension = matrices[0].rows
    if any(matrix.shape != (dimension, dimension) for matrix in matrices):
        raise ValueError("coupling matrices must have one common square shape")
    weights = tuple(
        _exact_expression(weight, f"scalar_weights[{index}]")
        for index, weight in enumerate(scalar_weights)
    )
    if len(weights) != len(matrices):
        raise ValueError("one scalar weight is required per coupling matrix")
    left = _exact_unitary(left_basis, "left_basis", dimension)
    right = _exact_unitary(right_basis, "right_basis", dimension)

    mass = sp.zeros(dimension)
    for weight, matrix in zip(weights, matrices, strict=True):
        mass += weight * matrix
    mass = sp.ImmutableMatrix(mass.applyfunc(sp.simplify))
    diagonal = sp.ImmutableMatrix(
        (left.adjoint() * mass * right).applyfunc(sp.simplify)
    )
    if not _is_zero_matrix(off_diagonal_part(diagonal)):
        raise ValueError("the supplied bases do not diagonalize the mass matrix")

    transformed = tuple(
        sp.ImmutableMatrix((left.adjoint() * matrix * right).applyfunc(sp.simplify))
        for matrix in matrices
    )
    off_diagonal = tuple(off_diagonal_part(matrix) for matrix in transformed)
    reconstructed = sp.zeros(dimension)
    for weight, matrix in zip(weights, transformed, strict=True):
        reconstructed += weight * matrix
    reconstructed = sp.ImmutableMatrix(reconstructed.applyfunc(sp.simplify))
    residual = sp.ImmutableMatrix((reconstructed - diagonal).applyfunc(sp.simplify))
    if not _is_zero_matrix(residual):
        raise ArithmeticError("mass-basis reconstruction failed")
    return MultiScalarMassBasisLedger(
        dimension=dimension,
        scalar_weights=weights,
        coupling_matrices=matrices,
        mass_matrix=mass,
        left_basis=left,
        right_basis=right,
        diagonal_mass_matrix=diagonal,
        mass_basis_couplings=transformed,
        off_diagonal_couplings=off_diagonal,
        reconstructed_diagonal_mass_matrix=reconstructed,
        reconstruction_residual=residual,
        all_couplings_diagonal=all(_is_zero_matrix(matrix) for matrix in off_diagonal),
    )


def takagi_multi_scalar_mass_basis_ledger(
    coupling_matrices: Iterable[Any],
    scalar_weights: Iterable[Any],
    takagi_basis: Any,
) -> MultiScalarMassBasisLedger:
    r"""Use the correct right basis for a declared complex-symmetric mass matrix.

    If ``M=U*D*U.T``, then ``U.H*M*conjugate(U)=D``.  Each source coupling must
    therefore transform as ``U.H*Y_a*conjugate(U)``.  Reusing ``U`` instead of
    ``conjugate(U)`` on the right is not the Takagi mass-basis transformation.
    """

    matrices = tuple(coupling_matrices)
    weights = tuple(scalar_weights)
    if not matrices:
        raise ValueError("at least one coupling matrix is required")
    exact_matrices = tuple(
        _exact_square_matrix(matrix, f"coupling_matrices[{index}]")
        for index, matrix in enumerate(matrices)
    )
    if any(matrix.shape != exact_matrices[0].shape for matrix in exact_matrices):
        raise ValueError("coupling matrices must have one common square shape")
    if len(weights) != len(exact_matrices):
        raise ValueError("one scalar weight is required per coupling matrix")
    exact_weights = tuple(
        _exact_expression(weight, f"scalar_weights[{index}]")
        for index, weight in enumerate(weights)
    )
    mass = sp.zeros(exact_matrices[0].rows)
    for weight, matrix in zip(exact_weights, exact_matrices, strict=True):
        mass += weight * matrix
    mass = sp.ImmutableMatrix(mass.applyfunc(sp.simplify))
    if not _is_zero_matrix(mass - mass.T):
        raise ValueError("the declared Takagi mass matrix must be symmetric")
    unitary = _exact_unitary(takagi_basis, "takagi_basis", mass.rows)
    return multi_scalar_mass_basis_ledger(
        exact_matrices,
        exact_weights,
        unitary,
        unitary.conjugate(),
    )
