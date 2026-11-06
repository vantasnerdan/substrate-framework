"""Diagonal-rephasing invariants of finite unitary matrices.

This module supplies matrix and group-action algebra only.  It does not assign
flavor, CKM, or physical CP interpretations to a unitary matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

from .matrix_decompositions import unitarity_residual


ComplexMatrix = NDArray[np.complex128]


def _unitary_matrix(matrix: Any, tolerance: float) -> ComplexMatrix:
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    value = np.asarray(matrix, dtype=np.complex128)
    if value.ndim != 2 or value.shape[0] == 0 or value.shape[0] != value.shape[1]:
        raise ValueError("matrix must be non-empty and square")
    if not np.all(np.isfinite(value)):
        raise ValueError("matrix entries must be finite")
    if unitarity_residual(value) > tolerance:
        raise ValueError("matrix must be unitary within tolerance")
    return value


@dataclass(frozen=True)
class RephasingCounts:
    """Generic dimensions for ``U(N)`` modulo two diagonal phase bases."""

    size: int
    unitary_parameters: int
    orthogonal_angles: int
    phase_parameters: int
    torus_parameters: int
    generic_kernel: int
    generic_orbit: int
    generic_quotient: int
    irreducible_phases: int


def generic_rephasing_counts(size: int) -> RephasingCounts:
    """Return generic parameter counts for a positive integer matrix size.

    The left/right diagonal torus has dimension ``2N``.  On the connected
    nonzero-support stratum its sole stabilizer is the common phase, so the
    orbit has dimension ``2N-1``.  Matrices with disconnected support can have
    larger stabilizers; use :func:`support_stabilizer_dimension` for them.
    """

    if isinstance(size, bool) or not isinstance(size, Integral) or size < 1:
        raise ValueError("size must be a positive integer")
    size = int(size)
    unitary = size**2
    angles = size * (size - 1) // 2
    phases = unitary - angles
    torus = 2 * size
    kernel = 1
    orbit = torus - kernel
    quotient = unitary - orbit
    irreducible = quotient - angles
    return RephasingCounts(
        size=size,
        unitary_parameters=unitary,
        orthogonal_angles=angles,
        phase_parameters=phases,
        torus_parameters=torus,
        generic_kernel=kernel,
        generic_orbit=orbit,
        generic_quotient=quotient,
        irreducible_phases=irreducible,
    )


def support_stabilizer_dimension(
    matrix: Any,
    *,
    tolerance: float = 1e-12,
) -> int:
    """Return the diagonal-rephasing stabilizer dimension from matrix support.

    Rows and columns form the vertices of a bipartite graph and every nonzero
    entry supplies an edge.  The stabilizer dimension is the number of
    connected components.  A unitary matrix has no isolated row or column.
    """

    value = _unitary_matrix(matrix, tolerance)
    size = value.shape[0]
    adjacency: list[set[int]] = [set() for _ in range(2 * size)]
    for row in range(size):
        for column in range(size):
            if abs(value[row, column]) > tolerance:
                right_vertex = size + column
                adjacency[row].add(right_vertex)
                adjacency[right_vertex].add(row)

    remaining = set(range(2 * size))
    components = 0
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            vertex = stack.pop()
            unseen = adjacency[vertex] & remaining
            remaining.difference_update(unseen)
            stack.extend(unseen)
    return components


def rephasing_orbit_dimension(
    matrix: Any,
    *,
    tolerance: float = 1e-12,
) -> int:
    """Return ``2N - stabilizer_dimension`` for the matrix support stratum."""

    value = _unitary_matrix(matrix, tolerance)
    return 2 * value.shape[0] - support_stabilizer_dimension(
        value, tolerance=tolerance
    )


def rephase_unitary(
    matrix: Any,
    left_phases: Sequence[float],
    right_phases: Sequence[float],
    *,
    tolerance: float = 1e-12,
) -> ComplexMatrix:
    """Apply ``D_left @ matrix @ D_right.H`` with real phase vectors."""

    value = _unitary_matrix(matrix, tolerance)
    left = np.asarray(left_phases, dtype=np.float64)
    right = np.asarray(right_phases, dtype=np.float64)
    if left.shape != (value.shape[0],) or right.shape != (value.shape[1],):
        raise ValueError("phase vectors must match the matrix size")
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ValueError("phases must be finite")
    left_diagonal = np.diag(np.exp(1j * left))
    right_diagonal = np.diag(np.exp(1j * right))
    return left_diagonal @ value @ right_diagonal.conj().T


def invariant_quartet(
    matrix: Any,
    first_row: int,
    second_row: int,
    first_column: int,
    second_column: int,
    *,
    tolerance: float = 1e-12,
) -> complex:
    """Return ``V_ij V_kl conjugate(V_il) conjugate(V_kj)``."""

    value = _unitary_matrix(matrix, tolerance)
    size = value.shape[0]
    indices = (first_row, second_row, first_column, second_column)
    if any(isinstance(index, bool) or not isinstance(index, Integral) for index in indices):
        raise TypeError("quartet indices must be integers")
    i, k, j, ell = (int(index) for index in indices)
    if any(index < 0 or index >= size for index in (i, k, j, ell)):
        raise IndexError("quartet index is outside the matrix")
    return complex(
        value[i, j]
        * value[k, ell]
        * np.conj(value[i, ell])
        * np.conj(value[k, j])
    )


def standard_three_angle_unitary(
    theta12: float,
    theta13: float,
    theta23: float,
    phase: float,
) -> ComplexMatrix:
    """Return the declared ``R23 @ R13(phase) @ R12`` unitary matrix."""

    values = np.asarray([theta12, theta13, theta23, phase], dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError("angles and phase must be finite")
    t12, t13, t23, delta = (float(value) for value in values)
    c12, s12 = np.cos(t12), np.sin(t12)
    c13, s13 = np.cos(t13), np.sin(t13)
    c23, s23 = np.cos(t23), np.sin(t23)
    rotation12 = np.array(
        [[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]], dtype=np.complex128
    )
    rotation13 = np.array(
        [
            [c13, 0, s13 * np.exp(-1j * delta)],
            [0, 1, 0],
            [-s13 * np.exp(1j * delta), 0, c13],
        ],
        dtype=np.complex128,
    )
    rotation23 = np.array(
        [[1, 0, 0], [0, c23, s23], [0, -s23, c23]],
        dtype=np.complex128,
    )
    return rotation23 @ rotation13 @ rotation12
