"""Convention-explicit finite-dimensional matrix decompositions.

The functions in this module are matrix algebra only.  In particular, a
relative unitary basis matrix is not assigned a particle, flavor, or CKM
interpretation here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray


ComplexMatrix = NDArray[np.complex128]


def _finite_complex_matrix(matrix: Any) -> ComplexMatrix:
    value = np.asarray(matrix, dtype=np.complex128)
    if value.ndim != 2 or 0 in value.shape:
        raise ValueError("matrix must be a non-empty two-dimensional array")
    if not np.all(np.isfinite(value)):
        raise ValueError("matrix entries must be finite")
    return value


def unitarity_residual(matrix: Any) -> float:
    """Return the largest two-sided unitarity residual of a square matrix."""

    value = _finite_complex_matrix(matrix)
    if value.shape[0] != value.shape[1]:
        raise ValueError("unitarity requires a square matrix")
    identity = np.eye(value.shape[0], dtype=np.complex128)
    left = np.max(np.abs(value.conj().T @ value - identity))
    right = np.max(np.abs(value @ value.conj().T - identity))
    return float(max(left, right))


@dataclass(frozen=True)
class BiunitaryDecomposition:
    """A full SVD in the convention ``left.H @ M @ right = diagonal``.

    ``left`` and ``right`` are square unitary matrices whose columns are the
    left and right singular-vector bases.  ``diagonal`` has the same shape as
    the input matrix and contains the nonnegative singular values.
    """

    left: ComplexMatrix
    singular_values: NDArray[np.float64]
    right: ComplexMatrix
    diagonal: ComplexMatrix

    def reconstruct(self) -> ComplexMatrix:
        """Reconstruct the original matrix as ``left @ diagonal @ right.H``."""

        return self.left @ self.diagonal @ self.right.conj().T

    def diagonalize(self, matrix: Any) -> ComplexMatrix:
        """Apply the stored bases to a shape-compatible matrix."""

        value = _finite_complex_matrix(matrix)
        if value.shape != self.diagonal.shape:
            raise ValueError("matrix shape does not match this decomposition")
        return self.left.conj().T @ value @ self.right


def biunitary_decomposition(matrix: Any) -> BiunitaryDecomposition:
    """Compute a full finite-dimensional singular-value decomposition.

    The routine supports square or rectangular complex matrices, including
    rank-deficient and repeated-singular-value cases.  Individual singular
    vectors are not canonical inside degenerate or null subspaces.
    """

    value = _finite_complex_matrix(matrix)
    left, singular_values, right_adjoint = np.linalg.svd(value, full_matrices=True)
    right = right_adjoint.conj().T
    diagonal = np.zeros(value.shape, dtype=np.complex128)
    count = min(value.shape)
    diagonal[np.arange(count), np.arange(count)] = singular_values
    return BiunitaryDecomposition(
        left=np.asarray(left, dtype=np.complex128),
        singular_values=np.asarray(singular_values, dtype=np.float64),
        right=np.asarray(right, dtype=np.complex128),
        diagonal=diagonal,
    )


def gram_eigenvalues(matrix: Any) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return descending spectra of ``M M.H`` and ``M.H M``.

    The nonzero entries of both spectra are the squared singular values.  The
    longer spectrum contains the additional zeros required by matrix shape.
    """

    value = _finite_complex_matrix(matrix)
    left = np.linalg.eigvalsh(value @ value.conj().T)[::-1]
    right = np.linalg.eigvalsh(value.conj().T @ value)[::-1]
    return np.asarray(left, dtype=np.float64), np.asarray(right, dtype=np.float64)


def relative_left_basis(
    first: Any,
    second: Any,
    *,
    tolerance: float = 1e-12,
) -> ComplexMatrix:
    """Return ``first.H @ second`` for two same-size unitary column bases."""

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    first_value = _finite_complex_matrix(first)
    second_value = _finite_complex_matrix(second)
    if first_value.shape != second_value.shape:
        raise ValueError("basis matrices must have the same shape")
    if first_value.shape[0] != first_value.shape[1]:
        raise ValueError("basis matrices must be square")
    if unitarity_residual(first_value) > tolerance:
        raise ValueError("first basis is not unitary within tolerance")
    if unitarity_residual(second_value) > tolerance:
        raise ValueError("second basis is not unitary within tolerance")
    return first_value.conj().T @ second_value


@dataclass(frozen=True)
class RealSymmetricRotation:
    """A proper rotation diagonalizing a real symmetric two-by-two matrix."""

    angle: float
    rotation: NDArray[np.float64]
    diagonalized: NDArray[np.float64]


def real_symmetric_rotation(
    matrix: Any,
    *,
    tolerance: float = 1e-12,
) -> RealSymmetricRotation:
    """Diagonalize a real symmetric 2-by-2 matrix with a proper rotation.

    For ``[[a,b],[b,d]]`` the returned convention is
    ``R=[[cos(theta),sin(theta)],[-sin(theta),cos(theta)]]`` with
    ``theta=atan2(2*b,d-a)/2`` and ``R.T @ matrix @ R`` diagonal.  A scalar
    multiple of the identity uses the canonical representative ``theta=0``;
    its diagonalizing basis is otherwise nonunique.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    value = _finite_complex_matrix(matrix)
    if value.shape != (2, 2):
        raise ValueError("matrix must be 2 by 2")
    if np.max(np.abs(value.imag)) > tolerance:
        raise ValueError("matrix must be real within tolerance")
    real = value.real
    if np.max(np.abs(real - real.T)) > tolerance:
        raise ValueError("matrix must be symmetric within tolerance")
    a, b, d = float(real[0, 0]), float(real[0, 1]), float(real[1, 1])
    angle = 0.5 * float(np.arctan2(2.0 * b, d - a))
    cosine, sine = float(np.cos(angle)), float(np.sin(angle))
    rotation = np.array([[cosine, sine], [-sine, cosine]], dtype=np.float64)
    diagonalized = rotation.T @ real @ rotation
    return RealSymmetricRotation(angle, rotation, diagonalized)
