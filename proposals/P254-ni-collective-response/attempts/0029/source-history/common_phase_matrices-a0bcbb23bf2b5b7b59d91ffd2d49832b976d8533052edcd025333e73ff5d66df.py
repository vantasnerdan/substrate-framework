"""Algebra of finite real matrices multiplied by removable global phases.

The routines in this module are matrix algebra only.  They do not identify an
input matrix with a Yukawa coupling, a relative basis with a CKM matrix, or a
vanishing algebraic invariant with a physical CP theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any

import numpy as np
from numpy.typing import NDArray


ComplexMatrix = NDArray[np.complex128]
RealMatrix = NDArray[np.float64]


def _finite_real_matrix(matrix: Any) -> RealMatrix:
    value = np.asarray(matrix, dtype=np.complex128)
    if value.ndim != 2 or 0 in value.shape:
        raise ValueError("matrix must be non-empty and two-dimensional")
    if not np.all(np.isfinite(value)):
        raise ValueError("matrix entries must be finite")
    if np.any(value.imag != 0.0):
        raise ValueError("matrix must be exactly real")
    return np.asarray(value.real, dtype=np.float64)


def _finite_phase(phase: Real) -> float:
    if isinstance(phase, bool) or not isinstance(phase, Real):
        raise TypeError("phase must be a finite real scalar")
    value = float(phase)
    if not np.isfinite(value):
        raise ValueError("phase must be a finite real scalar")
    return value


@dataclass(frozen=True)
class CommonPhaseGrams:
    """One globally phased real matrix and its two Hermitian Gram matrices."""

    phased_matrix: ComplexMatrix
    left_gram: ComplexMatrix
    right_gram: ComplexMatrix


def common_phase_grams(matrix: Any, phase: Real = 0.0) -> CommonPhaseGrams:
    """Return both Grams of ``exp(i*phase)*R`` for an exactly real matrix ``R``.

    Mathematically the returned left and right Grams are respectively
    ``R*R.T`` and ``R.T*R`` and are independent of ``phase``.  The function
    constructs the phased matrix explicitly so current floating-point
    consumers also exercise the cancellation rather than receiving a copied
    expected result.
    """

    real = _finite_real_matrix(matrix)
    angle = _finite_phase(phase)
    phased = np.asarray(np.exp(1j * angle) * real, dtype=np.complex128)
    return CommonPhaseGrams(
        phased_matrix=phased,
        left_gram=np.asarray(phased @ phased.conj().T, dtype=np.complex128),
        right_gram=np.asarray(phased.conj().T @ phased, dtype=np.complex128),
    )


@dataclass(frozen=True)
class RealGramRelativeBasis:
    """Real spectral bases and their relative orthogonal transformation."""

    first_gram: RealMatrix
    second_gram: RealMatrix
    first_eigenvalues: NDArray[np.float64]
    second_eigenvalues: NDArray[np.float64]
    first_basis: RealMatrix
    second_basis: RealMatrix
    relative_basis: RealMatrix
    commutator: RealMatrix


def real_gram_relative_basis(
    first: Any,
    second: Any,
    *,
    first_phase: Real = 0.0,
    second_phase: Real = 0.0,
    tolerance: float = 1e-12,
) -> RealGramRelativeBasis:
    """Diagonalize two common-phase left Grams in real orthogonal bases.

    The input matrices may be rectangular but must have the same row count.
    Eigenvectors inside repeated eigenspaces are noncanonical.  This routine
    chooses a real representative; it does not assert that an arbitrary
    complex basis chosen inside a degenerate subspace has real entries.
    """

    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    first_real = _finite_real_matrix(first)
    second_real = _finite_real_matrix(second)
    if first_real.shape[0] != second_real.shape[0]:
        raise ValueError("matrices must have the same row dimension")

    first_result = common_phase_grams(first_real, first_phase)
    second_result = common_phase_grams(second_real, second_phase)
    imaginary_residual = max(
        float(np.max(np.abs(first_result.left_gram.imag))),
        float(np.max(np.abs(second_result.left_gram.imag))),
    )
    if imaginary_residual > tolerance:
        raise ArithmeticError("common-phase Gram cancellation exceeded tolerance")

    first_gram = np.asarray(first_result.left_gram.real, dtype=np.float64)
    second_gram = np.asarray(second_result.left_gram.real, dtype=np.float64)
    first_values, first_basis = np.linalg.eigh(first_gram)
    second_values, second_basis = np.linalg.eigh(second_gram)
    relative = np.asarray(first_basis.T @ second_basis, dtype=np.float64)
    commutator = np.asarray(
        first_gram @ second_gram - second_gram @ first_gram,
        dtype=np.float64,
    )
    return RealGramRelativeBasis(
        first_gram=first_gram,
        second_gram=second_gram,
        first_eigenvalues=np.asarray(first_values, dtype=np.float64),
        second_eigenvalues=np.asarray(second_values, dtype=np.float64),
        first_basis=np.asarray(first_basis, dtype=np.float64),
        second_basis=np.asarray(second_basis, dtype=np.float64),
        relative_basis=relative,
        commutator=commutator,
    )


def odd_antisymmetric_trace(
    matrix: Any,
    power: int,
    *,
    tolerance: float = 1e-12,
) -> float:
    """Return the trace of an odd power of a real antisymmetric matrix.

    The exact value is zero.  The numeric return exposes floating-point
    regression residuals while validation prevents applying the identity to a
    matrix outside its hypotheses.
    """

    if isinstance(power, bool) or not isinstance(power, Integral):
        raise TypeError("power must be an odd positive integer")
    power = int(power)
    if power < 1 or power % 2 == 0:
        raise ValueError("power must be an odd positive integer")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    value = _finite_real_matrix(matrix)
    if value.shape[0] != value.shape[1]:
        raise ValueError("matrix must be square")
    if float(np.max(np.abs(value + value.T))) > tolerance:
        raise ValueError("matrix must be antisymmetric within tolerance")
    return float(np.trace(np.linalg.matrix_power(value, power)))
