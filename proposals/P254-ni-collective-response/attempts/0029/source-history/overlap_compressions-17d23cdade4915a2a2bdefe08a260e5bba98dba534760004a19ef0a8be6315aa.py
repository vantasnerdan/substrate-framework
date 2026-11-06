"""Exact finite compressions of real multiplication operators.

The objects here are Hilbert-space and matrix algebra.  They do not identify
the supplied modes with generations, the multiplier with a Yukawa field, or a
relative eigenbasis with a physical mixing matrix.  Such an interpretation
requires separately accepted fields, interactions, sector maps, spectra, and
current operators.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and result.is_real is not True:
        raise ValueError(f"{name} must be real")
    return result


def _positive_real(value: Any, name: str) -> sp.Expr:
    result = _real(value, name)
    if result.is_number and result.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return result


def _square_matrix(matrix: Any, name: str = "matrix") -> sp.Matrix:
    value = sp.Matrix(matrix)
    if value.rows == 0 or value.rows != value.cols:
        raise ValueError(f"{name} must be non-empty and square")
    return value


def _is_zero_matrix(matrix: sp.MatrixBase) -> bool:
    return all(sp.simplify(entry) == 0 for entry in matrix)


def _exact_integral(
    integrand: sp.Expr,
    coordinate: sp.Symbol,
    lower: sp.Expr,
    upper: sp.Expr,
) -> sp.Expr:
    result = sp.integrate(integrand, (coordinate, lower, upper))
    if result.has(sp.Integral):
        raise ValueError("the declared overlap did not integrate exactly")
    return sp.simplify(result)


@dataclass(frozen=True)
class MultiplicationCompression:
    """An exact compression in one declared orthonormal ordered basis."""

    coordinate: sp.Symbol
    lower: sp.Expr
    upper: sp.Expr
    modes: tuple[sp.Expr, ...]
    multiplier: sp.Expr
    gram_matrix: sp.ImmutableMatrix
    matrix: sp.ImmutableMatrix


@dataclass(frozen=True)
class QuarticAsymmetricCompressionLedger:
    """C-QBL-003 mode compression for an explicitly asymmetric profile."""

    amplitude: sp.Expr
    asymmetry: sp.Expr
    inverse_width: sp.Expr
    even_mode_norm: sp.Expr
    odd_mode_norm: sp.Expr
    even_profile_matrix: sp.ImmutableMatrix
    odd_profile_matrix: sp.ImmutableMatrix
    matrix: sp.ImmutableMatrix


@dataclass(frozen=True)
class SpectralMultiplicityLedger:
    """Exact eigenvalue multiplicities and additional degenerate freedom."""

    dimension: int
    eigenvalue_multiplicities: tuple[tuple[sp.Expr, int], ...]
    degenerate_subspace_dimensions: tuple[int, ...]
    additional_unitary_parameters: int


def multiplication_compression(
    modes: Sequence[Any],
    multiplier: Any,
    coordinate: sp.Symbol,
    *,
    lower: Any = -sp.oo,
    upper: Any = sp.oo,
) -> MultiplicationCompression:
    """Return ``Y_ij = integral conjugate(eta_i) Phi eta_j dx`` exactly.

    The ordered modes must be orthonormal on the declared interval and the
    multiplier must be provably real.  Complex conjugation is load-bearing:
    omitting it changes both the Gram matrix and compression for complex
    modes.  If the profile has essential range ``lower <= Phi <= upper``, then
    ``c.H*Y*c`` is the normalized profile expectation of the reconstructed
    mode and every eigenvalue of ``Y`` lies in that interval.  Unresolved
    symbolic integrals are rejected instead of being presented as verified
    closed forms.
    """

    if not isinstance(coordinate, sp.Symbol) or coordinate.is_real is not True:
        raise ValueError("coordinate must be a real SymPy symbol")
    basis = tuple(sp.sympify(mode) for mode in modes)
    if not basis:
        raise ValueError("modes must be non-empty")
    profile = sp.sympify(multiplier)
    if sp.simplify(sp.conjugate(profile) - profile) != 0:
        raise ValueError("multiplier must be provably real")
    lower_bound = sp.sympify(lower)
    upper_bound = sp.sympify(upper)

    gram = sp.Matrix(
        [
            [
                _exact_integral(
                    sp.conjugate(left) * right,
                    coordinate,
                    lower_bound,
                    upper_bound,
                )
                for right in basis
            ]
            for left in basis
        ]
    )
    if not _is_zero_matrix(gram - sp.eye(len(basis))):
        raise ValueError("modes must be orthonormal on the declared interval")

    matrix = sp.Matrix(
        [
            [
                _exact_integral(
                    sp.conjugate(left) * profile * right,
                    coordinate,
                    lower_bound,
                    upper_bound,
                )
                for right in basis
            ]
            for left in basis
        ]
    )
    if not _is_zero_matrix(matrix - matrix.adjoint()):
        raise ValueError("real multiplication compression must be Hermitian")
    return MultiplicationCompression(
        coordinate=coordinate,
        lower=lower_bound,
        upper=upper_bound,
        modes=basis,
        multiplier=profile,
        gram_matrix=sp.ImmutableMatrix(gram),
        matrix=sp.ImmutableMatrix(matrix),
    )


def unitary_similarity(matrix: Any, basis_change: Any) -> sp.ImmutableMatrix:
    """Return ``U.H M U`` after exact square-unitary validation."""

    value = _square_matrix(matrix)
    unitary = _square_matrix(basis_change, "basis_change")
    if value.shape != unitary.shape:
        raise ValueError("matrix and basis_change must have the same shape")
    identity = sp.eye(unitary.rows)
    if not _is_zero_matrix(unitary.adjoint() * unitary - identity) or not _is_zero_matrix(
        unitary * unitary.adjoint() - identity
    ):
        raise ValueError("basis_change must be exactly unitary")
    return sp.ImmutableMatrix(
        value.rows,
        value.cols,
        lambda row, column: sp.simplify(
            (unitary.adjoint() * value * unitary)[row, column]
        ),
    )


def parity_forces_zero(
    left_mode_parity: int,
    multiplier_parity: int,
    right_mode_parity: int,
) -> bool:
    """Return whether the whole-line integrand has odd total parity.

    Each parity must be ``+1`` or ``-1``.  In particular, an even multiplier
    has zero cross overlap between one even and one odd mode regardless of its
    even width parameter.
    """

    parities = (left_mode_parity, multiplier_parity, right_mode_parity)
    if any(isinstance(value, bool) or value not in (-1, 1) for value in parities):
        raise ValueError("parities must each be +1 or -1")
    return left_mode_parity * multiplier_parity * right_mode_parity == -1


def matrix_commutator(first: Any, second: Any) -> sp.ImmutableMatrix:
    """Return the exact commutator ``first*second-second*first``."""

    left = _square_matrix(first, "first")
    right = _square_matrix(second, "second")
    if left.shape != right.shape:
        raise ValueError("matrices must have the same shape")
    product = left * right - right * left
    return sp.ImmutableMatrix(
        product.rows,
        product.cols,
        lambda row, column: sp.simplify(product[row, column]),
    )


def commuting_hermitian(first: Any, second: Any) -> bool:
    """Test the exact simultaneous-unitary-diagonalization criterion.

    Two finite Hermitian matrices admit a common unitary eigenbasis exactly
    when this function returns true.  Matrix inequality alone is neither
    necessary nor sufficient for eigenbasis misalignment.
    """

    left = _square_matrix(first, "first")
    right = _square_matrix(second, "second")
    if left.shape != right.shape:
        raise ValueError("matrices must have the same shape")
    if not _is_zero_matrix(left - left.adjoint()):
        raise ValueError("first must be Hermitian")
    if not _is_zero_matrix(right - right.adjoint()):
        raise ValueError("second must be Hermitian")
    return _is_zero_matrix(left * right - right * left)


def real_symmetric_commutator_scalar(first: Any, second: Any) -> sp.Expr:
    """Return the load-bearing off-diagonal commutator scalar for 2 by 2.

    For ``[[a,b],[b,d]]`` and ``[[e,f],[f,h]]`` this is
    ``f*(a-d) + b*(h-e)``.  The matrices commute exactly when it vanishes.
    """

    left = _square_matrix(first, "first")
    right = _square_matrix(second, "second")
    if left.shape != (2, 2) or right.shape != (2, 2):
        raise ValueError("matrices must both be 2 by 2")
    for value, name in ((left, "first"), (right, "second")):
        if not _is_zero_matrix(value - value.T):
            raise ValueError(f"{name} must be symmetric")
        if any(sp.simplify(sp.conjugate(entry) - entry) != 0 for entry in value):
            raise ValueError(f"{name} must be provably real")
    return sp.simplify((left * right - right * left)[0, 1])


def spectral_multiplicity_ledger(matrix: Any) -> SpectralMultiplicityLedger:
    """Return exact multiplicities and extra ``U(m)`` degeneracy freedoms.

    The reported additional parameter count is ``sum(m**2-m)`` over repeated
    eigenspaces.  It quotients the ``m`` independent eigenvector phases already
    present without degeneracy and therefore isolates the extra freedom caused
    by each degenerate block.
    """

    value = _square_matrix(matrix)
    if not _is_zero_matrix(value - value.adjoint()):
        raise ValueError("matrix must be Hermitian")
    multiplicities = tuple(
        (sp.simplify(eigenvalue), int(multiplicity))
        for eigenvalue, multiplicity in value.eigenvals().items()
    )
    repeated = tuple(multiplicity for _, multiplicity in multiplicities if multiplicity > 1)
    return SpectralMultiplicityLedger(
        dimension=value.rows,
        eigenvalue_multiplicities=multiplicities,
        degenerate_subspace_dimensions=repeated,
        additional_unitary_parameters=sum(size**2 - size for size in repeated),
    )


def quartic_asymmetric_compression_ledger(
    amplitude: Any,
    asymmetry: Any,
    inverse_width: Any,
) -> QuarticAsymmetricCompressionLedger:
    """Compress ``A*sech(z)*(1+b*tanh(z))`` onto C-QBL-003 modes.

    The unnormalized ordered modes are the accepted even ``sech(z)**2`` and
    odd ``sech(z)*tanh(z)`` shapes, with ``z=kappa*x``.  Exact normalization
    gives

    ``[[9*pi*A/32, sqrt(2)*A*b/5],
       [sqrt(2)*A*b/5, 3*pi*A/16]]``.

    The common inverse width cancels.  The off-diagonal entry is nonzero only
    because the separately supplied ``b*tanh(z)`` term is odd; it is not
    generated by changing the width of a centered even multiplier.
    """

    coefficient = _real(amplitude, "amplitude")
    odd_weight = _real(asymmetry, "asymmetry")
    kappa = _positive_real(inverse_width, "inverse_width")
    even_norm = sp.simplify(4 / (3 * kappa))
    odd_norm = sp.simplify(2 / (3 * kappa))
    even_profile = sp.ImmutableMatrix(
        [
            [9 * sp.pi * coefficient / 32, 0],
            [0, 3 * sp.pi * coefficient / 16],
        ]
    )
    cross = sp.simplify(sp.sqrt(2) * coefficient * odd_weight / 5)
    odd_profile = sp.ImmutableMatrix([[0, cross], [cross, 0]])
    return QuarticAsymmetricCompressionLedger(
        amplitude=coefficient,
        asymmetry=odd_weight,
        inverse_width=kappa,
        even_mode_norm=even_norm,
        odd_mode_norm=odd_norm,
        even_profile_matrix=even_profile,
        odd_profile_matrix=odd_profile,
        matrix=sp.ImmutableMatrix(even_profile + odd_profile),
    )
