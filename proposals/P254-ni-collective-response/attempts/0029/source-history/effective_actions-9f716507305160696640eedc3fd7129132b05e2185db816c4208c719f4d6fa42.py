"""Exact conditional elimination of quadratic heavy fields.

The functions here implement finite-dimensional action algebra.  They do not
identify a field with a vector meson, supply a kinetic action, choose a mass or
coupling, or turn an invariant local term into an anomaly.  Those are separate
model premises.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp
from sympy.matrices.exceptions import NonInvertibleMatrixError


@dataclass(frozen=True)
class QuadraticElimination:
    r"""Exact stationary reduction of ``L=V^T K V/2+V^T J``."""

    kernel: sp.ImmutableMatrix
    source: sp.ImmutableMatrix
    inverse_kernel: sp.ImmutableMatrix
    stationary_field: sp.ImmutableMatrix
    stationarity_residual: sp.ImmutableMatrix
    effective_term: sp.Expr


@dataclass(frozen=True)
class EvenOddElimination:
    """Parity decomposition of a quadratic-source effective term.

    The kernel is assumed parity even.  ``even_source`` is unchanged and
    ``odd_source`` changes sign under the declared parity bookkeeping.  The
    two square terms are even and the cross term is odd.
    """

    elimination: QuadraticElimination
    even_source: sp.ImmutableMatrix
    odd_source: sp.ImmutableMatrix
    even_square: sp.Expr
    odd_square: sp.Expr
    odd_cross: sp.Expr
    parity_transformed_effective_term: sp.Expr


@dataclass(frozen=True)
class InverseKernelExpansion:
    """Truncated Neumann inverse with its exact multiplication residuals."""

    mass_kernel: sp.ImmutableMatrix
    derivative_kernel: sp.ImmutableMatrix
    max_order: int
    approximation: sp.ImmutableMatrix
    left_residual: sp.ImmutableMatrix
    right_residual: sp.ImmutableMatrix


def _immutable_simplified(matrix: Any) -> sp.ImmutableMatrix:
    value = sp.Matrix(matrix).applyfunc(sp.simplify)
    return sp.ImmutableMatrix(value)


def _validated_kernel(kernel: Any) -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]:
    value = sp.Matrix(kernel)
    if value.rows == 0 or value.rows != value.cols:
        raise ValueError("kernel must be a nonempty square matrix")
    if _immutable_simplified(value - value.T) != sp.zeros(value.rows):
        raise ValueError("kernel must be symmetric")
    try:
        inverse = value.inv()
    except NonInvertibleMatrixError as exc:
        raise ValueError("kernel must be invertible") from exc
    return _immutable_simplified(value), _immutable_simplified(inverse)


def _validated_column(source: Any, dimension: int, name: str) -> sp.ImmutableMatrix:
    value = sp.Matrix(source)
    if value.shape != (dimension, 1):
        raise ValueError(f"{name} must be a {dimension}-entry column")
    return _immutable_simplified(value)


def quadratic_source_action(
    field: Any,
    kernel: Any,
    source: Any,
) -> sp.Expr:
    r"""Return ``V^T K V/2+V^T J`` in the declared plus-source convention."""

    matrix, _ = _validated_kernel(kernel)
    vector = _validated_column(field, matrix.rows, "field")
    current = _validated_column(source, matrix.rows, "source")
    return sp.simplify(
        (vector.T * matrix * vector)[0] / 2 + (vector.T * current)[0]
    )


def eliminate_quadratic_field(kernel: Any, source: Any) -> QuadraticElimination:
    r"""Eliminate ``V`` exactly from ``L=V^T K V/2+V^T J``.

    For a symmetric invertible kernel, stationarity gives ``V_*=-K^-1 J``
    and the reduced term is ``-J^T K^-1 J/2``.  All objects are derived from
    the supplied kernel and source; no propagator sign or target coefficient
    is inserted separately.
    """

    matrix, inverse = _validated_kernel(kernel)
    current = _validated_column(source, matrix.rows, "source")
    stationary = _immutable_simplified(-inverse * current)
    residual = _immutable_simplified(matrix * stationary + current)
    effective = sp.simplify(-(current.T * inverse * current)[0] / 2)
    return QuadraticElimination(
        kernel=matrix,
        source=current,
        inverse_kernel=inverse,
        stationary_field=stationary,
        stationarity_residual=residual,
        effective_term=effective,
    )


def eliminate_even_odd_sources(
    kernel: Any,
    even_source: Any,
    odd_source: Any,
) -> EvenOddElimination:
    """Derive even squares and the odd cross term after heavy-field elimination."""

    matrix, inverse = _validated_kernel(kernel)
    even = _validated_column(even_source, matrix.rows, "even_source")
    odd = _validated_column(odd_source, matrix.rows, "odd_source")
    elimination = eliminate_quadratic_field(matrix, even + odd)
    even_square = sp.simplify(-(even.T * inverse * even)[0] / 2)
    odd_square = sp.simplify(-(odd.T * inverse * odd)[0] / 2)
    odd_cross = sp.simplify(
        -((even.T * inverse * odd)[0] + (odd.T * inverse * even)[0]) / 2
    )
    parity_transformed = sp.simplify(even_square + odd_square - odd_cross)
    if sp.simplify(
        elimination.effective_term - (even_square + odd_square + odd_cross)
    ) != 0:
        raise AssertionError("source decomposition failed")
    return EvenOddElimination(
        elimination=elimination,
        even_source=even,
        odd_source=odd,
        even_square=even_square,
        odd_square=odd_square,
        odd_cross=odd_cross,
        parity_transformed_effective_term=parity_transformed,
    )


def low_momentum_inverse_expansion(
    mass_kernel: Any,
    derivative_kernel: Any,
    max_order: int,
) -> InverseKernelExpansion:
    r"""Return a finite Neumann expansion of ``(M+D)^-1`` and its residual.

    The approximation is

    ``sum_(n=0)^N (-M^-1 D)^n M^-1``.

    It is a formal low-momentum expansion only when the caller supplies a
    power counting in which ``M^-1 D`` is small.  The returned left and right
    residual matrices prevent a finite truncation from masquerading as the
    exact inverse.
    """

    if not isinstance(max_order, int):
        raise TypeError("max_order must be an integer")
    if max_order < 0:
        raise ValueError("max_order must be nonnegative")
    mass, inverse_mass = _validated_kernel(mass_kernel)
    derivative = sp.Matrix(derivative_kernel)
    if derivative.shape != mass.shape:
        raise ValueError("derivative_kernel must have the mass-kernel shape")
    if _immutable_simplified(derivative - derivative.T) != sp.zeros(mass.rows):
        raise ValueError("derivative_kernel must be symmetric")
    derivative = _immutable_simplified(derivative)
    ratio = sp.Matrix(inverse_mass * derivative)
    approximation = sp.zeros(mass.rows)
    for order in range(max_order + 1):
        approximation += (-1) ** order * ratio**order * inverse_mass
    approximation = _immutable_simplified(approximation)
    full_kernel = sp.Matrix(mass + derivative)
    identity = sp.eye(mass.rows)
    left_residual = _immutable_simplified(full_kernel * approximation - identity)
    right_residual = _immutable_simplified(approximation * full_kernel - identity)
    return InverseKernelExpansion(
        mass_kernel=mass,
        derivative_kernel=derivative,
        max_order=max_order,
        approximation=approximation,
        left_residual=left_residual,
        right_residual=right_residual,
    )


def stationary_reduced_variation(
    explicit_variation: Any,
    stationarity_residual: Any,
    induced_field_variation: Any,
) -> sp.Expr:
    r"""Return the chain-rule variation after a stationary field substitution.

    For ``S_eff[phi]=S[phi,V_*(phi)]`` this represents
    ``delta_explicit S + (delta S/delta V)|_* . delta V_*``.  An actual zero
    stationarity residual removes the induced-field term; this does not fix
    or generate the remaining explicit variation.
    """

    residual = sp.Matrix(stationarity_residual)
    induced = sp.Matrix(induced_field_variation)
    if residual.cols != 1 or induced.shape != residual.shape:
        raise ValueError("residual and induced variation must be equal-size columns")
    return sp.simplify(sp.sympify(explicit_variation) + (residual.T * induced)[0])
