"""Exact carrier-factor and charge ledgers for a declared SU(2) doublet.

The standard two-dimensional SU(2) carrier, an independent projector carrier,
and a commuting Abelian charge are distinct mathematical objects.  Keeping
those factors separate prevents a rank-one projector on isospin space from
being relabelled as Lorentz chirality.  These helpers prove finite-dimensional
identities only; they do not identify physical states or construct gauge
dynamics, currents, anomalies, or interactions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


_PAULI = (
    sp.ImmutableMatrix([[0, 1], [1, 0]]),
    sp.ImmutableMatrix([[0, -sp.I], [sp.I, 0]]),
    sp.ImmutableMatrix([[1, 0], [0, -1]]),
)
_GENERATORS = tuple(
    sp.ImmutableMatrix((sp.Matrix(matrix) / 2).applyfunc(sp.simplify))
    for matrix in _PAULI
)
_CYCLIC = ((0, 1, 2), (1, 2, 0), (2, 0, 1))


def _immutable_simplified(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sp.Matrix(matrix).applyfunc(sp.simplify))


def _zero_matrix(matrix: sp.MatrixBase) -> bool:
    return _immutable_simplified(matrix) == sp.zeros(*matrix.shape)


def _exact_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.Matrix(value)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError(f"{name} must be nonempty")
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError(f"{name} must contain exact entries")
    return _immutable_simplified(matrix)


def _hermitian_projector(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = _exact_matrix(value, name)
    if matrix.rows != matrix.cols:
        raise ValueError(f"{name} must be square")
    if not _zero_matrix(matrix - matrix.H):
        raise ValueError(f"{name} must be Hermitian")
    if not _zero_matrix(matrix * matrix - matrix):
        raise ValueError(f"{name} must be idempotent")
    return matrix


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return sp.simplify(expression)


def _commutator_residuals(
    generators: Sequence[sp.MatrixBase],
) -> tuple[sp.ImmutableMatrix, ...]:
    return tuple(
        _immutable_simplified(
            generators[first] * generators[second]
            - generators[second] * generators[first]
            - sp.I * generators[target]
        )
        for first, second, target in _CYCLIC
    )


def _fundamental_commutant_basis() -> tuple[sp.ImmutableMatrix, ...]:
    elementary = (
        sp.ImmutableMatrix([[1, 0], [0, 0]]),
        sp.ImmutableMatrix([[0, 1], [0, 0]]),
        sp.ImmutableMatrix([[0, 0], [1, 0]]),
        sp.ImmutableMatrix([[0, 0], [0, 1]]),
    )
    rows: list[list[sp.Expr]] = []
    for generator in _GENERATORS:
        commutators = tuple(
            sp.Matrix(element * generator - generator * element)
            for element in elementary
        )
        for row in range(2):
            for column in range(2):
                rows.append(
                    [commutator[row, column] for commutator in commutators]
                )
    coefficient_matrix = sp.Matrix(rows)
    basis = []
    for vector in coefficient_matrix.nullspace():
        basis.append(
            _immutable_simplified(
                sp.Matrix(
                    [[vector[0], vector[1]], [vector[2], vector[3]]]
                )
            )
        )
    return tuple(basis)


@dataclass(frozen=True)
class SU2FundamentalLedger:
    """Exact standard fundamental-representation data."""

    generators: tuple[sp.ImmutableMatrix, ...]
    raising_operator: sp.ImmutableMatrix
    lowering_operator: sp.ImmutableMatrix
    casimir: sp.ImmutableMatrix
    commutator_residuals: tuple[sp.ImmutableMatrix, ...]
    commutant_basis: tuple[sp.ImmutableMatrix, ...]


@dataclass(frozen=True)
class SU2ChiralFactorLedger:
    """SU(2) action on an independent Hermitian-projector factor."""

    projector: sp.ImmutableMatrix
    complementary_projector: sp.ImmutableMatrix
    left_generators: tuple[sp.ImmutableMatrix, ...]
    right_generators: tuple[sp.ImmutableMatrix, ...]
    vector_generators: tuple[sp.ImmutableMatrix, ...]
    axial_generators: tuple[sp.ImmutableMatrix, ...]
    left_hermiticity_residuals: tuple[sp.ImmutableMatrix, ...]
    right_hermiticity_residuals: tuple[sp.ImmutableMatrix, ...]
    left_commutator_residuals: tuple[sp.ImmutableMatrix, ...]
    right_commutator_residuals: tuple[sp.ImmutableMatrix, ...]
    parity_left_to_right_residuals: tuple[sp.ImmutableMatrix, ...] | None
    parity_vector_even_residuals: tuple[sp.ImmutableMatrix, ...] | None
    parity_axial_odd_residuals: tuple[sp.ImmutableMatrix, ...] | None


@dataclass(frozen=True)
class SU2SameCarrierProjectorLedger:
    """Diagnostics for multiplying the irreducible generators by one projector."""

    projector: sp.ImmutableMatrix
    projector_rank: int
    projected_generators: tuple[sp.ImmutableMatrix, ...]
    hermiticity_residuals: tuple[sp.ImmutableMatrix, ...]
    commutator_residuals: tuple[sp.ImmutableMatrix, ...]
    projector_commutator_residuals: tuple[sp.ImmutableMatrix, ...]
    projector_in_fundamental_commutant: bool


@dataclass(frozen=True)
class SU2CommonChargeLedger:
    """Charge data for a common Abelian generator on one irreducible doublet."""

    abelian_eigenvalue: sp.Expr
    coefficient: sp.Expr
    abelian_generator: sp.ImmutableMatrix
    charge_operator: sp.ImmutableMatrix
    upper_eigenvalue: sp.Expr
    lower_eigenvalue: sp.Expr
    eigenvalue_separation: sp.Expr
    commutator_residuals: tuple[sp.ImmutableMatrix, ...]
    assigned_labels: tuple[sp.Expr, sp.Expr] | None
    assigned_label_residuals: tuple[sp.Expr, sp.Expr] | None
    labels_compatible: bool | None


def su2_fundamental_ledger() -> SU2FundamentalLedger:
    """Return the exact Pauli-half representation and its full commutant."""

    raising = _immutable_simplified(_GENERATORS[0] + sp.I * _GENERATORS[1])
    lowering = _immutable_simplified(_GENERATORS[0] - sp.I * _GENERATORS[1])
    casimir = _immutable_simplified(
        sum((generator * generator for generator in _GENERATORS), sp.zeros(2))
    )
    return SU2FundamentalLedger(
        generators=_GENERATORS,
        raising_operator=raising,
        lowering_operator=lowering,
        casimir=casimir,
        commutator_residuals=_commutator_residuals(_GENERATORS),
        commutant_basis=_fundamental_commutant_basis(),
    )


def su2_chiral_factor_ledger(
    projector: Any,
    *,
    parity_exchange: Any | None = None,
) -> SU2ChiralFactorLedger:
    """Return the exact independent-factor left/right SU(2) decomposition.

    ``projector`` acts on a carrier distinct from the standard isospin factor.
    If supplied, ``parity_exchange`` must be an exact unitary on that carrier
    which conjugates the projector to its complement.
    """

    left_projector = _hermitian_projector(projector, "projector")
    dimension = left_projector.rows
    identity = sp.ImmutableMatrix(sp.eye(dimension))
    right_projector = _immutable_simplified(identity - left_projector)
    left = tuple(
        _immutable_simplified(sp.kronecker_product(generator, left_projector))
        for generator in _GENERATORS
    )
    right = tuple(
        _immutable_simplified(sp.kronecker_product(generator, right_projector))
        for generator in _GENERATORS
    )
    vector = tuple(
        _immutable_simplified(left_generator + right_generator)
        for left_generator, right_generator in zip(left, right, strict=True)
    )
    axial = tuple(
        _immutable_simplified(left_generator - right_generator)
        for left_generator, right_generator in zip(left, right, strict=True)
    )

    parity_left_to_right = None
    parity_vector_even = None
    parity_axial_odd = None
    if parity_exchange is not None:
        exchange = _exact_matrix(parity_exchange, "parity_exchange")
        if exchange.shape != (dimension, dimension):
            raise ValueError("parity_exchange must match the projector")
        if not _zero_matrix(exchange.H * exchange - identity):
            raise ValueError("parity_exchange must be unitary")
        if not _zero_matrix(
            exchange * left_projector * exchange.H - right_projector
        ):
            raise ValueError("parity_exchange must exchange the projectors")
        full_exchange = _immutable_simplified(
            sp.kronecker_product(sp.eye(2), exchange)
        )
        parity_left_to_right = tuple(
            _immutable_simplified(
                full_exchange * left_generator * full_exchange.H
                - right_generator
            )
            for left_generator, right_generator in zip(left, right, strict=True)
        )
        parity_vector_even = tuple(
            _immutable_simplified(
                full_exchange * vector_generator * full_exchange.H
                - vector_generator
            )
            for vector_generator in vector
        )
        parity_axial_odd = tuple(
            _immutable_simplified(
                full_exchange * axial_generator * full_exchange.H
                + axial_generator
            )
            for axial_generator in axial
        )

    return SU2ChiralFactorLedger(
        projector=left_projector,
        complementary_projector=right_projector,
        left_generators=left,
        right_generators=right,
        vector_generators=vector,
        axial_generators=axial,
        left_hermiticity_residuals=tuple(
            _immutable_simplified(generator - generator.H) for generator in left
        ),
        right_hermiticity_residuals=tuple(
            _immutable_simplified(generator - generator.H) for generator in right
        ),
        left_commutator_residuals=_commutator_residuals(left),
        right_commutator_residuals=_commutator_residuals(right),
        parity_left_to_right_residuals=parity_left_to_right,
        parity_vector_even_residuals=parity_vector_even,
        parity_axial_odd_residuals=parity_axial_odd,
    )


def su2_same_carrier_projector_ledger(
    projector: Any,
) -> SU2SameCarrierProjectorLedger:
    """Diagnose the same-carrier product ``T_a*P`` exactly."""

    same_projector = _hermitian_projector(projector, "projector")
    if same_projector.shape != (2, 2):
        raise ValueError("a same-carrier projector must be two by two")
    projected = tuple(
        _immutable_simplified(generator * same_projector)
        for generator in _GENERATORS
    )
    projector_commutators = tuple(
        _immutable_simplified(
            same_projector * generator - generator * same_projector
        )
        for generator in _GENERATORS
    )
    return SU2SameCarrierProjectorLedger(
        projector=same_projector,
        projector_rank=int(same_projector.rank()),
        projected_generators=projected,
        hermiticity_residuals=tuple(
            _immutable_simplified(generator - generator.H)
            for generator in projected
        ),
        commutator_residuals=_commutator_residuals(projected),
        projector_commutator_residuals=projector_commutators,
        projector_in_fundamental_commutant=all(
            _zero_matrix(residual) for residual in projector_commutators
        ),
    )


def su2_common_charge_ledger(
    abelian_eigenvalue: Any,
    *,
    coefficient: Any = 1,
    assigned_labels: Sequence[Any] | None = None,
) -> SU2CommonChargeLedger:
    """Return the common-Abelian charge spectrum on the fundamental doublet."""

    eigenvalue = _exact_real(abelian_eigenvalue, "abelian_eigenvalue")
    charge_coefficient = _exact_real(coefficient, "coefficient")
    abelian = _immutable_simplified(eigenvalue * sp.eye(2))
    charge = _immutable_simplified(
        _GENERATORS[2] + charge_coefficient * abelian
    )
    upper = sp.simplify(charge_coefficient * eigenvalue + sp.Rational(1, 2))
    lower = sp.simplify(charge_coefficient * eigenvalue - sp.Rational(1, 2))

    label_tuple = None
    label_residuals = None
    labels_compatible = None
    if assigned_labels is not None:
        values = tuple(assigned_labels)
        if len(values) != 2:
            raise ValueError("assigned_labels must contain upper and lower labels")
        label_tuple = (
            _exact_real(values[0], "assigned_labels[0]"),
            _exact_real(values[1], "assigned_labels[1]"),
        )
        label_residuals = (
            sp.simplify(label_tuple[0] - upper),
            sp.simplify(label_tuple[1] - lower),
        )
        labels_compatible = all(residual == 0 for residual in label_residuals)

    return SU2CommonChargeLedger(
        abelian_eigenvalue=eigenvalue,
        coefficient=charge_coefficient,
        abelian_generator=abelian,
        charge_operator=charge,
        upper_eigenvalue=upper,
        lower_eigenvalue=lower,
        eigenvalue_separation=sp.simplify(upper - lower),
        commutator_residuals=tuple(
            _immutable_simplified(abelian * generator - generator * abelian)
            for generator in _GENERATORS
        ),
        assigned_labels=label_tuple,
        assigned_label_residuals=label_residuals,
        labels_compatible=labels_compatible,
    )
