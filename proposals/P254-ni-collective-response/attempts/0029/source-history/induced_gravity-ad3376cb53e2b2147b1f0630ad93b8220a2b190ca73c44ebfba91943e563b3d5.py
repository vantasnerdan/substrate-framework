"""Exact conditional ledgers for cutoff-induced Newton scaling.

The functions in this module perform dimension, additive inverse-coupling,
and log-identifiability bookkeeping. They do not derive a quantum field
spectrum, regulator, heat-kernel coefficient, cutoff ontology, vanishing bare
term, Einstein dynamics, optical-source normalization, or observed scale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp

from .dimensional_analysis import monomial_exponents
from .scale_constraints import LogConstraintDiagnostics, diagnose_log_constraints


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _nonzero_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and (
        expression.is_real is not True or expression.is_zero is not False
    ):
        raise ValueError(f"{name} must be nonzero and real")
    return expression


def _dimension_column(value: Sequence[Any], name: str) -> sp.ImmutableMatrix:
    column = sp.Matrix(value)
    if column.shape != (3, 1):
        raise ValueError(f"{name} must be an M,L,T dimension column")
    return sp.ImmutableMatrix(column)


@dataclass(frozen=True)
class NewtonDimensionLedger:
    """Unique ``(a,c,hbar)`` monomial powers for ``G`` and ``1/G``."""

    dimension_matrix: sp.ImmutableMatrix
    newton_dimension: sp.ImmutableMatrix
    inverse_newton_dimension: sp.ImmutableMatrix
    newton_exponents: sp.ImmutableMatrix
    inverse_newton_exponents: sp.ImmutableMatrix


@dataclass(frozen=True)
class InducedInverseNewtonLedger:
    """A declared cutoff contribution and an independent additive baseline."""

    cutoff_length: sp.Expr
    signal_speed: sp.Expr
    action_scale: sp.Expr
    induced_coefficient: sp.Expr
    cutoff_energy: sp.Expr
    induced_inverse_newton: sp.Expr
    baseline_inverse_newton: sp.Expr
    total_inverse_newton: sp.Expr
    pure_induced_newton: sp.Expr


@dataclass(frozen=True)
class GravitySourceNormalizationLedger:
    """Dimension requirements for mapping Newton ``G`` to a source coupling."""

    operator_dimension: sp.ImmutableMatrix
    source_dimension: sp.ImmutableMatrix
    required_coupling_dimension: sp.ImmutableMatrix
    newton_dimension: sp.ImmutableMatrix
    normalization_dimension: sp.ImmutableMatrix
    dimensionless_normalization_allowed: bool


def newton_dimension_ledger() -> NewtonDimensionLedger:
    """Return exact M,L,T dimensions for columns ``(a,c,hbar)``.

    The target ``[G]=M^-1 L^3 T^-2`` has powers ``(2,3,-1)``; its reciprocal
    has ``(-2,-3,1)``. A dimensionless multiplier remains outside this solve.
    """

    matrix = sp.Matrix(
        [
            [0, 0, 1],
            [1, 1, 2],
            [0, -1, -1],
        ]
    )
    newton = sp.Matrix([-1, 3, -2])
    inverse = -newton
    return NewtonDimensionLedger(
        dimension_matrix=sp.ImmutableMatrix(matrix),
        newton_dimension=sp.ImmutableMatrix(newton),
        inverse_newton_dimension=sp.ImmutableMatrix(inverse),
        newton_exponents=sp.ImmutableMatrix(monomial_exponents(matrix, newton)),
        inverse_newton_exponents=sp.ImmutableMatrix(
            monomial_exponents(matrix, inverse)
        ),
    )


def induced_inverse_newton_ledger(
    cutoff_length: Any,
    induced_coefficient: Any,
    signal_speed: Any,
    action_scale: Any,
    *,
    baseline_inverse_newton: Any = 0,
) -> InducedInverseNewtonLedger:
    """Return a declared leading cutoff shift plus an additive baseline.

    The declared leading contribution is
    ``Delta(1/G)=s*hbar/(a**2*c**3)``. The algebraically equivalent energy
    form uses ``E_cut=hbar*c/a`` and
    ``Delta(1/G)=s*E_cut**2/(hbar*c**5)``. The baseline is an independent
    bare, renormalized, logarithmic, or finite contribution.
    """

    length = _positive(cutoff_length, "cutoff_length")
    speed = _positive(signal_speed, "signal_speed")
    action = _positive(action_scale, "action_scale")
    coefficient = _nonzero_real(induced_coefficient, "induced_coefficient")
    baseline = sp.sympify(baseline_inverse_newton)
    cutoff_energy = sp.simplify(action * speed / length)
    shift = sp.simplify(coefficient * action / (length**2 * speed**3))
    total = sp.simplify(baseline + shift)
    return InducedInverseNewtonLedger(
        cutoff_length=length,
        signal_speed=speed,
        action_scale=action,
        induced_coefficient=coefficient,
        cutoff_energy=cutoff_energy,
        induced_inverse_newton=shift,
        baseline_inverse_newton=baseline,
        total_inverse_newton=total,
        pure_induced_newton=sp.simplify(1 / shift),
    )


def effective_newton_from_inverse(total_inverse_newton: Any) -> sp.Expr:
    """Invert a separately established positive total inverse coupling."""

    inverse = _positive(total_inverse_newton, "total_inverse_newton")
    return sp.simplify(1 / inverse)


def cutoff_length_from_pure_induced_newton(
    newton_constant: Any,
    induced_coefficient: Any,
    signal_speed: Any,
    action_scale: Any,
) -> sp.Expr:
    """Infer ``a`` after supplying ``G,s,c,hbar`` and a zero baseline.

    This inverse is conditional data reduction, not a prediction of the
    cutoff length or coefficient.
    """

    newton = _positive(newton_constant, "newton_constant")
    coefficient = _positive(induced_coefficient, "induced_coefficient")
    speed = _positive(signal_speed, "signal_speed")
    action = _positive(action_scale, "action_scale")
    return sp.simplify(sp.sqrt(coefficient * action * newton / speed**3))


def induced_scaling_log_constraint(
    reduced_newton_ratio: Any,
    *,
    provenance: str,
) -> LogConstraintDiagnostics:
    """Diagnose ``2*log(a_ratio)-log(s_ratio)=log(G_reduced_ratio)``.

    All ratios are dimensionless relative to declared compatible references.
    The row ``(2,-1)`` has nullspace generated by ``(1,2)``; one supplied
    Newton ratio identifies neither input ratio when both remain free.
    """

    ratio = _positive(reduced_newton_ratio, "reduced_newton_ratio")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("provenance must be a non-empty string")
    return diagnose_log_constraints(
        [[2, -1]],
        [sp.log(ratio)],
        provenance=[provenance],
    )


def induced_scaling_null_rescaling(
    cutoff_ratio: Any,
    coefficient_ratio: Any,
    factor: Any,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Apply ``a_ratio->rho*a_ratio`` and ``s_ratio->rho**2*s_ratio``."""

    length = _positive(cutoff_ratio, "cutoff_ratio")
    coefficient = _positive(coefficient_ratio, "coefficient_ratio")
    rho = _positive(factor, "factor")
    changed_length = sp.simplify(rho * length)
    changed_coefficient = sp.simplify(rho**2 * coefficient)
    invariant = sp.simplify(length**2 / coefficient)
    return changed_length, changed_coefficient, invariant


def gravity_source_normalization_ledger(
    operator_dimension: Sequence[Any],
    source_dimension: Sequence[Any],
) -> GravitySourceNormalizationLedger:
    """Return dimensions required by ``operator = kappa*source``.

    Dimension columns use M,L,T order. A direct dimensionless multiplier
    ``kappa=alpha*G`` is allowed only when the required source-coupling
    dimension equals ``[G]``. Otherwise ``alpha`` must carry the returned
    normalization dimension; a numeral such as ``8*pi`` cannot supply it.
    """

    operator = _dimension_column(operator_dimension, "operator_dimension")
    source = _dimension_column(source_dimension, "source_dimension")
    required = sp.ImmutableMatrix(operator - source)
    newton = newton_dimension_ledger().newton_dimension
    normalization = sp.ImmutableMatrix(required - newton)
    return GravitySourceNormalizationLedger(
        operator_dimension=operator,
        source_dimension=source,
        required_coupling_dimension=required,
        newton_dimension=newton,
        normalization_dimension=normalization,
        dimensionless_normalization_allowed=normalization == sp.zeros(3, 1),
    )


def normalized_gravity_source_coupling(
    newton_constant: Any,
    normalization_factor: Any,
) -> sp.Expr:
    """Return ``alpha*G`` for a separately supplied dimensioned ``alpha``."""

    newton = _positive(newton_constant, "newton_constant")
    normalization = sp.sympify(normalization_factor)
    if normalization.is_number and normalization.is_real is not True:
        raise ValueError("normalization_factor must be real")
    return sp.simplify(normalization * newton)
