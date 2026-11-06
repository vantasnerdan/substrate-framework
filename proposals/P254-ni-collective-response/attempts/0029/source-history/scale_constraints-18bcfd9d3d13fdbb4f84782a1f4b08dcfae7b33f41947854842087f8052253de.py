"""Exact diagnostics for positive monomial scale constraints.

The functions in this module are deliberately conditional.  They classify a
supplied exponent matrix, right-hand side, covariance, or interval collection;
they do not derive physical constraint rows, independence, or observations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp

from substrate_framework.linear_systems import (
    LinearSystemDiagnostics,
    diagnose_linear_system,
)


@dataclass(frozen=True)
class LogConstraintDiagnostics:
    """Exact structural ledger for ``design * log_scales = rhs``."""

    design: sp.ImmutableMatrix
    rhs: sp.ImmutableMatrix
    provenance: tuple[str, ...]
    linear: LinearSystemDiagnostics
    nullspace: tuple[sp.ImmutableMatrix, ...]
    left_nullspace: tuple[sp.ImmutableMatrix, ...]
    compatibility_residuals: tuple[sp.Expr, ...]
    coefficient_ranks_by_row: tuple[int, ...]
    augmented_ranks_by_row: tuple[int, ...]
    coefficient_informative_rows: tuple[bool, ...]
    augmented_informative_rows: tuple[bool, ...]
    coordinate_identifiable: tuple[bool, ...] | None


@dataclass(frozen=True)
class GeneralizedLeastSquaresLedger:
    """Exact GLS result for a declared positive-definite covariance."""

    design: sp.ImmutableMatrix
    rhs: sp.ImmutableMatrix
    provenance: tuple[str, ...]
    covariance: sp.ImmutableMatrix
    covariance_provenance: str
    precision: sp.ImmutableMatrix
    normal_matrix: sp.ImmutableMatrix
    estimator: sp.ImmutableMatrix
    fitted: sp.ImmutableMatrix
    residual: sp.ImmutableMatrix
    residual_projector: sp.ImmutableMatrix
    normal_residual: sp.ImmutableMatrix
    chi_squared: sp.Expr
    degrees_of_freedom: int


@dataclass(frozen=True)
class ClosedIntervalIntersection:
    """Intersection status of supplied closed intervals on one coordinate."""

    lower: sp.Expr
    upper: sp.Expr
    provenance: tuple[str, ...]
    feasible: bool
    point_identified: bool
    contradiction: bool
    active_lower_indices: tuple[int, ...]
    active_upper_indices: tuple[int, ...]


def _column(values: Any, *, rows: int, name: str) -> sp.Matrix:
    result = sp.Matrix(values)
    if result.cols != 1 or result.rows != rows:
        raise ValueError(f"{name} must be a column matching the design rows")
    return result


def _provenance(values: Sequence[str], *, rows: int) -> tuple[str, ...]:
    result = tuple(values)
    if len(result) != rows:
        raise ValueError("provenance must name every design row")
    if any(not isinstance(item, str) or not item.strip() for item in result):
        raise ValueError("provenance labels must be non-empty strings")
    return result


def _require_exact(name: str, *values: Any) -> None:
    for value in values:
        expressions = value if isinstance(value, sp.MatrixBase) else sp.Matrix(value)
        if any(sp.sympify(entry).has(sp.Float) for entry in expressions):
            raise ValueError(
                f"{name} must use exact real values rather than floating inputs"
            )
        if any(sp.sympify(entry).is_real is not True for entry in expressions):
            raise ValueError(f"{name} must be provably real")


def diagnose_log_constraints(
    design: Any,
    rhs: Any,
    *,
    provenance: Sequence[str],
) -> LogConstraintDiagnostics:
    """Diagnose an exact log-linear constraint system.

    A coordinate is identifiable precisely when every coefficient-null vector
    has zero in that coordinate.  This can hold in a rank-deficient system and
    does not require more equations than unknowns.  Identifiability is reported
    only for a consistent system.

    The left-null residuals are compatibility conditions: for every
    ``ell`` with ``ell.T * design == 0``, consistency requires
    ``ell.T * rhs == 0``.  Incremental ranks respect the supplied row order and
    are structural ledgers, not proofs that two physical data sources are
    statistically or causally independent.
    """

    matrix = sp.Matrix(design)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError("design must be non-empty")
    values = _column(rhs, rows=matrix.rows, name="rhs")
    _require_exact("design and rhs", matrix, values)
    labels = _provenance(provenance, rows=matrix.rows)
    linear = diagnose_linear_system(matrix, values)

    nullspace = tuple(sp.ImmutableMatrix(vector) for vector in matrix.nullspace())
    left_nullspace = tuple(
        sp.ImmutableMatrix(vector) for vector in matrix.T.nullspace()
    )
    compatibility = tuple(
        sp.simplify((vector.T * values)[0]) for vector in left_nullspace
    )

    coefficient_ranks: list[int] = []
    augmented_ranks: list[int] = []
    coefficient_informative: list[bool] = []
    augmented_informative: list[bool] = []
    old_coefficient_rank = 0
    old_augmented_rank = 0
    for stop in range(1, matrix.rows + 1):
        prefix = matrix[:stop, :]
        prefix_values = values[:stop, :]
        coefficient_rank = int(prefix.rank())
        augmented_rank = int(prefix.row_join(prefix_values).rank())
        coefficient_ranks.append(coefficient_rank)
        augmented_ranks.append(augmented_rank)
        coefficient_informative.append(coefficient_rank > old_coefficient_rank)
        augmented_informative.append(augmented_rank > old_augmented_rank)
        old_coefficient_rank = coefficient_rank
        old_augmented_rank = augmented_rank

    identifiable = None
    if linear.consistent:
        identifiable = tuple(
            all(sp.simplify(vector[index]) == 0 for vector in nullspace)
            for index in range(matrix.cols)
        )

    return LogConstraintDiagnostics(
        design=sp.ImmutableMatrix(matrix),
        rhs=sp.ImmutableMatrix(values),
        provenance=labels,
        linear=linear,
        nullspace=nullspace,
        left_nullspace=left_nullspace,
        compatibility_residuals=compatibility,
        coefficient_ranks_by_row=tuple(coefficient_ranks),
        augmented_ranks_by_row=tuple(augmented_ranks),
        coefficient_informative_rows=tuple(coefficient_informative),
        augmented_informative_rows=tuple(augmented_informative),
        coordinate_identifiable=identifiable,
    )


def positive_monomial_log_system(
    exponents: Any,
    observable_ratios: Any,
    coefficient_ratios: Any,
    *,
    provenance: Sequence[str],
) -> LogConstraintDiagnostics:
    r"""Convert positive dimensionless monomials into a log-linear system.

    The declared equations are

    ``observable_i = coefficient_i * product_j(scale_j**exponents_ij)``.

    Every observable, coefficient, and scale must already be expressed as a
    dimensionless positive ratio to a declared reference.  This function can
    check positivity assumptions visible to SymPy, but physical dimensions and
    the provenance of the supplied equations remain caller obligations.
    """

    matrix = sp.Matrix(exponents)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError("exponents must be non-empty")
    observables = _column(
        observable_ratios, rows=matrix.rows, name="observable_ratios"
    )
    coefficients = _column(
        coefficient_ratios, rows=matrix.rows, name="coefficient_ratios"
    )
    for name, values in (
        ("observable_ratios", observables),
        ("coefficient_ratios", coefficients),
    ):
        if any(sp.sympify(value).is_positive is not True for value in values):
            raise ValueError(f"{name} must be declared positive")
    rhs = sp.Matrix(
        [sp.log(sp.simplify(observables[index] / coefficients[index]))
         for index in range(matrix.rows)]
    )
    return diagnose_log_constraints(matrix, rhs, provenance=provenance)


def shift_log_references(
    system: LogConstraintDiagnostics,
    shifts: Any,
) -> LogConstraintDiagnostics:
    r"""Change log-coordinate references without changing physical residuals.

    If ``x_new = x_old - delta``, then the same equations are represented by
    ``design * x_new = rhs - design * delta``.
    """

    delta = _column(shifts, rows=system.design.cols, name="shifts")
    _require_exact("shifts", delta)
    shifted_rhs = sp.Matrix(system.rhs) - sp.Matrix(system.design) * delta
    return diagnose_log_constraints(
        system.design,
        shifted_rhs,
        provenance=system.provenance,
    )


def generalized_least_squares(
    design: Any,
    rhs: Any,
    covariance: Any,
    *,
    provenance: Sequence[str],
    covariance_provenance: str,
) -> GeneralizedLeastSquaresLedger:
    """Return the exact unique GLS ledger for supplied log-space data.

    The covariance must be symmetric and provably positive definite, and the
    design must have full column rank.  ``chi_squared`` is a supplied-model
    residual statistic only; this routine assigns no distribution, p-value, or
    physical independence to it.
    """

    matrix = sp.Matrix(design)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError("design must be non-empty")
    values = _column(rhs, rows=matrix.rows, name="rhs")
    labels = _provenance(provenance, rows=matrix.rows)
    if not isinstance(covariance_provenance, str) or not covariance_provenance.strip():
        raise ValueError("covariance_provenance must be a non-empty string")
    sigma = sp.Matrix(covariance)
    if sigma.shape != (matrix.rows, matrix.rows):
        raise ValueError("covariance must be square and match the design rows")
    _require_exact("design, rhs, and covariance", matrix, values, sigma)
    if sigma != sigma.T:
        raise ValueError("covariance must be symmetric")
    if sigma.is_positive_definite is not True:
        raise ValueError("covariance must be provably positive definite")
    if matrix.rank() != matrix.cols:
        raise ValueError("design must have full column rank for a unique GLS estimate")

    precision = sigma.inv()
    normal = matrix.T * precision * matrix
    estimator = normal.inv() * matrix.T * precision * values
    fitted = matrix * estimator
    residual = values - fitted
    projector = sp.eye(matrix.rows) - matrix * normal.inv() * matrix.T * precision
    normal_residual = sp.simplify(matrix.T * precision * residual)
    chi_squared = sp.simplify((residual.T * precision * residual)[0])

    return GeneralizedLeastSquaresLedger(
        design=sp.ImmutableMatrix(matrix),
        rhs=sp.ImmutableMatrix(values),
        provenance=labels,
        covariance=sp.ImmutableMatrix(sigma),
        covariance_provenance=covariance_provenance,
        precision=sp.ImmutableMatrix(precision),
        normal_matrix=sp.ImmutableMatrix(normal),
        estimator=sp.ImmutableMatrix(estimator),
        fitted=sp.ImmutableMatrix(fitted),
        residual=sp.ImmutableMatrix(residual),
        residual_projector=sp.ImmutableMatrix(projector),
        normal_residual=sp.ImmutableMatrix(normal_residual),
        chi_squared=chi_squared,
        degrees_of_freedom=matrix.rows - int(matrix.rank()),
    )


def _compare(left: sp.Expr, right: sp.Expr) -> int:
    """Return -1/0/1 for an exact, decidable real ordering."""

    difference = sp.simplify(left - right)
    if difference == 0:
        return 0
    if difference.is_negative is True:
        return -1
    if difference.is_positive is True:
        return 1
    raise ValueError("interval endpoint ordering must be exactly decidable")


def intersect_closed_intervals(
    intervals: Sequence[tuple[Any, Any]],
    *,
    provenance: Sequence[str],
) -> ClosedIntervalIntersection:
    """Intersect non-empty exact closed intervals for one log coordinate."""

    if not intervals:
        raise ValueError("at least one interval is required")
    parsed = tuple((sp.sympify(lower), sp.sympify(upper)) for lower, upper in intervals)
    labels = _provenance(provenance, rows=len(parsed))
    _require_exact("intervals", sp.Matrix(parsed))
    for lower, upper in parsed:
        if _compare(lower, upper) > 0:
            raise ValueError("each interval must satisfy lower <= upper")

    lower = parsed[0][0]
    upper = parsed[0][1]
    for candidate_lower, candidate_upper in parsed[1:]:
        if _compare(candidate_lower, lower) > 0:
            lower = candidate_lower
        if _compare(candidate_upper, upper) < 0:
            upper = candidate_upper

    ordering = _compare(lower, upper)
    feasible = ordering <= 0
    return ClosedIntervalIntersection(
        lower=lower,
        upper=upper,
        provenance=labels,
        feasible=feasible,
        point_identified=ordering == 0,
        contradiction=ordering > 0,
        active_lower_indices=tuple(
            index for index, (value, _) in enumerate(parsed) if _compare(value, lower) == 0
        ),
        active_upper_indices=tuple(
            index for index, (_, value) in enumerate(parsed) if _compare(value, upper) == 0
        ),
    )
