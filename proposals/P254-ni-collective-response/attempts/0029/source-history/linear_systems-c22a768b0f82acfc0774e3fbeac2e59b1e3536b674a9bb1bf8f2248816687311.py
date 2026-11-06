"""Exact semantic diagnostics for finite linear systems."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


@dataclass(frozen=True)
class LinearSystemDiagnostics:
    """Rank-based logical status of ``coefficients * x = rhs``."""

    equations: int
    unknowns: int
    coefficient_rank: int
    augmented_rank: int
    consistent: bool
    solution_dimension: int | None
    unique: bool
    underdetermined: bool
    overdetermined_by_count: bool
    coefficient_row_dependencies: int


def diagnose_linear_system(
    coefficients: Any,
    rhs: Any,
) -> LinearSystemDiagnostics:
    """Classify a finite exact linear system by coefficient and augmented rank.

    ``overdetermined_by_count`` records only that there are more equations than
    unknowns. It does not imply consistency or uniqueness. An inconsistent
    system has no solution dimension.
    """

    matrix = sp.Matrix(coefficients)
    values = sp.Matrix(rhs)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError("coefficients must be non-empty")
    if values.cols != 1 or values.rows != matrix.rows:
        raise ValueError("rhs must be a column matching the coefficient rows")

    rank = int(matrix.rank())
    augmented_rank = int(matrix.row_join(values).rank())
    consistent = rank == augmented_rank
    solution_dimension = matrix.cols - rank if consistent else None
    return LinearSystemDiagnostics(
        equations=matrix.rows,
        unknowns=matrix.cols,
        coefficient_rank=rank,
        augmented_rank=augmented_rank,
        consistent=consistent,
        solution_dimension=solution_dimension,
        unique=consistent and rank == matrix.cols,
        underdetermined=consistent and rank < matrix.cols,
        overdetermined_by_count=matrix.rows > matrix.cols,
        coefficient_row_dependencies=matrix.rows - rank,
    )
