"""Status-gated three-factor gauge-only boundary-running evidence.

This module solves a declared inverse problem.  It does not select the gauge
coefficients, boundary ratios, low-scale constraints, or physical meaning of
the factors, and it does not add terms omitted by the supplied beta ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import sympy as sp
from scipy.optimize import least_squares

from .gauge_beta import GaugeCoefficientLedger
from .numerics import (
    IVPEvidence,
    NumericalFailure,
    SolverTolerances,
    solve_ivp_evidence,
)


@dataclass(frozen=True)
class ThreeFactorBoundaryProblem:
    """A supplied two-constraint inverse problem for three gauge factors."""

    coefficients: GaugeCoefficientLedger
    reference_scale: Any
    boundary_ratios: tuple[Any, Any, Any]
    low_constraint_matrix: tuple[
        tuple[Any, Any, Any],
        tuple[Any, Any, Any],
    ]
    low_constraint_targets: tuple[Any, Any]
    readout_weights: tuple[Any, Any, Any]
    readout_normalization: Any


@dataclass(frozen=True)
class ExactOneLoopBoundarySolution:
    """Exact affine solution obtained when the two-loop matrix is zero."""

    boundary_amplitude: sp.Expr
    scaled_log_span: sp.Expr
    log_scale_span: sp.Expr
    high_scale: sp.Expr
    low_inverse_couplings: tuple[sp.Expr, sp.Expr, sp.Expr]
    constraint_residuals: tuple[sp.Expr, sp.Expr]
    readout: sp.Expr


@dataclass(frozen=True)
class BoundaryRunningSolution:
    """Converged shooting and IVP diagnostics for one declared problem."""

    problem: ThreeFactorBoundaryProblem
    matrix_scale: float
    integration_method: str
    shooting_method: str
    tolerances: SolverTolerances
    root_tolerance: float
    residual_tolerance: float
    boundary_amplitude: float
    log_scale_span: float
    high_scale: float
    low_inverse_couplings: tuple[float, float, float]
    constraint_residuals: tuple[float, float]
    max_abs_constraint_residual: float
    readout: float
    minimum_inverse_coupling: float
    root_function_evaluations: int
    root_cost: float
    root_optimality: float
    integration: IVPEvidence


def _finite_float(value: Any, name: str) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be real and finite") from error
    if not np.isfinite(numeric):
        raise ValueError(f"{name} must be real and finite")
    return numeric


def _positive_float(value: Any, name: str) -> float:
    numeric = _finite_float(value, name)
    if numeric <= 0.0:
        raise ValueError(f"{name} must be positive")
    return numeric


def _problem_arrays(
    problem: ThreeFactorBoundaryProblem,
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    if not isinstance(problem, ThreeFactorBoundaryProblem):
        raise TypeError("problem must be a ThreeFactorBoundaryProblem")
    if not isinstance(problem.coefficients, GaugeCoefficientLedger):
        raise TypeError("coefficients must be a GaugeCoefficientLedger")
    if len(problem.coefficients.factors) != 3:
        raise ValueError("the boundary solver requires exactly three gauge factors")
    reference_scale = _positive_float(problem.reference_scale, "reference_scale")
    if len(problem.boundary_ratios) != 3:
        raise ValueError("boundary_ratios must contain three values")
    boundary = np.asarray(
        [
            _positive_float(value, f"boundary_ratios[{index}]")
            for index, value in enumerate(problem.boundary_ratios)
        ],
        dtype=np.float64,
    )
    constraints = np.asarray(problem.low_constraint_matrix, dtype=np.float64)
    if constraints.shape != (2, 3) or not np.all(np.isfinite(constraints)):
        raise ValueError("low_constraint_matrix must be a finite 2 by 3 matrix")
    if np.linalg.matrix_rank(constraints) != 2:
        raise ValueError("low_constraint_matrix must have row rank two")
    if len(problem.low_constraint_targets) != 2:
        raise ValueError("low_constraint_targets must contain two values")
    targets = np.asarray(
        [
            _finite_float(value, f"low_constraint_targets[{index}]")
            for index, value in enumerate(problem.low_constraint_targets)
        ],
        dtype=np.float64,
    )
    if np.any(targets == 0.0):
        raise ValueError("low_constraint_targets must be nonzero for scaled shooting")
    if len(problem.readout_weights) != 3:
        raise ValueError("readout_weights must contain three values")
    readout = np.asarray(
        [
            _finite_float(value, f"readout_weights[{index}]")
            for index, value in enumerate(problem.readout_weights)
        ],
        dtype=np.float64,
    )
    if not np.any(readout):
        raise ValueError("readout_weights must contain a nonzero value")
    normalization = _positive_float(
        problem.readout_normalization,
        "readout_normalization",
    )
    return reference_scale, boundary, constraints, targets, readout, normalization


def _coefficient_arrays(
    ledger: GaugeCoefficientLedger,
) -> tuple[np.ndarray, np.ndarray]:
    one_loop = np.asarray([float(value) for value in ledger.one_loop], dtype=np.float64)
    two_loop = np.asarray(
        [[float(value) for value in row] for row in ledger.two_loop_gauge_matrix],
        dtype=np.float64,
    )
    if one_loop.shape != (3,) or two_loop.shape != (3, 3):
        raise ValueError("coefficient ledger must contain a three-factor vector and matrix")
    if not np.all(np.isfinite(one_loop)) or not np.all(np.isfinite(two_loop)):
        raise ValueError("coefficient ledger must be finite when evaluated numerically")
    return one_loop, two_loop


def inverse_coupling_downward_rhs(
    coefficients: GaugeCoefficientLedger,
    inverse_couplings: Any,
    *,
    matrix_scale: Any = 1.0,
) -> np.ndarray:
    r"""Return ``d alpha_i^-1 / d log(Lambda/mu)`` gauge-only running.

    For the C-RGE-005 convention this is
    ``b_i/(2*pi) + k*sum_j B_ij/(8*pi^2*alpha_j^-1)``.
    """

    one_loop, two_loop = _coefficient_arrays(coefficients)
    state = np.asarray(inverse_couplings, dtype=np.float64)
    if state.shape != (3,) or not np.all(np.isfinite(state)):
        raise ValueError("inverse_couplings must be a finite three-vector")
    if np.any(state <= 0.0):
        raise NumericalFailure("inverse couplings left the positive domain")
    scale = _finite_float(matrix_scale, "matrix_scale")
    return one_loop / (2.0 * np.pi) + scale * (two_loop @ (1.0 / state)) / (
        8.0 * np.pi**2
    )


def direct_coupling_downward_rhs(
    coefficients: GaugeCoefficientLedger,
    couplings: Any,
    *,
    matrix_scale: Any = 1.0,
) -> np.ndarray:
    r"""Return direct-``g`` running toward lower scales for review oracles."""

    one_loop, two_loop = _coefficient_arrays(coefficients)
    state = np.asarray(couplings, dtype=np.float64)
    if state.shape != (3,) or not np.all(np.isfinite(state)):
        raise ValueError("couplings must be a finite three-vector")
    if np.any(state <= 0.0):
        raise NumericalFailure("couplings left the positive domain")
    scale = _finite_float(matrix_scale, "matrix_scale")
    loop = 16.0 * np.pi**2
    return -(
        one_loop * state**3 / loop
        + scale * state**3 * (two_loop @ state**2) / loop**2
    )


def exact_one_loop_boundary_solution(
    problem: ThreeFactorBoundaryProblem,
) -> ExactOneLoopBoundarySolution:
    """Solve the zero-matrix affine inverse problem with exact arithmetic."""

    _problem_arrays(problem)

    def exact(value: Any, name: str) -> sp.Expr:
        expression = sp.sympify(value)
        if expression.has(sp.Float):
            raise ValueError(f"{name} must be exact for the exact one-loop route")
        if expression.is_real is not True:
            raise ValueError(f"{name} must be provably real")
        return expression

    reference = exact(problem.reference_scale, "reference_scale")
    boundary = sp.Matrix(
        [exact(value, f"boundary_ratios[{index}]") for index, value in enumerate(problem.boundary_ratios)]
    )
    constraints = sp.Matrix(
        [
            [exact(value, f"low_constraint_matrix[{row}][{column}]") for column, value in enumerate(values)]
            for row, values in enumerate(problem.low_constraint_matrix)
        ]
    )
    targets = sp.Matrix(
        [
            exact(value, f"low_constraint_targets[{index}]")
            for index, value in enumerate(problem.low_constraint_targets)
        ]
    )
    beta = sp.Matrix(problem.coefficients.one_loop)
    design = sp.Matrix.hstack(constraints * boundary, constraints * beta)
    if design.det() == 0:
        raise ValueError("the exact one-loop boundary design must be nonsingular")
    amplitude, scaled_span = tuple(sp.simplify(value) for value in design.inv() * targets)
    if amplitude.is_positive is not True or scaled_span.is_positive is not True:
        raise ValueError("the exact solution requires positive amplitude and log span")
    low = sp.simplify(amplitude * boundary + scaled_span * beta)
    residuals = tuple(sp.simplify(value) for value in constraints * low - targets)
    weights = sp.Matrix(
        [exact(value, f"readout_weights[{index}]") for index, value in enumerate(problem.readout_weights)]
    )
    normalization = exact(problem.readout_normalization, "readout_normalization")
    if normalization.is_positive is not True:
        raise ValueError("readout_normalization must be positive")
    log_span = sp.simplify(2 * sp.pi * scaled_span)
    return ExactOneLoopBoundarySolution(
        boundary_amplitude=amplitude,
        scaled_log_span=scaled_span,
        log_scale_span=log_span,
        high_scale=sp.simplify(reference * sp.exp(log_span)),
        low_inverse_couplings=tuple(sp.simplify(value) for value in low),  # type: ignore[arg-type]
        constraint_residuals=residuals,  # type: ignore[arg-type]
        readout=sp.simplify((weights.dot(low)) / normalization),
    )


def _integrate_inverse_couplings(
    problem: ThreeFactorBoundaryProblem,
    boundary_amplitude: float,
    log_scale_span: float,
    *,
    matrix_scale: float,
    tolerances: SolverTolerances,
    method: str,
    sample_count: int,
) -> IVPEvidence:
    _reference, boundary, _constraints, _targets, _readout, _normalization = (
        _problem_arrays(problem)
    )
    amplitude = _positive_float(boundary_amplitude, "boundary_amplitude")
    span = _positive_float(log_scale_span, "log_scale_span")
    if sample_count < 2:
        raise ValueError("sample_count must be at least two")
    sample_times = np.linspace(0.0, span, int(sample_count))
    return solve_ivp_evidence(
        lambda _downward_log, state: inverse_coupling_downward_rhs(
            problem.coefficients,
            state,
            matrix_scale=matrix_scale,
        ),
        (0.0, span),
        amplitude * boundary,
        sample_times=sample_times,
        tolerances=tolerances,
        method=method,
    )


def solve_three_factor_boundary_running(
    problem: ThreeFactorBoundaryProblem,
    *,
    initial_boundary_amplitude: Any,
    initial_log_scale_span: Any,
    matrix_scale: Any = 1.0,
    tolerances: SolverTolerances = SolverTolerances(rtol=1.0e-10, atol=1.0e-12),
    integration_method: str = "DOP853",
    root_tolerance: float = 1.0e-11,
    residual_tolerance: float = 1.0e-8,
    sample_count: int = 257,
) -> BoundaryRunningSolution:
    """Solve for a positive boundary amplitude and logarithmic scale span.

    The shooting variables are logarithms, so every trial amplitude and span is
    positive.  Failed IVP trials receive a finite penalty; the returned solution
    exists only when least-squares and the final IVP both succeed and the two
    unscaled constraint residuals pass ``residual_tolerance``.
    """

    reference, _boundary, constraints, targets, readout, normalization = (
        _problem_arrays(problem)
    )
    scale = _finite_float(matrix_scale, "matrix_scale")
    initial_amplitude = _positive_float(
        initial_boundary_amplitude,
        "initial_boundary_amplitude",
    )
    initial_span = _positive_float(initial_log_scale_span, "initial_log_scale_span")
    root_tol = _positive_float(root_tolerance, "root_tolerance")
    residual_tol = _positive_float(residual_tolerance, "residual_tolerance")
    target_scales = np.maximum(np.abs(targets), 1.0)

    def scaled_residuals(log_parameters: np.ndarray) -> np.ndarray:
        amplitude, span = np.exp(log_parameters)
        if not np.all(np.isfinite((amplitude, span))):
            return np.full(2, 1.0e6, dtype=np.float64)
        try:
            evidence = _integrate_inverse_couplings(
                problem,
                float(amplitude),
                float(span),
                matrix_scale=scale,
                tolerances=tolerances,
                method=integration_method,
                sample_count=2,
            )
        except (NumericalFailure, ValueError, FloatingPointError):
            return np.full(2, 1.0e6, dtype=np.float64)
        low = evidence.state[:, -1]
        return (constraints @ low - targets) / target_scales

    root_result = least_squares(
        scaled_residuals,
        np.log([initial_amplitude, initial_span]),
        xtol=root_tol,
        ftol=root_tol,
        gtol=root_tol,
        max_nfev=300,
    )
    if not root_result.success or not np.all(np.isfinite(root_result.x)):
        raise NumericalFailure(f"shooting solve failed: {root_result.message}")
    amplitude, span = (float(value) for value in np.exp(root_result.x))
    integration = _integrate_inverse_couplings(
        problem,
        amplitude,
        span,
        matrix_scale=scale,
        tolerances=tolerances,
        method=integration_method,
        sample_count=sample_count,
    )
    low = integration.state[:, -1]
    residuals = constraints @ low - targets
    max_residual = float(np.max(np.abs(residuals)))
    if max_residual > residual_tol:
        raise NumericalFailure(
            f"shooting residual {max_residual:.3e} exceeds {residual_tol:.3e}"
        )
    minimum_inverse = float(np.min(integration.state))
    if minimum_inverse <= 0.0:
        raise NumericalFailure("inverse-coupling trajectory left the positive domain")
    high_scale = reference * float(np.exp(span))
    if not np.isfinite(high_scale):
        raise NumericalFailure("high scale overflowed")
    return BoundaryRunningSolution(
        problem=problem,
        matrix_scale=scale,
        integration_method=integration_method,
        shooting_method="scipy.optimize.least_squares(log amplitude, log span)",
        tolerances=tolerances,
        root_tolerance=root_tol,
        residual_tolerance=residual_tol,
        boundary_amplitude=amplitude,
        log_scale_span=span,
        high_scale=high_scale,
        low_inverse_couplings=tuple(float(value) for value in low),  # type: ignore[arg-type]
        constraint_residuals=tuple(float(value) for value in residuals),  # type: ignore[arg-type]
        max_abs_constraint_residual=max_residual,
        readout=float(readout @ low / normalization),
        minimum_inverse_coupling=minimum_inverse,
        root_function_evaluations=int(root_result.nfev),
        root_cost=float(root_result.cost),
        root_optimality=float(root_result.optimality),
        integration=integration,
    )
