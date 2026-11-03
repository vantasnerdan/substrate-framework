"""Finite-wall numerical evidence for the averaged spherical scalar model.

This module solves the exact reduced equations from
``spherical_einstein_scalar``.  Its Robin condition is an asymptotic
finite-radius approximation.  A converged object is numerical evidence for
that boundary-value problem, not a proof of an infinite-domain solution or a
pointwise solution of the unaveraged Einstein--scalar PDE.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import CubicSpline
from scipy.optimize import root
from scipy.special import expit, jv

from .numerics import (
    BVPEvidence,
    IVPEvidence,
    NumericalFailure,
    SolverTolerances,
    solve_bvp_evidence,
    solve_ivp_evidence,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class StaticSphericalScalarBVPSolution:
    """Converged adaptive-collocation evidence for one finite-wall branch."""

    radius: FloatArray
    state: FloatArray
    frequency: float
    fitted_parameter: float
    central_amplitude: float
    dimensionless_coupling: float
    origin_epsilon: float
    outer_radius: float
    initial_mesh_points: int
    tolerance: float
    collocation_rms_residuals: FloatArray
    iterations: int
    boundary_residual_max_abs: float
    off_grid_relative_ode_residual: float
    minimum_radial_metric_function: float
    outer_condition: str
    completed: bool

    @property
    def amplitude(self) -> FloatArray:
        return self.state[0]

    @property
    def radial_derivative(self) -> FloatArray:
        return self.state[1]

    @property
    def mass(self) -> FloatArray:
        return self.state[2]

    @property
    def lapse_exponent(self) -> FloatArray:
        return self.state[3]

    @property
    def adaptive_nodes(self) -> int:
        return int(self.radius.size)

    @property
    def max_collocation_rms_residual(self) -> float:
        return float(np.max(self.collocation_rms_residuals, initial=0.0))

    @property
    def outer_mass(self) -> float:
        return float(self.mass[-1])

    @property
    def central_lapse_exponent(self) -> float:
        phi_second = _origin_coefficients(
            self.central_amplitude,
            float(self.lapse_exponent[0]),
            self.frequency,
            self.dimensionless_coupling,
        )[3]
        return float(self.lapse_exponent[0] - 0.5 * phi_second * self.origin_epsilon**2)


@dataclass(frozen=True)
class StaticSphericalScalarShootingSolution:
    """Independent origin-to-wall DOP853 shooting evidence."""

    ivp: IVPEvidence
    frequency: float
    central_lapse_exponent: float
    central_amplitude: float
    dimensionless_coupling: float
    origin_epsilon: float
    outer_radius: float
    root_success: bool
    root_function_evaluations: int
    boundary_residuals: FloatArray
    minimum_radial_metric_function: float


def _positive_finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _nonnegative_finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be nonnegative and finite")
    return result


def _frequency_from_parameter(parameter: float) -> float:
    frequency = float(expit(parameter))
    if not 0.0 < frequency < 1.0:
        raise NumericalFailure("fitted frequency left the open sub-gap interval")
    return frequency


def _parameter_from_frequency(frequency: float) -> float:
    value = _positive_finite(frequency, "frequency")
    if value >= 1.0:
        raise ValueError("frequency must be strictly below one")
    return float(np.log(value / (1.0 - value)))


def _origin_coefficients(
    central_amplitude: float,
    central_lapse_exponent: float,
    frequency: float,
    dimensionless_coupling: float,
) -> tuple[float, float, float, float]:
    amplitude = float(central_amplitude)
    phi_zero = float(central_lapse_exponent)
    omega = float(frequency)
    alpha = float(dimensionless_coupling)
    inverse_lapse_squared = float(np.exp(-2.0 * phi_zero))
    density = (
        0.25 * omega**2 * amplitude**2 * inverse_lapse_squared
        + 1.0
        - float(jv(0, amplitude))
    )
    pressure = (
        0.25 * omega**2 * amplitude**2 * inverse_lapse_squared
        - (1.0 - float(jv(0, amplitude)))
    )
    amplitude_second = (
        2.0 * float(jv(1, amplitude))
        - omega**2 * amplitude * inverse_lapse_squared
    ) / 3.0
    mass_cubic = alpha * density / 6.0
    phi_second = alpha * (density / 6.0 + pressure / 2.0)
    values = (amplitude_second, mass_cubic, density, phi_second)
    if not np.all(np.isfinite(values)):
        raise NumericalFailure("regular-origin coefficients became non-finite")
    return values


def regular_origin_numeric_state(
    central_amplitude: float,
    central_lapse_exponent: float,
    frequency: float,
    dimensionless_coupling: float,
    origin_epsilon: float,
) -> FloatArray:
    """Return the second-order/cubic regular Taylor state at ``epsilon``."""

    amplitude = _positive_finite(central_amplitude, "central_amplitude")
    phi_zero = float(central_lapse_exponent)
    if not np.isfinite(phi_zero):
        raise ValueError("central_lapse_exponent must be finite")
    omega = _positive_finite(frequency, "frequency")
    if omega >= 1.0:
        raise ValueError("frequency must be strictly below one")
    alpha = _nonnegative_finite(dimensionless_coupling, "dimensionless_coupling")
    epsilon = _positive_finite(origin_epsilon, "origin_epsilon")
    amplitude_second, mass_cubic, _density, phi_second = _origin_coefficients(
        amplitude, phi_zero, omega, alpha
    )
    return np.asarray(
        [
            amplitude + 0.5 * amplitude_second * epsilon**2,
            amplitude_second * epsilon,
            mass_cubic * epsilon**3,
            phi_zero + 0.5 * phi_second * epsilon**2,
        ],
        dtype=np.float64,
    )


def static_spherical_scalar_rhs(
    radius: ArrayLike,
    state: ArrayLike,
    frequency: float,
    dimensionless_coupling: float,
) -> FloatArray:
    """Evaluate the first-order reduced ODE without clipping the metric."""

    coordinate = np.asarray(radius, dtype=np.float64)
    values = np.asarray(state, dtype=np.float64)
    scalar_input = coordinate.ndim == 0
    coordinate = np.atleast_1d(coordinate)
    if values.ndim == 1:
        values = values[:, None]
    if values.shape != (4, coordinate.size):
        raise ValueError("state must have shape (4, radius.size)")
    if np.any(~np.isfinite(coordinate)) or np.any(coordinate <= 0.0):
        raise ValueError("radius must contain positive finite values")
    if np.any(~np.isfinite(values)):
        raise NumericalFailure("state contains non-finite values")
    omega = _positive_finite(frequency, "frequency")
    if omega >= 1.0:
        raise ValueError("frequency must be strictly below one")
    alpha = _nonnegative_finite(dimensionless_coupling, "dimensionless_coupling")

    amplitude, derivative, mass, phi = values
    radial_metric = 1.0 - 2.0 * mass / coordinate
    if np.any(radial_metric <= 0.0) or np.any(~np.isfinite(radial_metric)):
        raise NumericalFailure("radial metric function became nonpositive")
    inverse_lapse_squared = np.exp(-2.0 * phi)
    potential = 1.0 - jv(0, amplitude)
    density = (
        0.25 * omega**2 * amplitude**2 * inverse_lapse_squared
        + 0.25 * radial_metric * derivative**2
        + potential
    )
    pressure = (
        0.25 * radial_metric * derivative**2
        + 0.25 * omega**2 * amplitude**2 * inverse_lapse_squared
        - potential
    )
    mass_prime = 0.5 * alpha * coordinate**2 * density
    phi_prime = (
        mass + 0.5 * alpha * coordinate**3 * pressure
    ) / (coordinate * (coordinate - 2.0 * mass))
    radial_metric_prime = 2.0 * mass / coordinate**2 - 2.0 * mass_prime / coordinate
    amplitude_second = (
        (2.0 * jv(1, amplitude) - omega**2 * inverse_lapse_squared * amplitude)
        / radial_metric
        - derivative
        * (
            phi_prime
            + radial_metric_prime / (2.0 * radial_metric)
            + 2.0 / coordinate
        )
    )
    result = np.vstack([derivative, amplitude_second, mass_prime, phi_prime])
    if np.any(~np.isfinite(result)):
        raise NumericalFailure("reduced ODE returned non-finite values")
    return result[:, 0] if scalar_input else np.asarray(result, dtype=np.float64)


def finite_wall_boundary_residual(
    left_state: ArrayLike,
    right_state: ArrayLike,
    frequency: float,
    *,
    central_amplitude: float,
    dimensionless_coupling: float,
    origin_epsilon: float,
    outer_radius: float,
) -> FloatArray:
    """Return origin-series, Robin-tail and Schwarzschild-match residuals."""

    left = np.asarray(left_state, dtype=np.float64)
    right = np.asarray(right_state, dtype=np.float64)
    if left.shape != (4,) or right.shape != (4,):
        raise ValueError("boundary states must each contain four values")
    amplitude = _positive_finite(central_amplitude, "central_amplitude")
    alpha = _nonnegative_finite(dimensionless_coupling, "dimensionless_coupling")
    epsilon = _positive_finite(origin_epsilon, "origin_epsilon")
    outer = _positive_finite(outer_radius, "outer_radius")
    omega = _positive_finite(frequency, "frequency")
    if omega >= 1.0:
        raise ValueError("frequency must be strictly below one")

    phi_zero_estimate = float(left[3])
    amplitude_second, mass_cubic, _density, _phi_second = _origin_coefficients(
        amplitude, phi_zero_estimate, omega, alpha
    )
    radial_metric_outer = 1.0 - 2.0 * right[2] / outer
    if radial_metric_outer <= 0.0:
        raise NumericalFailure("outer radial metric function became nonpositive")
    tail_squared = (
        1.0 - omega**2 / radial_metric_outer
    ) / radial_metric_outer
    if tail_squared <= 0.0:
        raise NumericalFailure("finite-wall tail is not evanescent")
    tail_rate = float(np.sqrt(tail_squared))
    return np.asarray(
        [
            left[0] - 0.5 * epsilon * left[1] - amplitude,
            left[1] - epsilon * amplitude_second,
            left[2] - mass_cubic * epsilon**3,
            right[1] + (tail_rate + 1.0 / outer) * right[0],
            right[3] - 0.5 * np.log(radial_metric_outer),
        ],
        dtype=np.float64,
    )


def _solution_diagnostics(
    evidence: BVPEvidence,
    frequency: float,
    alpha: float,
    boundary_residual: FloatArray,
) -> tuple[float, float, float]:
    coordinate = evidence.coordinate
    state = evidence.state
    radial_metric = 1.0 - 2.0 * state[2] / coordinate
    minimum_metric = float(np.min(radial_metric))
    if minimum_metric <= 0.0:
        raise NumericalFailure("converged solution contains a nonpositive metric function")
    midpoints = 0.5 * (coordinate[:-1] + coordinate[1:])
    spline = CubicSpline(coordinate, state, axis=1)
    sampled_state = np.asarray(spline(midpoints), dtype=np.float64)
    spline_derivative = np.asarray(spline(midpoints, 1), dtype=np.float64)
    rhs = static_spherical_scalar_rhs(midpoints, sampled_state, frequency, alpha)
    scale = 1.0 + np.max(np.abs(rhs), axis=0)
    relative_residual = float(
        np.max(np.abs(spline_derivative - rhs) / scale[None, :], initial=0.0)
    )
    boundary_norm = float(np.max(np.abs(boundary_residual), initial=0.0))
    if not np.isfinite(relative_residual) or not np.isfinite(boundary_norm):
        raise NumericalFailure("solution diagnostics became non-finite")
    return boundary_norm, relative_residual, minimum_metric


def solve_static_spherical_scalar_bvp(
    *,
    central_amplitude: float,
    dimensionless_coupling: float,
    origin_epsilon: float,
    outer_radius: float,
    initial_mesh_points: int,
    tolerance: float,
    max_nodes: int = 100_000,
    frequency_guess: float = 0.9,
    initial_solution: StaticSphericalScalarBVPSolution | None = None,
) -> StaticSphericalScalarBVPSolution:
    """Solve one finite-wall averaged Einstein--scalar branch point."""

    amplitude = _positive_finite(central_amplitude, "central_amplitude")
    alpha = _nonnegative_finite(dimensionless_coupling, "dimensionless_coupling")
    epsilon = _positive_finite(origin_epsilon, "origin_epsilon")
    outer = _positive_finite(outer_radius, "outer_radius")
    if epsilon >= outer:
        raise ValueError("origin_epsilon must be below outer_radius")
    points = int(initial_mesh_points)
    if points != initial_mesh_points or points < 20:
        raise ValueError("initial_mesh_points must be an integer of at least 20")
    tol = _positive_finite(tolerance, "tolerance")
    omega_guess = _positive_finite(frequency_guess, "frequency_guess")
    if omega_guess >= 1.0:
        raise ValueError("frequency_guess must be strictly below one")

    coordinate = np.linspace(epsilon, outer, points, dtype=np.float64)
    if initial_solution is None:
        decay = float(np.sqrt(1.0 - omega_guess**2))
        profile = amplitude * np.exp(-decay * coordinate)
        state_guess = np.vstack(
            [
                profile,
                -decay * profile,
                np.zeros_like(coordinate),
                np.zeros_like(coordinate),
            ]
        )
    else:
        source = initial_solution
        normalized = (coordinate - epsilon) / (outer - epsilon)
        source_coordinate = source.origin_epsilon + normalized * (
            source.outer_radius - source.origin_epsilon
        )
        source_spline = CubicSpline(source.radius, source.state, axis=1)
        state_guess = np.asarray(source_spline(source_coordinate), dtype=np.float64)
        omega_guess = source.frequency

    def equations(x: FloatArray, state: FloatArray, fitted: FloatArray) -> FloatArray:
        return static_spherical_scalar_rhs(
            x, state, _frequency_from_parameter(float(fitted[0])), alpha
        )

    def boundary(left: FloatArray, right: FloatArray, fitted: FloatArray) -> FloatArray:
        return finite_wall_boundary_residual(
            left,
            right,
            _frequency_from_parameter(float(fitted[0])),
            central_amplitude=amplitude,
            dimensionless_coupling=alpha,
            origin_epsilon=epsilon,
            outer_radius=outer,
        )

    evidence = solve_bvp_evidence(
        equations,
        boundary,
        coordinate,
        state_guess,
        initial_parameters=[_parameter_from_frequency(omega_guess)],
        tolerance=tol,
        max_nodes=max_nodes,
    )
    if evidence.parameters is None or evidence.parameters.size != 1:
        raise NumericalFailure("parameterized BVP returned no fitted frequency")
    fitted_parameter = float(evidence.parameters[0])
    frequency = _frequency_from_parameter(fitted_parameter)
    boundary_values = finite_wall_boundary_residual(
        evidence.state[:, 0],
        evidence.state[:, -1],
        frequency,
        central_amplitude=amplitude,
        dimensionless_coupling=alpha,
        origin_epsilon=epsilon,
        outer_radius=outer,
    )
    boundary_norm, off_grid, minimum_metric = _solution_diagnostics(
        evidence, frequency, alpha, boundary_values
    )
    return StaticSphericalScalarBVPSolution(
        radius=np.asarray(evidence.coordinate, dtype=np.float64),
        state=np.asarray(evidence.state, dtype=np.float64),
        frequency=frequency,
        fitted_parameter=fitted_parameter,
        central_amplitude=amplitude,
        dimensionless_coupling=alpha,
        origin_epsilon=epsilon,
        outer_radius=outer,
        initial_mesh_points=points,
        tolerance=tol,
        collocation_rms_residuals=np.asarray(
            evidence.rms_residuals, dtype=np.float64
        ),
        iterations=evidence.iterations,
        boundary_residual_max_abs=boundary_norm,
        off_grid_relative_ode_residual=off_grid,
        minimum_radial_metric_function=minimum_metric,
        outer_condition="approximate_finite_wall_evanescent_Robin",
        completed=True,
    )


def shoot_static_spherical_scalar_bvp(
    reference: StaticSphericalScalarBVPSolution,
    *,
    tolerances: SolverTolerances = SolverTolerances(
        rtol=1.0e-10, atol=1.0e-12, max_step=0.05
    ),
    root_tolerance: float = 1.0e-10,
) -> StaticSphericalScalarShootingSolution:
    """Independently solve the same wall data by DOP853 plus two-root shooting."""

    root_tol = _positive_finite(root_tolerance, "root_tolerance")
    amplitude = reference.central_amplitude
    alpha = reference.dimensionless_coupling
    epsilon = reference.origin_epsilon
    outer = reference.outer_radius
    initial_parameter = np.asarray(
        [
            _parameter_from_frequency(reference.frequency),
            reference.central_lapse_exponent,
        ],
        dtype=np.float64,
    )

    def integrate(parameters: FloatArray) -> tuple[IVPEvidence, FloatArray]:
        omega = _frequency_from_parameter(float(parameters[0]))
        phi_zero = float(parameters[1])
        initial_state = regular_origin_numeric_state(
            amplitude, phi_zero, omega, alpha, epsilon
        )
        evidence = solve_ivp_evidence(
            lambda x, state: static_spherical_scalar_rhs(x, state, omega, alpha),
            (epsilon, outer),
            initial_state,
            tolerances=tolerances,
            method="DOP853",
        )
        wall = finite_wall_boundary_residual(
            evidence.state[:, 0],
            evidence.state[:, -1],
            omega,
            central_amplitude=amplitude,
            dimensionless_coupling=alpha,
            origin_epsilon=epsilon,
            outer_radius=outer,
        )
        return evidence, np.asarray(wall[3:], dtype=np.float64)

    def residual(parameters: FloatArray) -> FloatArray:
        _evidence, wall = integrate(parameters)
        return wall

    solved = root(residual, initial_parameter, method="hybr", tol=root_tol)
    if not solved.success:
        raise NumericalFailure(f"shooting root failed: {solved.message}")
    ivp, wall_residual = integrate(np.asarray(solved.x, dtype=np.float64))
    frequency = _frequency_from_parameter(float(solved.x[0]))
    metric = 1.0 - 2.0 * ivp.state[2] / ivp.time
    minimum_metric = float(np.min(metric))
    if minimum_metric <= 0.0:
        raise NumericalFailure("shooting solution contains a nonpositive metric function")
    return StaticSphericalScalarShootingSolution(
        ivp=ivp,
        frequency=frequency,
        central_lapse_exponent=float(solved.x[1]),
        central_amplitude=amplitude,
        dimensionless_coupling=alpha,
        origin_epsilon=epsilon,
        outer_radius=outer,
        root_success=bool(solved.success),
        root_function_evaluations=int(solved.nfev),
        boundary_residuals=np.asarray(wall_residual, dtype=np.float64),
        minimum_radial_metric_function=minimum_metric,
    )


def compare_static_spherical_scalar_solutions(
    reference: StaticSphericalScalarBVPSolution,
    comparison: StaticSphericalScalarBVPSolution,
    *,
    samples: int = 500,
) -> float:
    """Return the preregistered normalized max state difference on overlap."""

    count = int(samples)
    if count != samples or count < 20:
        raise ValueError("samples must be an integer of at least 20")
    left = max(reference.origin_epsilon, comparison.origin_epsilon)
    right = min(reference.outer_radius, comparison.outer_radius)
    if left >= right:
        raise ValueError("solutions have no common radial interval")
    coordinate = np.linspace(left, right, count, dtype=np.float64)
    reference_state = CubicSpline(reference.radius, reference.state, axis=1)(coordinate)
    comparison_state = CubicSpline(comparison.radius, comparison.state, axis=1)(coordinate)
    scale = 1.0 + np.max(np.abs(reference_state), axis=1)
    return float(
        np.max(np.abs(comparison_state - reference_state) / scale[:, None])
    )
