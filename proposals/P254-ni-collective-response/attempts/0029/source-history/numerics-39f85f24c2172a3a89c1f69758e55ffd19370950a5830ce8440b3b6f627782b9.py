"""Reusable SciPy-backed numerical evidence for physics claims.

These helpers centralize solver failure handling and the evidence that campaign
verifiers need to inspect.  They deliberately do not assign epistemic status:
an IVP or PDE solution remains resolution-bounded numerical evidence until the
claim-specific convergence, conservation, sensitivity, and comparison gates
have passed.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_bvp, solve_ivp
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter

FloatArray = NDArray[np.float64]


class NumericalFailure(RuntimeError):
    """Raised when a numerical method cannot produce valid evidence."""


@dataclass(frozen=True)
class SolverTolerances:
    """Explicit adaptive-solver tolerances shared by IVP and PDE checks."""

    rtol: float = 1.0e-9
    atol: float = 1.0e-12
    max_step: float = np.inf

    def __post_init__(self) -> None:
        if not np.isfinite(self.rtol) or self.rtol <= 0.0:
            raise ValueError("rtol must be positive and finite")
        if not np.isfinite(self.atol) or self.atol <= 0.0:
            raise ValueError("atol must be positive and finite")
        if self.max_step <= 0.0 or np.isnan(self.max_step):
            raise ValueError("max_step must be positive")


@dataclass(frozen=True)
class IVPEvidence:
    """The solution and diagnostic evidence returned by an IVP integration."""

    time: FloatArray
    state: FloatArray
    method: str
    function_evaluations: int
    invariant_values: FloatArray | None
    max_abs_invariant_drift: float | None


@dataclass(frozen=True)
class BVPEvidence:
    """The converged mesh, state, parameters, and residuals for a BVP."""

    coordinate: FloatArray
    state: FloatArray
    rms_residuals: FloatArray
    iterations: int
    parameters: FloatArray | None = None

    @property
    def max_rms_residual(self) -> float:
        """Return the largest interval residual, or zero for an empty mesh."""

        return float(np.max(self.rms_residuals, initial=0.0))


@dataclass(frozen=True)
class RefinementEvidence:
    """Errors and empirical orders from a claim-defined refinement study."""

    resolutions: tuple[int, ...]
    errors: tuple[float, ...]
    observed_orders: tuple[float | None, ...]

    @property
    def errors_strictly_decrease(self) -> bool:
        """Whether every refinement reduced the claim-defined error."""

        return all(fine < coarse for coarse, fine in zip(self.errors, self.errors[1:]))

    @property
    def final_error(self) -> float:
        return self.errors[-1]


def _real_vector(values: ArrayLike, *, name: str) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional real array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _trapezoid_implementation() -> Callable[..., Any]:
    implementation = getattr(np, "trapezoid", None)
    if implementation is None:
        implementation = getattr(np, "trapz", None)
    if implementation is None:
        raise NumericalFailure("NumPy provides neither trapezoid nor trapz")
    return implementation


def trapezoid_integral(values: ArrayLike, coordinate: ArrayLike) -> float:
    """Integrate sampled values across supported NumPy API generations.

    NumPy 2 renamed ``trapz`` to ``trapezoid`` and later removed the legacy
    alias, while this package also supports NumPy 1.26.  Resolve the current
    name first and use the legacy name only when the installed version needs
    it.  Keeping the dispatch here prevents scientific modules from carrying
    inconsistent version probes.
    """

    result = float(_trapezoid_implementation()(values, coordinate))
    if not np.isfinite(result):
        raise NumericalFailure("trapezoidal integration returned a non-finite value")
    return result


def trapezoid_integral_axis(
    values: ArrayLike,
    coordinate: ArrayLike,
    *,
    axis: int = -1,
) -> FloatArray:
    """Integrate an array along one axis with the shared NumPy dispatch."""

    samples = np.asarray(values, dtype=np.float64)
    points = _real_vector(coordinate, name="coordinate")
    if samples.ndim == 0:
        raise ValueError("values must have at least one dimension")
    normalized_axis = axis if axis >= 0 else samples.ndim + axis
    if normalized_axis < 0 or normalized_axis >= samples.ndim:
        raise ValueError(
            f"axis must select one of {samples.ndim} dimensions; got {axis}"
        )
    if samples.shape[normalized_axis] != points.size:
        raise ValueError("coordinate length must match the integration axis")
    if not np.all(np.isfinite(samples)):
        raise ValueError("values must contain only finite values")
    result = np.asarray(
        _trapezoid_implementation()(samples, points, axis=normalized_axis),
        dtype=np.float64,
    )
    if not np.all(np.isfinite(result)):
        raise NumericalFailure("trapezoidal integration returned non-finite values")
    return result


def local_polynomial_time_derivative(
    time: ArrayLike,
    values: ArrayLike,
    derivative_order: int,
    *,
    window_duration: float,
    polynomial_order: int = 5,
) -> FloatArray:
    """Differentiate a uniformly sampled trace by local polynomial fitting.

    ``window_duration`` is converted to the nearest symmetric odd-point window,
    so refinement studies can preserve a physical support width while changing
    the sample interval. SciPy's Savitzky--Golay implementation supplies the
    local least-squares polynomial derivative. Endpoint extrapolations are
    returned, but a scientific claim must exclude a justified boundary window
    rather than interpreting them as interior evidence.
    """

    samples, data, spacing = _uniform_time_data(time, values)
    order = _nonnegative_integer(derivative_order, "derivative_order")
    degree = _nonnegative_integer(polynomial_order, "polynomial_order")
    if degree < order:
        raise ValueError("polynomial_order must be at least derivative_order")
    duration = float(window_duration)
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("window_duration must be positive and finite")
    half_intervals = max(1, int(round(duration / (2.0 * spacing))))
    window_points = 2 * half_intervals + 1
    if window_points <= degree:
        raise ValueError("window_duration gives too few points for polynomial_order")
    if window_points > samples.size:
        raise ValueError("window_duration exceeds the sampled trace")
    derivative = savgol_filter(
        data,
        window_length=window_points,
        polyorder=degree,
        deriv=order,
        delta=spacing,
        axis=-1,
        mode="interp",
    )
    if not np.all(np.isfinite(derivative)):
        raise NumericalFailure("local polynomial derivative became non-finite")
    return np.asarray(derivative, dtype=np.float64)


def interpolating_spline_time_derivative(
    time: ArrayLike,
    values: ArrayLike,
    derivative_order: int,
    *,
    spline_degree: int = 5,
) -> FloatArray:
    """Differentiate sampled data with an independent interpolating B-spline.

    Unlike :func:`local_polynomial_time_derivative`, this route performs no
    local least-squares smoothing. Agreement between the two therefore probes
    whether a reported derivative is set by the resolved trace rather than by
    one estimator's filtering choice. Endpoint values still require exclusion
    in claim-level evidence.
    """

    samples, data, _spacing = _uniform_time_data(time, values)
    order = _nonnegative_integer(derivative_order, "derivative_order")
    degree = _nonnegative_integer(spline_degree, "spline_degree")
    if degree < 1:
        raise ValueError("spline_degree must be positive")
    if order > degree:
        raise ValueError("derivative_order cannot exceed spline_degree")
    if samples.size <= degree:
        raise ValueError("sampled trace has too few points for spline_degree")
    spline = make_interp_spline(samples, data, k=degree, axis=-1)
    derivative = spline.derivative(order)(samples)
    if not np.all(np.isfinite(derivative)):
        raise NumericalFailure("interpolating spline derivative became non-finite")
    return np.asarray(derivative, dtype=np.float64)


def _uniform_time_data(
    time: ArrayLike,
    values: ArrayLike,
) -> tuple[FloatArray, FloatArray, float]:
    samples = _real_vector(time, name="time")
    if samples.size < 3 or np.any(np.diff(samples) <= 0.0):
        raise ValueError("time must contain at least three strictly increasing samples")
    steps = np.diff(samples)
    spacing = float(np.mean(steps))
    if not np.allclose(steps, spacing, rtol=1.0e-10, atol=1.0e-13):
        raise ValueError("time must be uniformly sampled")
    data = np.asarray(values, dtype=np.float64)
    if data.ndim == 0 or data.shape[-1] != samples.size:
        raise ValueError("values must have time along its last axis")
    if not np.all(np.isfinite(data)):
        raise ValueError("values must contain only finite samples")
    return samples, data, spacing


def _nonnegative_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or int(value) != value or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def solve_ivp_evidence(
    rhs: Callable[[float, FloatArray], ArrayLike],
    time_span: tuple[float, float],
    initial_state: ArrayLike,
    *,
    sample_times: ArrayLike | None = None,
    tolerances: SolverTolerances = SolverTolerances(),
    method: str = "DOP853",
    invariant: Callable[[FloatArray], float] | None = None,
    **solver_options: Any,
) -> IVPEvidence:
    """Solve a real-valued IVP and retain diagnostics needed by a verifier.

    ``invariant`` is evaluated on every returned sample.  The helper reports
    absolute drift; the claim verifier must set a physically justified pass
    threshold and should also refine tolerances and sampling independently.
    """

    state0 = _real_vector(initial_state, name="initial_state")
    times = None if sample_times is None else _real_vector(sample_times, name="sample_times")
    if not np.all(np.isfinite(time_span)) or time_span[1] <= time_span[0]:
        raise ValueError("time_span must contain finite, increasing endpoints")

    result = solve_ivp(
        rhs,
        time_span,
        state0,
        method=method,
        t_eval=times,
        rtol=tolerances.rtol,
        atol=tolerances.atol,
        max_step=tolerances.max_step,
        **solver_options,
    )
    if not result.success:
        raise NumericalFailure(f"IVP integration failed: {result.message}")
    if not np.all(np.isfinite(result.y)):
        raise NumericalFailure("IVP integration returned non-finite state values")

    invariant_values: FloatArray | None = None
    invariant_drift: float | None = None
    if invariant is not None:
        invariant_values = np.asarray(
            [invariant(result.y[:, index]) for index in range(result.y.shape[1])],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(invariant_values)):
            raise NumericalFailure("invariant evaluation returned non-finite values")
        invariant_drift = float(np.max(np.abs(invariant_values - invariant_values[0])))

    return IVPEvidence(
        time=np.asarray(result.t, dtype=np.float64),
        state=np.asarray(result.y, dtype=np.float64),
        method=method,
        function_evaluations=int(result.nfev),
        invariant_values=invariant_values,
        max_abs_invariant_drift=invariant_drift,
    )


def solve_method_of_lines(
    spatial_rhs: Callable[[float, FloatArray], ArrayLike],
    time_span: tuple[float, float],
    initial_state: ArrayLike,
    **options: Any,
) -> IVPEvidence:
    """Integrate a caller-discretized PDE in time with IVP evidence capture.

    The spatial discretization, boundary conditions, mesh, and refinement study
    remain explicit responsibilities of the claim implementation.  This avoids
    presenting a generic time integrator as a proof of a particular PDE model.
    """

    return solve_ivp_evidence(spatial_rhs, time_span, initial_state, **options)


def solve_bvp_evidence(
    equations: Callable[..., ArrayLike],
    boundary_residual: Callable[..., ArrayLike],
    coordinate: ArrayLike,
    state_guess: ArrayLike,
    *,
    initial_parameters: ArrayLike | None = None,
    tolerance: float = 1.0e-6,
    max_nodes: int = 10_000,
    **solver_options: Any,
) -> BVPEvidence:
    """Solve a two-point BVP and retain collocation and parameter evidence.

    When ``initial_parameters`` is supplied, ``equations`` and
    ``boundary_residual`` must accept SciPy's additional parameter argument.
    Existing parameter-free callers retain their two-argument signatures.
    """

    x = _real_vector(coordinate, name="coordinate")
    guess = np.asarray(state_guess, dtype=np.float64)
    if guess.ndim != 2 or guess.shape[1] != x.size:
        raise ValueError("state_guess must have shape (state dimension, coordinate size)")
    if not np.all(np.isfinite(guess)):
        raise ValueError("state_guess must contain only finite values")
    if tolerance <= 0.0 or not np.isfinite(tolerance):
        raise ValueError("tolerance must be positive and finite")
    if max_nodes < x.size:
        raise ValueError("max_nodes cannot be smaller than the initial mesh")
    parameter_guess = None
    if initial_parameters is not None:
        parameter_guess = _real_vector(initial_parameters, name="initial_parameters")

    result = solve_bvp(
        equations,
        boundary_residual,
        x,
        guess,
        p=parameter_guess,
        tol=tolerance,
        max_nodes=max_nodes,
        **solver_options,
    )
    if not result.success:
        raise NumericalFailure(f"BVP solve failed: {result.message}")
    if not np.all(np.isfinite(result.y)):
        raise NumericalFailure("BVP solve returned non-finite state values")

    return BVPEvidence(
        coordinate=np.asarray(result.x, dtype=np.float64),
        state=np.asarray(result.y, dtype=np.float64),
        rms_residuals=np.asarray(result.rms_residuals, dtype=np.float64),
        iterations=int(result.niter),
        parameters=(
            None
            if result.p is None
            else np.asarray(result.p, dtype=np.float64)
        ),
    )


def refinement_study(
    resolutions: Sequence[int],
    solve: Callable[[int], Any],
    error: Callable[[Any], float],
) -> RefinementEvidence:
    """Run a resolution study using a claim-defined solution and error metric.

    Resolution is interpreted as inverse characteristic spacing, so empirical
    order is ``log(error_coarse/error_fine) / log(N_fine/N_coarse)``.  The
    caller must justify that interpretation for its discretization and choose
    an error metric appropriate to the physical claim.
    """

    levels = tuple(int(value) for value in resolutions)
    if len(levels) < 3:
        raise ValueError("a refinement study requires at least three resolutions")
    if any(level <= 0 for level in levels) or any(
        fine <= coarse for coarse, fine in zip(levels, levels[1:])
    ):
        raise ValueError("resolutions must be positive and strictly increasing")

    errors = tuple(float(error(solve(level))) for level in levels)
    if any(value < 0.0 or not np.isfinite(value) for value in errors):
        raise NumericalFailure("refinement error must be finite and nonnegative")

    orders: list[float | None] = []
    for coarse_n, fine_n, coarse_error, fine_error in zip(
        levels, levels[1:], errors, errors[1:]
    ):
        if coarse_error == 0.0 or fine_error == 0.0:
            orders.append(None)
        else:
            orders.append(
                float(
                    np.log(coarse_error / fine_error)
                    / np.log(float(fine_n) / float(coarse_n))
                )
            )

    return RefinementEvidence(levels, errors, tuple(orders))
