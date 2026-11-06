"""Conditional one-dimensional reflected first-passage constructions.

The functions in this module do not infer a physical bath, mobility, barrier,
or event channel.  They solve or sample a separately declared overdamped Itô
diffusion on a finite interval with a reflecting left endpoint and an absorbing
right endpoint.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import IntegrationWarning, quad

from .numerics import NumericalFailure

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]
ScalarPotential = Callable[[float], float]


@dataclass(frozen=True)
class FirstPassageQuadratureEvidence:
    """Adaptive-quadrature evidence for one mean first-passage time."""

    initial_coordinate: float
    reflecting_boundary: float
    absorbing_boundary: float
    thermal_scale: float
    friction: float
    mean_first_passage_time: float
    outer_absolute_error: float
    maximum_inner_absolute_error: float
    inner_evaluations: int
    epsabs: float
    epsrel: float

    @property
    def inverse_mean_first_passage_time(self) -> float:
        """Return ``1/E[T]``; this is not automatically a constant hazard."""

        if self.mean_first_passage_time == 0.0:
            return math.inf
        return 1.0 / self.mean_first_passage_time


@dataclass(frozen=True)
class CensoredFirstPassageSummary:
    """Finite-horizon statistics with completed and censored paths separated."""

    horizon: float
    trajectory_count: int
    completed_count: int
    completion_fraction: float
    completed_only_mean: float | None
    inverse_completed_only_mean: float | None
    restricted_mean: float
    inverse_restricted_mean: float


@dataclass(frozen=True)
class ReflectedEscapeEnsemble:
    """Euler-Maruyama paths and the exact discrete conventions that made them."""

    event_times: FloatArray
    completed: BoolArray
    final_coordinates: FloatArray
    reflecting_boundary: float
    absorbing_boundary: float
    initial_coordinate: float
    thermal_scale: float
    friction: float
    time_step: float
    horizon: float
    trajectory_count: int
    seed: int
    steps_executed: int
    reflection_rule: str = "absolute_overshoot"

    @property
    def summary(self) -> CensoredFirstPassageSummary:
        """Return the explicit completed-only and restricted-mean statistics."""

        return summarize_censored_first_passage(
            self.event_times,
            self.completed,
            self.horizon,
        )


def _finite_real(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_finite(value: float, name: str) -> float:
    result = _finite_real(value, name)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def reflected_absorbing_backward_residual(
    thermal_scale: float,
    friction: float,
    potential_gradient: float,
    first_derivative: float,
    second_derivative: float,
) -> float:
    r"""Return the residual of ``Theta*tau''-U'*tau'=-gamma``."""

    theta = _positive_finite(thermal_scale, "thermal_scale")
    gamma = _positive_finite(friction, "friction")
    gradient = _finite_real(potential_gradient, "potential_gradient")
    first = _finite_real(first_derivative, "first_derivative")
    second = _finite_real(second_derivative, "second_derivative")
    return theta * second - gradient * first + gamma


def reflected_absorbing_mfpt(
    potential: ScalarPotential,
    initial_coordinate: float,
    reflecting_boundary: float,
    absorbing_boundary: float,
    thermal_scale: float,
    friction: float,
    *,
    epsabs: float = 1.0e-11,
    epsrel: float = 1.0e-11,
    limit: int = 200,
) -> FirstPassageQuadratureEvidence:
    r"""Evaluate the exact reflected-to-absorbing MFPT integral.

    For the declared interior generator

    ``L f = -(U'/gamma) f' + (Theta/gamma) f''``,

    reflection at ``a`` and absorption at ``b``, the backward solution is

    ``tau(x)=(gamma/Theta) int_x^b exp(U(y)/Theta)
    int_a^y exp(-U(z)/Theta) dz dy``.

    The potential is shifted internally by ``U(a)``.  This improves numerical
    conditioning and makes the exact invariance under an additive potential
    constant explicit.  The caller remains responsible for establishing the
    regularity and physical provenance of the supplied potential.
    """

    if not callable(potential):
        raise TypeError("potential must be callable")
    initial = _finite_real(initial_coordinate, "initial_coordinate")
    reflecting = _finite_real(reflecting_boundary, "reflecting_boundary")
    absorbing = _finite_real(absorbing_boundary, "absorbing_boundary")
    theta = _positive_finite(thermal_scale, "thermal_scale")
    gamma = _positive_finite(friction, "friction")
    absolute_tolerance = _positive_finite(epsabs, "epsabs")
    relative_tolerance = _positive_finite(epsrel, "epsrel")
    if absorbing <= reflecting:
        raise ValueError("absorbing_boundary must exceed reflecting_boundary")
    if initial < reflecting or initial > absorbing:
        raise ValueError("initial_coordinate must lie in the closed interval")
    if isinstance(limit, bool) or int(limit) != limit or limit < 20:
        raise ValueError("limit must be an integer at least 20")
    if initial == absorbing:
        return FirstPassageQuadratureEvidence(
            initial_coordinate=initial,
            reflecting_boundary=reflecting,
            absorbing_boundary=absorbing,
            thermal_scale=theta,
            friction=gamma,
            mean_first_passage_time=0.0,
            outer_absolute_error=0.0,
            maximum_inner_absolute_error=0.0,
            inner_evaluations=0,
            epsabs=absolute_tolerance,
            epsrel=relative_tolerance,
        )

    reference = _finite_real(potential(reflecting), "potential(reflecting_boundary)")
    maximum_inner_error = 0.0
    inner_evaluations = 0

    def shifted_potential(coordinate: float) -> float:
        value = _finite_real(potential(float(coordinate)), "potential(coordinate)")
        return (value - reference) / theta

    def inner_integral(outer_coordinate: float) -> float:
        nonlocal maximum_inner_error, inner_evaluations
        value, error = quad(
            lambda coordinate: math.exp(-shifted_potential(coordinate)),
            reflecting,
            outer_coordinate,
            epsabs=absolute_tolerance,
            epsrel=relative_tolerance,
            limit=int(limit),
        )
        inner_evaluations += 1
        maximum_inner_error = max(maximum_inner_error, float(error))
        if not math.isfinite(value) or not math.isfinite(error):
            raise NumericalFailure("inner first-passage quadrature became non-finite")
        return float(value)

    def outer_integrand(coordinate: float) -> float:
        try:
            result = math.exp(shifted_potential(coordinate)) * inner_integral(coordinate)
        except OverflowError as exc:
            raise NumericalFailure("first-passage exponential overflowed") from exc
        if not math.isfinite(result):
            raise NumericalFailure("first-passage integrand became non-finite")
        return result

    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        try:
            outer_value, outer_error = quad(
                outer_integrand,
                initial,
                absorbing,
                epsabs=absolute_tolerance,
                epsrel=relative_tolerance,
                limit=int(limit),
            )
        except IntegrationWarning as exc:
            raise NumericalFailure("first-passage quadrature did not converge") from exc

    mean_time = gamma * float(outer_value) / theta
    scaled_outer_error = gamma * float(outer_error) / theta
    if not math.isfinite(mean_time) or mean_time <= 0.0:
        raise NumericalFailure("first-passage quadrature did not return a positive finite time")
    return FirstPassageQuadratureEvidence(
        initial_coordinate=initial,
        reflecting_boundary=reflecting,
        absorbing_boundary=absorbing,
        thermal_scale=theta,
        friction=gamma,
        mean_first_passage_time=mean_time,
        outer_absolute_error=scaled_outer_error,
        maximum_inner_absolute_error=maximum_inner_error,
        inner_evaluations=inner_evaluations,
        epsabs=absolute_tolerance,
        epsrel=relative_tolerance,
    )


def free_reflected_absorbing_mfpt(
    interval_length: float,
    thermal_scale: float,
    friction: float,
) -> float:
    r"""Return ``gamma*L**2/(2*Theta)`` for zero potential and a start at reflection."""

    length = _positive_finite(interval_length, "interval_length")
    theta = _positive_finite(thermal_scale, "thermal_scale")
    gamma = _positive_finite(friction, "friction")
    return gamma * length**2 / (2.0 * theta)


def linear_potential_reflected_absorbing_mfpt(
    force: float,
    interval_length: float,
    thermal_scale: float,
    friction: float,
) -> float:
    r"""Return the exact start-point MFPT for ``U=F*(x-a)``.

    With ``s=F*L/Theta``, the result is
    ``gamma*Theta*(exp(s)-1-s)/F**2``.  A series protects the removable
    zero-force limit from cancellation.
    """

    slope = _finite_real(force, "force")
    length = _positive_finite(interval_length, "interval_length")
    theta = _positive_finite(thermal_scale, "thermal_scale")
    gamma = _positive_finite(friction, "friction")
    if slope == 0.0:
        return free_reflected_absorbing_mfpt(length, theta, gamma)
    scaled = slope * length / theta
    if abs(scaled) < 1.0e-3:
        numerator = sum(scaled**power / math.factorial(power) for power in range(2, 10))
    else:
        try:
            numerator = math.expm1(scaled) - scaled
        except OverflowError as exc:
            raise NumericalFailure("linear-potential MFPT overflowed") from exc
    result = gamma * theta * numerator / slope**2
    if not math.isfinite(result) or result <= 0.0:
        raise NumericalFailure("linear-potential MFPT is not positive and finite")
    return result


def quadratic_barrier_potential(coordinate: float, barrier_energy: float) -> float:
    r"""Return the declared scaled barrier ``E*(2*rho-rho**2)``."""

    rho = _finite_real(coordinate, "coordinate")
    energy = _positive_finite(barrier_energy, "barrier_energy")
    return energy * (2.0 * rho - rho**2)


def quadratic_barrier_gradient(coordinate: float, barrier_energy: float) -> float:
    r"""Return ``dU/drho=2*E*(1-rho)`` for the scaled barrier."""

    rho = _finite_real(coordinate, "coordinate")
    energy = _positive_finite(barrier_energy, "barrier_energy")
    return 2.0 * energy * (1.0 - rho)


def summarize_censored_first_passage(
    event_times: ArrayLike,
    completed: ArrayLike,
    horizon: float,
) -> CensoredFirstPassageSummary:
    """Separate completed-only and administratively restricted time means."""

    cutoff = _positive_finite(horizon, "horizon")
    times = np.asarray(event_times, dtype=np.float64)
    flags = np.asarray(completed, dtype=np.bool_)
    if times.ndim != 1 or flags.ndim != 1 or times.shape != flags.shape or times.size == 0:
        raise ValueError("event_times and completed must be equal nonempty vectors")
    completed_times = times[flags]
    if not np.all(np.isfinite(completed_times)):
        raise ValueError("completed event times must be finite")
    if np.any(completed_times < 0.0) or np.any(completed_times > cutoff):
        raise ValueError("completed event times must lie within the horizon")
    if np.any(np.isfinite(times[~flags])):
        raise ValueError("censored event times must be NaN")

    count = int(times.size)
    completed_count = int(np.count_nonzero(flags))
    fraction = completed_count / count
    restricted = np.where(flags, times, cutoff)
    restricted_mean = float(np.mean(restricted))
    completed_mean: float | None = None
    inverse_completed: float | None = None
    if completed_count:
        completed_mean = float(np.mean(completed_times))
        if completed_mean <= 0.0:
            raise ValueError("completed event-time mean must be positive")
        inverse_completed = 1.0 / completed_mean
    return CensoredFirstPassageSummary(
        horizon=cutoff,
        trajectory_count=count,
        completed_count=completed_count,
        completion_fraction=fraction,
        completed_only_mean=completed_mean,
        inverse_completed_only_mean=inverse_completed,
        restricted_mean=restricted_mean,
        inverse_restricted_mean=1.0 / restricted_mean,
    )


def thresholded_completed_only_rate(
    summary: CensoredFirstPassageSummary,
    minimum_completion_fraction: float = 0.05,
) -> float:
    """Reproduce a declared completion-threshold classifier, not a physical zero."""

    threshold = _finite_real(minimum_completion_fraction, "minimum_completion_fraction")
    if threshold < 0.0 or threshold > 1.0:
        raise ValueError("minimum_completion_fraction must lie in [0,1]")
    if summary.completion_fraction < threshold:
        return 0.0
    if summary.inverse_completed_only_mean is None:
        return 0.0
    return summary.inverse_completed_only_mean


def simulate_reflected_euler_maruyama(
    potential_gradient: ScalarPotential,
    initial_coordinate: float,
    reflecting_boundary: float,
    absorbing_boundary: float,
    thermal_scale: float,
    friction: float,
    *,
    trajectory_count: int,
    time_step: float,
    horizon: float,
    seed: int,
) -> ReflectedEscapeEnsemble:
    """Sample a declared reflected Euler-Maruyama discretization.

    A proposal crossing below the left boundary is mapped to its absolute
    overshoot about that boundary.  Absorption is tested after this reflection.
    This is a discrete numerical convention whose continuum interpretation
    requires timestep refinement.
    """

    if not callable(potential_gradient):
        raise TypeError("potential_gradient must be callable")
    initial = _finite_real(initial_coordinate, "initial_coordinate")
    reflecting = _finite_real(reflecting_boundary, "reflecting_boundary")
    absorbing = _finite_real(absorbing_boundary, "absorbing_boundary")
    theta = _positive_finite(thermal_scale, "thermal_scale")
    gamma = _positive_finite(friction, "friction")
    step = _positive_finite(time_step, "time_step")
    cutoff = _positive_finite(horizon, "horizon")
    if absorbing <= reflecting:
        raise ValueError("absorbing_boundary must exceed reflecting_boundary")
    if initial < reflecting or initial >= absorbing:
        raise ValueError("initial_coordinate must lie in the half-open interval")
    if isinstance(trajectory_count, bool) or int(trajectory_count) != trajectory_count or trajectory_count <= 0:
        raise ValueError("trajectory_count must be a positive integer")
    if isinstance(seed, bool) or int(seed) != seed or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    step_count_float = cutoff / step
    step_count = int(round(step_count_float))
    if not math.isclose(step_count_float, step_count, rel_tol=1.0e-12, abs_tol=1.0e-12):
        raise ValueError("horizon must be an integer multiple of time_step")

    count = int(trajectory_count)
    random = np.random.default_rng(int(seed))
    coordinates = np.full(count, initial, dtype=np.float64)
    completed = np.zeros(count, dtype=np.bool_)
    event_times = np.full(count, np.nan, dtype=np.float64)
    noise_scale = math.sqrt(2.0 * theta * step / gamma)
    steps_executed = 0

    for step_index in range(1, step_count + 1):
        live_indices = np.flatnonzero(~completed)
        if live_indices.size == 0:
            break
        live_coordinates = coordinates[live_indices]
        try:
            gradients = np.asarray(
                potential_gradient(live_coordinates),  # type: ignore[arg-type]
                dtype=np.float64,
            )
        except (TypeError, ValueError):
            gradients = np.asarray(
                [potential_gradient(float(value)) for value in live_coordinates],
                dtype=np.float64,
            )
        if gradients.ndim == 0:
            gradients = np.full(live_indices.shape, float(gradients), dtype=np.float64)
        if gradients.shape != live_indices.shape or not np.all(np.isfinite(gradients)):
            raise NumericalFailure("potential_gradient returned invalid samples")
        proposals = (
            coordinates[live_indices]
            - gradients * step / gamma
            + noise_scale * random.standard_normal(live_indices.size)
        )
        proposals = reflecting + np.abs(proposals - reflecting)
        if not np.all(np.isfinite(proposals)):
            raise NumericalFailure("Euler-Maruyama coordinates became non-finite")
        coordinates[live_indices] = proposals
        hits = proposals >= absorbing
        if np.any(hits):
            hit_indices = live_indices[hits]
            completed[hit_indices] = True
            event_times[hit_indices] = step_index * step
        steps_executed = step_index

    return ReflectedEscapeEnsemble(
        event_times=event_times,
        completed=completed,
        final_coordinates=coordinates,
        reflecting_boundary=reflecting,
        absorbing_boundary=absorbing,
        initial_coordinate=initial,
        thermal_scale=theta,
        friction=gamma,
        time_step=step,
        horizon=cutoff,
        trajectory_count=count,
        seed=int(seed),
        steps_executed=steps_executed,
    )
