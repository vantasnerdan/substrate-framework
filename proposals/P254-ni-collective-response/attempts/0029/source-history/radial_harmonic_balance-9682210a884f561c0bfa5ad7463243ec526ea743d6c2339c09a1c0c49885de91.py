"""Finite-harmonic radial sine-Gordon balance with explicit tail semantics.

For the dimensionless radial equation

``u_tt - u_rr - 2*u_r/r + sin(u) = 0``

and an odd-cosine ansatz ``u = sum_n a_n(r) cos(n*omega*t)``, projection gives

``a_n'' + 2*a_n'/r + (n*omega)**2*a_n - S_n[a] = 0``.

The fundamental can be evanescent when ``omega < 1``.  Higher harmonics are
not automatically localized: every channel with ``n*omega > 1`` has an
oscillatory ``1/r`` tail.  The solver below therefore labels its outer
conditions.  A Dirichlet condition on a radiative channel defines a finite-box
standing wave; it is not evidence for a finite-energy infinite-domain
breather.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import expit

from .numerics import BVPEvidence, NumericalFailure, solve_bvp_evidence

FloatArray = NDArray[np.float64]
TailBehavior = Literal["evanescent", "threshold", "radiative"]
OuterPolicy = Literal["channel_aware", "dirichlet_box"]


@dataclass(frozen=True)
class HarmonicTailChannel:
    """Linear far-field classification of one temporal harmonic."""

    harmonic: int
    angular_frequency: float
    behavior: TailBehavior
    radial_rate: float


@dataclass(frozen=True)
class RadialHarmonicBalanceSolution:
    """A converged finite-harmonic, finite-radius BVP branch point."""

    radius: FloatArray
    harmonics: tuple[int, ...]
    amplitudes: FloatArray
    radial_derivatives: FloatArray
    frequency: float
    central_fundamental: float
    origin_epsilon: float
    outer_radius: float
    temporal_samples: int
    outer_policy: OuterPolicy
    outer_conditions: tuple[str, ...]
    collocation_rms_residuals: FloatArray
    iterations: int
    fitted_parameter: float
    completed: bool

    @property
    def max_collocation_rms_residual(self) -> float:
        """Largest interval residual reported by the collocation solver."""

        return float(np.max(self.collocation_rms_residuals, initial=0.0))

    @property
    def central_amplitudes(self) -> FloatArray:
        """Second-order extrapolation of every even radial mode to ``r=0``."""

        return np.asarray(
            self.amplitudes[:, 0]
            - 0.5 * self.origin_epsilon * self.radial_derivatives[:, 0],
            dtype=np.float64,
        )

    @property
    def tail_channels(self) -> tuple[HarmonicTailChannel, ...]:
        """Return the final solution's linear far-field channel labels."""

        return classify_harmonic_tail_channels(self.harmonics, self.frequency)


def odd_harmonics(max_harmonic: int) -> tuple[int, ...]:
    """Return all positive odd harmonics through ``max_harmonic``."""

    maximum = int(max_harmonic)
    if maximum != max_harmonic or maximum < 1 or maximum % 2 == 0:
        raise ValueError("max_harmonic must be a positive odd integer")
    return tuple(range(1, maximum + 1, 2))


def reconstruct_radial_harmonics(
    amplitudes: ArrayLike,
    harmonics: tuple[int, ...],
    phase: ArrayLike,
) -> FloatArray:
    """Reconstruct ``u(r,phase)`` with phase rows and radial columns."""

    coefficients = _amplitude_matrix(amplitudes, harmonics)
    phases = _finite_vector(phase, "phase")
    basis = np.cos(np.outer(np.asarray(harmonics, dtype=np.float64), phases))
    return np.asarray(basis.T @ coefficients, dtype=np.float64)


def project_sine_harmonics(
    amplitudes: ArrayLike,
    harmonics: tuple[int, ...],
    *,
    target_harmonics: tuple[int, ...] | None = None,
    temporal_samples: int = 256,
) -> FloatArray:
    """Project ``sin(sum a_n cos(n*tau))`` onto cosine harmonics.

    The periodic trapezoidal/DFT normalization is ``2/N``.  No NumPy
    integration alias is required, so the calculation is stable across the
    NumPy 1.26 and 2.x API generations supported by the framework.
    """

    source = _validate_harmonics(harmonics)
    targets = source if target_harmonics is None else _validate_harmonics(
        target_harmonics, require_fundamental=False
    )
    coefficients = _amplitude_matrix(amplitudes, source)
    samples = _temporal_sample_count(temporal_samples, source + targets)
    tau = 2.0 * np.pi * np.arange(samples, dtype=np.float64) / samples
    source_basis = np.cos(np.outer(np.asarray(source, dtype=np.float64), tau))
    target_basis = np.cos(np.outer(np.asarray(targets, dtype=np.float64), tau))
    field = source_basis.T @ coefficients
    return np.asarray(
        (2.0 / samples) * target_basis @ np.sin(field), dtype=np.float64
    )


def nonlinear_projection_remainder(
    amplitudes: ArrayLike,
    harmonics: tuple[int, ...],
    *,
    temporal_samples: int = 512,
) -> FloatArray:
    """Return the unretained time-dependent part of ``sin(u)``.

    Rows are temporal phase samples and columns are radial samples.  This is
    the load-bearing truncation residual after the retained projected ODEs have
    converged; a small frequency change does not bound it.
    """

    modes = _validate_harmonics(harmonics)
    coefficients = _amplitude_matrix(amplitudes, modes)
    samples = _temporal_sample_count(temporal_samples, modes)
    tau = 2.0 * np.pi * np.arange(samples, dtype=np.float64) / samples
    basis = np.cos(np.outer(np.asarray(modes, dtype=np.float64), tau))
    field = basis.T @ coefficients
    retained = (2.0 / samples) * basis @ np.sin(field)
    return np.asarray(np.sin(field) - basis.T @ retained, dtype=np.float64)


def sampled_harmonic_balance_residual(
    radius: ArrayLike,
    amplitudes: ArrayLike,
    harmonics: tuple[int, ...],
    frequency: float,
    *,
    temporal_samples: int = 512,
) -> FloatArray:
    """Evaluate retained projected ODE residuals from sampled amplitudes.

    Radial derivatives are reconstructed independently with second-order
    finite differences.  This is an off-collocation diagnostic, not the BVP
    solver's own residual estimate.
    """

    modes = _validate_harmonics(harmonics)
    coordinate = _radial_coordinate(radius)
    coefficients = _amplitude_matrix(amplitudes, modes, coordinate.size)
    omega = _subgap_frequency_value(frequency)
    first = np.gradient(coefficients, coordinate, axis=1, edge_order=2)
    second = np.gradient(first, coordinate, axis=1, edge_order=2)
    laplacian = np.empty_like(coefficients)
    if coordinate[0] == 0.0:
        laplacian[:, 0] = 3.0 * second[:, 0]
        laplacian[:, 1:] = (
            second[:, 1:] + 2.0 * first[:, 1:] / coordinate[None, 1:]
        )
    else:
        laplacian = second + 2.0 * first / coordinate[None, :]
    nonlinear = project_sine_harmonics(
        coefficients,
        modes,
        temporal_samples=temporal_samples,
    )
    squared_frequencies = np.square(np.asarray(modes, dtype=np.float64) * omega)
    return np.asarray(
        laplacian + squared_frequencies[:, None] * coefficients - nonlinear,
        dtype=np.float64,
    )


def classify_harmonic_tail_channels(
    harmonics: tuple[int, ...],
    frequency: float,
    *,
    threshold_tolerance: float = 1.0e-12,
) -> tuple[HarmonicTailChannel, ...]:
    """Classify far-field channels of the linearized radial equation.

    ``radial_rate`` is the positive decay rate for an evanescent channel, zero
    at threshold, and the positive radial wavenumber for a radiative channel.
    """

    modes = _validate_harmonics(harmonics, require_fundamental=False)
    omega = _positive_finite(frequency, "frequency")
    tolerance = _positive_finite(threshold_tolerance, "threshold_tolerance")
    channels: list[HarmonicTailChannel] = []
    for harmonic in modes:
        angular = harmonic * omega
        difference = 1.0 - angular**2
        if difference > tolerance:
            behavior: TailBehavior = "evanescent"
            rate = np.sqrt(difference)
        elif difference < -tolerance:
            behavior = "radiative"
            rate = np.sqrt(-difference)
        else:
            behavior = "threshold"
            rate = 0.0
        channels.append(
            HarmonicTailChannel(
                harmonic=harmonic,
                angular_frequency=float(angular),
                behavior=behavior,
                radial_rate=float(rate),
            )
        )
    return tuple(channels)


def solve_radial_harmonic_balance(
    harmonics: tuple[int, ...],
    *,
    central_fundamental: float,
    outer_radius: float,
    frequency_guess: float,
    radial_points: int = 400,
    temporal_samples: int = 192,
    tolerance: float = 1.0e-7,
    max_nodes: int = 50_000,
    origin_epsilon: float = 1.0e-3,
    outer_policy: OuterPolicy = "channel_aware",
    initial_solution: RadialHarmonicBalanceSolution | None = None,
) -> RadialHarmonicBalanceSolution:
    """Solve a finite-harmonic radial BVP with a fitted sub-gap frequency.

    The central fundamental amplitude is a declared branch coordinate; it is
    not determined by the equation.  ``channel_aware`` imposes a decaying
    Robin condition on evanescent modes and a Dirichlet wall on radiative
    modes.  ``dirichlet_box`` places every retained mode at zero on the wall.
    Either policy remains a finite-radius model whenever a radiative channel
    is present.
    """

    modes = _validate_harmonics(harmonics)
    amplitude = _positive_finite(central_fundamental, "central_fundamental")
    outer = _positive_finite(outer_radius, "outer_radius")
    epsilon = _positive_finite(origin_epsilon, "origin_epsilon")
    if epsilon >= outer:
        raise ValueError("origin_epsilon must be below outer_radius")
    omega_guess = _subgap_frequency_value(frequency_guess)
    points = int(radial_points)
    if points != radial_points or points < 20:
        raise ValueError("radial_points must be an integer of at least 20")
    samples = _temporal_sample_count(temporal_samples, modes)
    if outer_policy not in {"channel_aware", "dirichlet_box"}:
        raise ValueError("outer_policy must be channel_aware or dirichlet_box")

    tau = 2.0 * np.pi * np.arange(samples, dtype=np.float64) / samples
    basis = np.cos(np.outer(np.asarray(modes, dtype=np.float64), tau))

    def nonlinear_projection(coefficients: FloatArray) -> FloatArray:
        field = basis.T @ coefficients
        return np.asarray(
            (2.0 / samples) * basis @ np.sin(field), dtype=np.float64
        )

    if initial_solution is None:
        coordinate = np.linspace(epsilon, outer, points)
        state_guess = np.zeros((2 * len(modes), points), dtype=np.float64)
        inverse_width = np.sqrt(1.0 - omega_guess**2)
        state_guess[0] = amplitude * np.exp(-inverse_width * coordinate)
        state_guess[1] = -inverse_width * state_guess[0]
    else:
        if not np.isclose(initial_solution.outer_radius, outer, rtol=0.0, atol=1.0e-12):
            raise ValueError("initial_solution must use the requested outer_radius")
        if not np.isclose(
            initial_solution.origin_epsilon, epsilon, rtol=0.0, atol=1.0e-15
        ):
            raise ValueError("initial_solution must use the requested origin_epsilon")
        if any(mode not in modes for mode in initial_solution.harmonics):
            raise ValueError("initial_solution harmonics must be a subset of harmonics")
        coordinate = np.asarray(initial_solution.radius, dtype=np.float64)
        state_guess = np.zeros(
            (2 * len(modes), coordinate.size), dtype=np.float64
        )
        for old_index, mode in enumerate(initial_solution.harmonics):
            new_index = modes.index(mode)
            state_guess[2 * new_index] = initial_solution.amplitudes[old_index]
            state_guess[2 * new_index + 1] = (
                initial_solution.radial_derivatives[old_index]
            )
        omega_guess = initial_solution.frequency

    def equations(
        coordinate_values: FloatArray,
        state: FloatArray,
        fitted: FloatArray,
    ) -> FloatArray:
        omega = _frequency_from_parameter(float(fitted[0]))
        coefficients = state[0::2]
        nonlinear = nonlinear_projection(coefficients)
        derivative = np.empty_like(state)
        derivative[0::2] = state[1::2]
        for index, harmonic in enumerate(modes):
            derivative[2 * index + 1] = (
                -2.0 * state[2 * index + 1] / coordinate_values
                - (harmonic * omega) ** 2 * state[2 * index]
                + nonlinear[index]
            )
        return derivative

    def boundary_residual(
        left: FloatArray,
        right: FloatArray,
        fitted: FloatArray,
    ) -> FloatArray:
        omega = _frequency_from_parameter(float(fitted[0]))
        nonlinear_left = nonlinear_projection(left[0::2, None])[:, 0]
        residuals: list[float] = []
        for index, harmonic in enumerate(modes):
            curvature = (
                nonlinear_left[index] - (harmonic * omega) ** 2 * left[2 * index]
            ) / 3.0
            residuals.append(left[2 * index + 1] - epsilon * curvature)
        channels = classify_harmonic_tail_channels(modes, omega)
        for index, channel in enumerate(channels):
            if outer_policy == "channel_aware" and channel.behavior == "evanescent":
                residuals.append(
                    right[2 * index + 1]
                    + (channel.radial_rate + 1.0 / outer) * right[2 * index]
                )
            else:
                residuals.append(right[2 * index])
        residuals.append(left[0] - 0.5 * epsilon * left[1] - amplitude)
        return np.asarray(residuals, dtype=np.float64)

    parameter_guess = _parameter_from_frequency(omega_guess)
    evidence: BVPEvidence = solve_bvp_evidence(
        equations,
        boundary_residual,
        coordinate,
        state_guess,
        initial_parameters=[parameter_guess],
        tolerance=tolerance,
        max_nodes=max_nodes,
    )
    if evidence.parameters is None or evidence.parameters.size != 1:
        raise NumericalFailure("parameterized BVP returned no fitted frequency")
    fitted_parameter = float(evidence.parameters[0])
    frequency = _frequency_from_parameter(fitted_parameter)
    channels = classify_harmonic_tail_channels(modes, frequency)
    conditions = tuple(
        "decaying_robin"
        if outer_policy == "channel_aware" and channel.behavior == "evanescent"
        else "dirichlet_box"
        for channel in channels
    )
    return RadialHarmonicBalanceSolution(
        radius=np.asarray(evidence.coordinate, dtype=np.float64),
        harmonics=modes,
        amplitudes=np.asarray(evidence.state[0::2], dtype=np.float64),
        radial_derivatives=np.asarray(evidence.state[1::2], dtype=np.float64),
        frequency=frequency,
        central_fundamental=amplitude,
        origin_epsilon=epsilon,
        outer_radius=outer,
        temporal_samples=samples,
        outer_policy=outer_policy,
        outer_conditions=conditions,
        collocation_rms_residuals=np.asarray(
            evidence.rms_residuals, dtype=np.float64
        ),
        iterations=evidence.iterations,
        fitted_parameter=fitted_parameter,
        completed=True,
    )


def _frequency_from_parameter(parameter: float) -> float:
    frequency = float(expit(parameter))
    if not 0.0 < frequency < 1.0:
        raise NumericalFailure("fitted frequency left the open sub-gap interval")
    return frequency


def _parameter_from_frequency(frequency: float) -> float:
    omega = _subgap_frequency_value(frequency)
    return float(np.log(omega / (1.0 - omega)))


def _subgap_frequency_value(value: float) -> float:
    frequency = _positive_finite(value, "frequency")
    if frequency >= 1.0:
        raise ValueError("frequency must be strictly below the linear threshold one")
    return frequency


def _validate_harmonics(
    harmonics: tuple[int, ...],
    *,
    require_fundamental: bool = True,
) -> tuple[int, ...]:
    modes = tuple(int(value) for value in harmonics)
    if not modes or modes != harmonics:
        raise ValueError("harmonics must be a nonempty tuple of integers")
    if any(mode <= 0 or mode % 2 == 0 for mode in modes):
        raise ValueError("harmonics must be positive odd integers")
    if tuple(sorted(set(modes))) != modes:
        raise ValueError("harmonics must be strictly increasing and unique")
    if require_fundamental and modes[0] != 1:
        raise ValueError("the solved harmonic set must begin with the fundamental")
    return modes


def _amplitude_matrix(
    amplitudes: ArrayLike,
    harmonics: tuple[int, ...],
    radial_size: int | None = None,
) -> FloatArray:
    array = np.asarray(amplitudes, dtype=np.float64)
    if array.ndim == 1:
        array = array[:, None]
    if array.ndim != 2 or array.shape[0] != len(harmonics):
        raise ValueError("amplitudes must have one row per harmonic")
    if radial_size is not None and array.shape[1] != radial_size:
        raise ValueError("amplitudes must have one column per radial sample")
    if not np.all(np.isfinite(array)):
        raise ValueError("amplitudes must contain only finite values")
    return array


def _temporal_sample_count(
    temporal_samples: int,
    harmonics: tuple[int, ...],
) -> int:
    samples = int(temporal_samples)
    if samples != temporal_samples or samples < 4 * (max(harmonics) + 1):
        raise ValueError(
            "temporal_samples must be an integer at least four times the largest harmonic plus one"
        )
    return samples


def _radial_coordinate(values: ArrayLike) -> FloatArray:
    radius = _finite_vector(values, "radius")
    if radius.size < 5 or radius[0] < 0.0 or np.any(np.diff(radius) <= 0.0):
        raise ValueError(
            "radius must contain at least five strictly increasing nonnegative samples"
        )
    return radius


def _finite_vector(values: ArrayLike, name: str) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite one-dimensional array")
    return array


def _positive_finite(value: float, name: str) -> float:
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number
