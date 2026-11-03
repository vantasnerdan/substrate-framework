"""Numerical evidence for the dimensionless 3+1 radial sine-Gordon IVP.

The module solves

``u_tt = u_rr + 2 u_r / r - sin(u)``

with even regularity at the origin and a homogeneous Dirichlet condition at a
finite outer radius.  Its results are finite-grid, finite-time simulation
evidence; no API here asserts an exact or eternally periodic breather.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import find_peaks

from .numerics import (
    NumericalFailure,
    SolverTolerances,
    solve_method_of_lines,
    trapezoid_integral,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class RadialEvolution:
    """Recorded diagnostics from one radial sine-Gordon evolution."""

    radius: FloatArray
    time: FloatArray
    center: FloatArray
    core_energy: FloatArray
    total_energy: FloatArray
    core_energy_radius_moment: FloatArray
    total_energy_radius_moment: FloatArray
    final_field: FloatArray
    final_velocity: FloatArray
    spacing: float
    timestep: float | None
    outer_radius: float
    core_radius: float
    method: str
    function_evaluations: int | None
    completed: bool
    max_monitor_amplitude: float


@dataclass(frozen=True)
class FrequencyEvidence:
    """Two estimates of a settled trace's fundamental angular frequency."""

    spectral_omega: float
    crossing_omega: float
    fft_bin_width: float
    crossing_cycles: int
    window_start: float
    window_end: float

    @property
    def estimator_disagreement(self) -> float:
        """Absolute disagreement between the FFT and crossing estimates."""

        return abs(self.spectral_omega - self.crossing_omega)


@dataclass(frozen=True)
class PeakFrequencyEvidence:
    """Time-domain frequency evidence from successive prominent maxima."""

    angular_frequency: float
    mean_period: float
    relative_period_standard_deviation: float
    cycles: int
    window_start: float
    window_end: float


def gaussian_radial_seed(radius: ArrayLike, amplitude: float, width: float) -> FloatArray:
    """Return ``amplitude * exp(-(r / width)**2)`` on a radial grid."""

    r = _radial_coordinate(radius)
    if not np.isfinite(amplitude):
        raise ValueError("amplitude must be finite")
    if not np.isfinite(width) or width <= 0.0:
        raise ValueError("width must be positive and finite")
    return np.asarray(amplitude * np.exp(-np.square(r / width)), dtype=np.float64)


def radial_laplacian(
    field: ArrayLike,
    spacing: float,
    *,
    geometric_coefficient: float = 2.0,
) -> FloatArray:
    """Second-order ``u_rr + geometric_coefficient*u_r/r`` on a uniform grid.

    The coordinate is understood to be ``r_j = j * spacing``.  Even radial
    regularity gives the origin limit ``(1 + geometric_coefficient) u_rr(0)``;
    the corresponding stencil is
    ``2 * (1 + geometric_coefficient) * (u[1] - u[0]) / spacing**2``.
    The returned outer entry is zero because the evolution fixes that point as
    a Dirichlet boundary rather than integrating its acceleration.
    """

    u = _field_vector(field)
    h = _positive_finite(spacing, "spacing")
    if not np.isfinite(geometric_coefficient):
        raise ValueError("geometric_coefficient must be finite")
    radius = h * np.arange(u.size, dtype=np.float64)
    laplacian = np.empty_like(u)
    laplacian[1:-1] = (
        (u[2:] - 2.0 * u[1:-1] + u[:-2]) / h**2
        + geometric_coefficient
        * (u[2:] - u[:-2])
        / (2.0 * h * radius[1:-1])
    )
    laplacian[0] = (
        2.0
        * (1.0 + geometric_coefficient)
        * (u[1] - u[0])
        / h**2
    )
    laplacian[-1] = 0.0
    return laplacian


def radial_gradient(field: ArrayLike, spacing: float) -> FloatArray:
    """Return a second-order radial derivative with regular origin value zero."""

    u = _field_vector(field)
    h = _positive_finite(spacing, "spacing")
    gradient = np.empty_like(u)
    gradient[0] = 0.0
    gradient[1:-1] = (u[2:] - u[:-2]) / (2.0 * h)
    gradient[-1] = (3.0 * u[-1] - 4.0 * u[-2] + u[-3]) / (2.0 * h)
    return gradient


def radial_sine_gordon_energy(
    field: ArrayLike,
    velocity: ArrayLike,
    radius: ArrayLike,
) -> float:
    """Evaluate ``4*pi*integral r^2 T00 dr`` by trapezoidal quadrature.

    The shared numerical helper selects NumPy's current ``trapezoid`` name and
    falls back to ``trapz`` only on supported older NumPy versions.
    """

    u = _field_vector(field)
    v = _matching_vector(velocity, u.size, "velocity")
    r = _radial_coordinate(radius, expected_size=u.size)
    density = _radial_sine_gordon_energy_density(u, v, r)
    return 4.0 * np.pi * trapezoid_integral(np.square(r) * density, r)


def radial_sine_gordon_energy_radius_moment(
    field: ArrayLike,
    velocity: ArrayLike,
    radius: ArrayLike,
) -> float:
    """Evaluate ``integral T00*r^2*d^3x = 4*pi*integral r^4*T00 dr``."""

    u = _field_vector(field)
    v = _matching_vector(velocity, u.size, "velocity")
    r = _radial_coordinate(radius, expected_size=u.size)
    density = _radial_sine_gordon_energy_density(u, v, r)
    return 4.0 * np.pi * trapezoid_integral(np.power(r, 4) * density, r)


def estimate_angular_frequency(
    time: ArrayLike,
    trace: ArrayLike,
    *,
    window_start: float,
    detrend_linear: bool = False,
) -> FrequencyEvidence:
    """Estimate angular frequency by a windowed FFT and rising crossings.

    The FFT uses a Hann window and quadratic peak interpolation.  The separate
    crossing estimate averages complete positive-slope periods of the
    mean-subtracted trace.  A claim verifier must compare their disagreement
    with the FFT resolution and repeat the estimate on another time window.
    """

    t = _one_dimensional_finite(time, "time")
    y = _matching_vector(trace, t.size, "trace")
    if t.size < 16 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must contain at least 16 strictly increasing samples")
    steps = np.diff(t)
    sample_step = float(np.mean(steps))
    if not np.allclose(steps, sample_step, rtol=1.0e-8, atol=1.0e-12):
        raise ValueError("frequency estimation requires uniform sampling")
    selected = t >= window_start
    if np.count_nonzero(selected) < 16:
        raise ValueError("frequency window must contain at least 16 samples")
    tw = t[selected]
    centered = np.array(y[selected], copy=True)
    if detrend_linear:
        centered -= np.polyval(np.polyfit(tw, centered, 1), tw)
    else:
        centered -= float(np.mean(centered))
    if np.max(np.abs(centered)) == 0.0:
        raise ValueError("trace has no oscillatory component")

    windowed = centered * np.hanning(centered.size)
    spectrum = np.abs(np.fft.rfft(windowed))
    frequencies = 2.0 * np.pi * np.fft.rfftfreq(centered.size, sample_step)
    peak = int(np.argmax(spectrum[1:]) + 1)
    interpolated_peak = float(peak)
    if 1 <= peak < spectrum.size - 1:
        left, middle, right = np.log(np.maximum(spectrum[peak - 1 : peak + 2], 1.0e-300))
        curvature = left - 2.0 * middle + right
        if curvature != 0.0:
            interpolated_peak += float(0.5 * (left - right) / curvature)
    bin_width = float(frequencies[1] - frequencies[0])
    spectral_omega = interpolated_peak * bin_width

    rising = np.flatnonzero((centered[:-1] <= 0.0) & (centered[1:] > 0.0))
    crossing_times: list[float] = []
    for index in rising:
        y0, y1 = centered[index], centered[index + 1]
        fraction = float(-y0 / (y1 - y0))
        crossing_times.append(float(tw[index] + fraction * (tw[index + 1] - tw[index])))
    if len(crossing_times) < 3:
        raise NumericalFailure("fewer than two complete crossing periods were resolved")
    crossing_omega = float(2.0 * np.pi / np.mean(np.diff(crossing_times)))
    return FrequencyEvidence(
        spectral_omega=spectral_omega,
        crossing_omega=crossing_omega,
        fft_bin_width=bin_width,
        crossing_cycles=len(crossing_times) - 1,
        window_start=float(tw[0]),
        window_end=float(tw[-1]),
    )


def estimate_peak_angular_frequency(
    time: ArrayLike,
    trace: ArrayLike,
    *,
    window_start: float,
    minimum_period: float,
    prominence_fraction: float = 0.1,
    detrend_linear: bool = True,
) -> PeakFrequencyEvidence:
    """Estimate angular frequency from prominent time-domain maxima.

    ``minimum_period`` is a declared physical separation scale that prevents
    high harmonics or sampling noise from being counted as extra cycles.  The
    function optionally removes only a least-squares linear drift and reports
    the period scatter rather than hiding it in a single mean.
    """

    t = _one_dimensional_finite(time, "time")
    y = _matching_vector(trace, t.size, "trace")
    if t.size < 16 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must contain at least 16 strictly increasing samples")
    steps = np.diff(t)
    sample_step = float(np.mean(steps))
    if not np.allclose(steps, sample_step, rtol=1.0e-8, atol=1.0e-12):
        raise ValueError("peak frequency estimation requires uniform sampling")
    period_floor = _positive_finite(minimum_period, "minimum_period")
    if (
        not np.isfinite(prominence_fraction)
        or prominence_fraction <= 0.0
        or prominence_fraction >= 1.0
    ):
        raise ValueError("prominence_fraction must lie strictly between zero and one")
    selected = t >= window_start
    if np.count_nonzero(selected) < 16:
        raise ValueError("frequency window must contain at least 16 samples")
    tw = t[selected]
    analyzed = np.array(y[selected], copy=True)
    if detrend_linear:
        analyzed -= np.polyval(np.polyfit(tw, analyzed, 1), tw)
    else:
        analyzed -= float(np.mean(analyzed))
    peak_to_peak = float(np.ptp(analyzed))
    if peak_to_peak == 0.0:
        raise ValueError("trace has no oscillatory component")
    minimum_samples = max(1, int(np.floor(period_floor / sample_step)))
    peaks, _ = find_peaks(
        analyzed,
        distance=minimum_samples,
        prominence=prominence_fraction * peak_to_peak,
    )
    if peaks.size < 3:
        raise NumericalFailure("fewer than two complete peak periods were resolved")
    periods = np.diff(tw[peaks])
    mean_period = float(np.mean(periods))
    relative_scatter = float(np.std(periods) / mean_period)
    return PeakFrequencyEvidence(
        angular_frequency=float(2.0 * np.pi / mean_period),
        mean_period=mean_period,
        relative_period_standard_deviation=relative_scatter,
        cycles=int(periods.size),
        window_start=float(tw[0]),
        window_end=float(tw[-1]),
    )


def evolve_radial_sine_gordon_leapfrog(
    *,
    amplitude: float,
    width: float,
    spacing: float,
    outer_radius: float,
    final_time: float,
    core_radius: float = 20.0,
    courant: float = 0.4,
    sample_interval: float = 0.2,
    damping_width: float = 0.0,
    damping_strength: float = 1.0,
    geometric_coefficient: float = 2.0,
) -> RadialEvolution:
    """Evolve a Gaussian-at-rest IVP using a centered explicit leapfrog.

    A centered damping term is applied only in the optional outer sponge.  The
    recorded state and centered velocity share the recorded timestamp.
    Non-finite values raise :class:`NumericalFailure` instead of returning a
    nominally successful trajectory.
    """

    h, radius, core_mask, monitor_index = _evolution_grid(
        spacing, outer_radius, core_radius, damping_width
    )
    duration = _positive_finite(final_time, "final_time")
    cfl = _positive_finite(courant, "courant")
    if cfl >= 1.0:
        raise ValueError("courant must be below one for the explicit radial scheme")
    requested_dt = cfl * h
    steps = int(np.ceil(duration / requested_dt))
    dt = duration / steps
    record_stride = max(1, int(round(_positive_finite(sample_interval, "sample_interval") / dt)))
    if not np.isfinite(damping_strength) or damping_strength < 0.0:
        raise ValueError("damping_strength must be finite and nonnegative")

    field_previous = gaussian_radial_seed(radius, amplitude, width)
    field_previous[-1] = 0.0
    acceleration = radial_laplacian(
        field_previous, h, geometric_coefficient=geometric_coefficient
    ) - np.sin(field_previous)
    field_current = field_previous + 0.5 * dt**2 * acceleration
    field_current[-1] = 0.0

    damping = np.zeros_like(radius)
    if damping_width > 0.0:
        start = outer_radius - damping_width
        selected = radius > start
        damping[selected] = damping_strength * np.square(
            (radius[selected] - start) / damping_width
        )

    times: list[float] = []
    centers: list[float] = []
    core_energies: list[float] = []
    total_energies: list[float] = []
    core_energy_radius_moments: list[float] = []
    total_energy_radius_moments: list[float] = []
    max_monitor = max(abs(float(field_previous[monitor_index])), abs(float(field_current[monitor_index])))
    final_velocity = np.zeros_like(field_current)

    for current_step in range(1, steps):
        acceleration = radial_laplacian(
            field_current, h, geometric_coefficient=geometric_coefficient
        ) - np.sin(field_current)
        numerator = (
            2.0 * field_current
            - (1.0 - 0.5 * damping * dt) * field_previous
            + dt**2 * acceleration
        )
        field_next = numerator / (1.0 + 0.5 * damping * dt)
        field_next[-1] = 0.0
        if not np.all(np.isfinite(field_next)):
            raise NumericalFailure(
                f"leapfrog produced non-finite values at step {current_step + 1}"
            )
        centered_velocity = (field_next - field_previous) / (2.0 * dt)
        max_monitor = max(max_monitor, abs(float(field_current[monitor_index])))
        if current_step % record_stride == 0:
            _record_diagnostics(
                radius,
                field_current,
                centered_velocity,
                current_step * dt,
                core_mask,
                times,
                centers,
                core_energies,
                total_energies,
                core_energy_radius_moments,
                total_energy_radius_moments,
            )
        field_previous, field_current = field_current, field_next
        final_velocity = centered_velocity

    if len(times) < 2:
        raise NumericalFailure("evolution recorded fewer than two diagnostic samples")
    return RadialEvolution(
        radius=radius,
        time=np.asarray(times, dtype=np.float64),
        center=np.asarray(centers, dtype=np.float64),
        core_energy=np.asarray(core_energies, dtype=np.float64),
        total_energy=np.asarray(total_energies, dtype=np.float64),
        core_energy_radius_moment=np.asarray(
            core_energy_radius_moments, dtype=np.float64
        ),
        total_energy_radius_moment=np.asarray(
            total_energy_radius_moments, dtype=np.float64
        ),
        # ``final_velocity`` is centered on ``field_previous`` after the last
        # swap; return the state at that same time rather than mixing levels.
        final_field=field_previous,
        final_velocity=final_velocity,
        spacing=h,
        timestep=dt,
        outer_radius=float(radius[-1]),
        core_radius=core_radius,
        method="centered-leapfrog",
        function_evaluations=None,
        completed=True,
        max_monitor_amplitude=max_monitor,
    )


def evolve_radial_sine_gordon_mol(
    *,
    amplitude: float,
    width: float,
    spacing: float,
    outer_radius: float,
    final_time: float,
    core_radius: float = 20.0,
    sample_interval: float = 0.2,
    tolerances: SolverTolerances = SolverTolerances(
        rtol=2.0e-7, atol=2.0e-9, max_step=0.1
    ),
    method: str = "DOP853",
) -> RadialEvolution:
    """Independently integrate the closed-box semidiscrete IVP with SciPy."""

    h, radius, core_mask, monitor_index = _evolution_grid(
        spacing, outer_radius, core_radius, 0.0
    )
    duration = _positive_finite(final_time, "final_time")
    sample = _positive_finite(sample_interval, "sample_interval")
    sample_count = int(np.floor(duration / sample)) + 1
    sample_times = np.linspace(0.0, duration, sample_count)
    field0 = gaussian_radial_seed(radius, amplitude, width)
    field0[-1] = 0.0
    interior_size = radius.size - 1
    state0 = np.concatenate((field0[:-1], np.zeros(interior_size)))

    def rhs(_time: float, state: FloatArray) -> FloatArray:
        field = np.zeros_like(radius)
        field[:-1] = state[:interior_size]
        velocity = state[interior_size:]
        acceleration = radial_laplacian(field, h)[:-1] - np.sin(field[:-1])
        return np.concatenate((velocity, acceleration))

    solution = solve_method_of_lines(
        rhs,
        (0.0, duration),
        state0,
        sample_times=sample_times,
        tolerances=tolerances,
        method=method,
    )
    times: list[float] = []
    centers: list[float] = []
    core_energies: list[float] = []
    total_energies: list[float] = []
    core_energy_radius_moments: list[float] = []
    total_energy_radius_moments: list[float] = []
    max_monitor = 0.0
    final_field = np.zeros_like(radius)
    final_velocity = np.zeros_like(radius)
    for index, time_value in enumerate(solution.time):
        field = np.zeros_like(radius)
        velocity = np.zeros_like(radius)
        field[:-1] = solution.state[:interior_size, index]
        velocity[:-1] = solution.state[interior_size:, index]
        _record_diagnostics(
            radius,
            field,
            velocity,
            float(time_value),
            core_mask,
            times,
            centers,
            core_energies,
            total_energies,
            core_energy_radius_moments,
            total_energy_radius_moments,
        )
        max_monitor = max(max_monitor, abs(float(field[monitor_index])))
        final_field, final_velocity = field, velocity

    return RadialEvolution(
        radius=radius,
        time=np.asarray(times, dtype=np.float64),
        center=np.asarray(centers, dtype=np.float64),
        core_energy=np.asarray(core_energies, dtype=np.float64),
        total_energy=np.asarray(total_energies, dtype=np.float64),
        core_energy_radius_moment=np.asarray(
            core_energy_radius_moments, dtype=np.float64
        ),
        total_energy_radius_moment=np.asarray(
            total_energy_radius_moments, dtype=np.float64
        ),
        final_field=final_field,
        final_velocity=final_velocity,
        spacing=h,
        timestep=None,
        outer_radius=float(radius[-1]),
        core_radius=core_radius,
        method=solution.method,
        function_evaluations=solution.function_evaluations,
        completed=True,
        max_monitor_amplitude=max_monitor,
    )


def _record_diagnostics(
    radius: FloatArray,
    field: FloatArray,
    velocity: FloatArray,
    time_value: float,
    core_mask: NDArray[np.bool_],
    times: list[float],
    centers: list[float],
    core_energies: list[float],
    total_energies: list[float],
    core_energy_radius_moments: list[float],
    total_energy_radius_moments: list[float],
) -> None:
    times.append(time_value)
    centers.append(float(field[0]))
    total_energies.append(radial_sine_gordon_energy(field, velocity, radius))
    core_energies.append(
        radial_sine_gordon_energy(field[core_mask], velocity[core_mask], radius[core_mask])
    )
    total_energy_radius_moments.append(
        radial_sine_gordon_energy_radius_moment(field, velocity, radius)
    )
    core_energy_radius_moments.append(
        radial_sine_gordon_energy_radius_moment(
            field[core_mask], velocity[core_mask], radius[core_mask]
        )
    )


def _radial_sine_gordon_energy_density(
    field: FloatArray,
    velocity: FloatArray,
    radius: FloatArray,
) -> FloatArray:
    spacing = _uniform_spacing(radius)
    gradient = radial_gradient(field, spacing)
    return np.asarray(
        0.5 * np.square(velocity)
        + 0.5 * np.square(gradient)
        + 1.0
        - np.cos(field),
        dtype=np.float64,
    )


def _evolution_grid(
    spacing: float,
    outer_radius: float,
    core_radius: float,
    damping_width: float,
) -> tuple[float, FloatArray, NDArray[np.bool_], int]:
    h = _positive_finite(spacing, "spacing")
    outer = _positive_finite(outer_radius, "outer_radius")
    core = _positive_finite(core_radius, "core_radius")
    if core >= outer:
        raise ValueError("core_radius must be smaller than outer_radius")
    if not np.isfinite(damping_width) or damping_width < 0.0 or damping_width >= outer:
        raise ValueError("damping_width must be finite, nonnegative, and below outer_radius")
    intervals = int(round(outer / h))
    if intervals < 3 or not np.isclose(intervals * h, outer, rtol=0.0, atol=1.0e-12):
        raise ValueError("outer_radius must be an integer multiple of spacing")
    radius = h * np.arange(intervals + 1, dtype=np.float64)
    core_mask = radius <= core
    monitor_radius = outer - damping_width if damping_width > 0.0 else outer - 2.0 * h
    monitor_index = int(np.argmin(np.abs(radius - monitor_radius)))
    return h, radius, core_mask, monitor_index


def _positive_finite(value: float, name: str) -> float:
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _one_dimensional_finite(values: ArrayLike, name: str) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite one-dimensional array")
    return array


def _field_vector(values: ArrayLike) -> FloatArray:
    array = _one_dimensional_finite(values, "field")
    if array.size < 3:
        raise ValueError("field must contain at least three grid points")
    return array


def _matching_vector(values: ArrayLike, expected_size: int, name: str) -> FloatArray:
    array = _one_dimensional_finite(values, name)
    if array.size != expected_size:
        raise ValueError(f"{name} must have the same size as the field")
    return array


def _radial_coordinate(values: ArrayLike, expected_size: int | None = None) -> FloatArray:
    radius = _one_dimensional_finite(values, "radius")
    if expected_size is not None and radius.size != expected_size:
        raise ValueError("radius must have the same size as the field")
    if radius.size < 3 or radius[0] != 0.0 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("radius must start at zero and be strictly increasing")
    return radius


def _uniform_spacing(radius: FloatArray) -> float:
    differences = np.diff(radius)
    spacing = float(differences[0])
    if not np.allclose(differences, spacing, rtol=1.0e-10, atol=1.0e-13):
        raise ValueError("radius must be uniformly spaced")
    return spacing
