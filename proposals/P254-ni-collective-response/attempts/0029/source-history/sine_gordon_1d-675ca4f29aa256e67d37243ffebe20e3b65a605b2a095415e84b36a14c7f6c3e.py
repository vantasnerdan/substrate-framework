"""Finite-time numerical evolution for the normalized 1+1 sine-Gordon PDE.

The exact framework convention is ``phi_tt-phi_xx+sin(phi)=0``.  This module
provides three deliberately separate numerical surfaces:

* a periodic centered spatial discretization with leapfrog and adaptive
  method-of-lines time integrators, useful for exact-solution convergence and
  independent-method checks;
* a homogeneous-Dirichlet interval with localized bulk forcing and declared
  damping, including source-work and dissipation ledgers;
* a half-line-style finite interval with prescribed left Neumann data and a
  first-order outgoing condition at the right endpoint, useful for auditing
  driven-boundary experiments.

Every result is finite-grid and finite-time evidence.  Endpoint field
differences are named as such because they become topological charges only
when the relevant endpoints satisfy the hypotheses of the accepted
topological-current claim.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .numerics import (
    NumericalFailure,
    SolverTolerances,
    solve_method_of_lines,
    trapezoid_integral,
)

FloatArray = NDArray[np.float64]
ScalarTrace = Callable[[float], float]
BulkSource = Callable[[float], ArrayLike]


@dataclass(frozen=True)
class PeriodicSineGordonEvolution:
    """Sampled periodic evolution and energy diagnostics."""

    coordinate: FloatArray
    time: FloatArray
    field: FloatArray
    velocity: FloatArray
    energy: FloatArray
    method: str
    spatial_step: float
    time_step: float | None
    function_evaluations: int | None


@dataclass(frozen=True)
class DrivenSineGordonEvolution:
    """Finite-time driven-Neumann evolution with explicit flux diagnostics."""

    coordinate: FloatArray
    final_field: FloatArray
    final_velocity: FloatArray
    boundary_time: FloatArray
    boundary_field: FloatArray
    boundary_velocity: FloatArray
    boundary_coordinate_derivative: FloatArray
    initial_energy: float
    final_energy: float
    cumulative_left_work: float
    cumulative_right_flux: float
    energy_balance_residual: float
    endpoint_field_difference: float
    endpoint_charge_coordinate: float
    bulk_start: float
    bulk_endpoint_field_difference: float
    bulk_endpoint_charge_coordinate: float
    boundary_sign_correlation: float
    spatial_step: float
    time_step: float | None
    method: str
    function_evaluations: int | None


@dataclass(frozen=True)
class BulkDrivenSineGordonEvolution:
    """Finite-time homogeneous-Dirichlet evolution with bulk forcing.

    ``cumulative_source_work`` approximates ``integral dt dx J*phi_t`` and
    ``cumulative_damping_loss`` approximates the positive integral
    ``integral dt dx gamma*phi_t**2``.  A claim must refine their time
    sampling when the balance residual is load bearing.
    """

    coordinate: FloatArray
    time: FloatArray
    center_field: FloatArray
    core_energy: FloatArray
    total_energy: FloatArray
    final_field: FloatArray
    final_velocity: FloatArray
    initial_energy: float
    cumulative_source_work: float
    cumulative_damping_loss: float
    energy_balance_residual: float
    core_radius: float
    spatial_step: float
    time_step: float | None
    method: str
    function_evaluations: int | None


@dataclass(frozen=True)
class BreatherTraceFit:
    """Constrained fit of a center trace to the exact rest-breather family."""

    angular_frequency: float
    phase_at_window_start: float
    relative_rms_error: float
    max_abs_error: float
    fitted_energy: float
    cycles_in_window: float


@dataclass(frozen=True)
class BreatherSnapshotFit:
    """Constrained phase-space fit to an exact centered rest breather."""

    angular_frequency: float
    phase: float
    joint_relative_l2_error: float
    field_relative_l2_error: float
    velocity_relative_l2_error: float
    gradient_relative_l2_error: float
    fitted_energy: float
    fit_radius: float


def moving_breather_samples(
    coordinate: ArrayLike,
    time: float,
    omega: float,
    velocity: float = 0.0,
    center: float = 0.0,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return exact moving-breather ``(phi,phi_t,phi_x)`` samples.

    The phase convention is the Lorentz transform of the accepted rest
    breather.  ``velocity`` is the center velocity, so a negative value moves
    the center toward decreasing coordinate.
    """

    x = _real_vector(coordinate, "coordinate")
    instant = _finite_scalar(time, "time")
    frequency = _open_unit_scalar(omega, "omega")
    speed = _open_signed_unit_scalar(velocity, "velocity")
    origin = _finite_scalar(center, "center")

    inverse_width = np.sqrt(1.0 - frequency**2)
    gamma = 1.0 / np.sqrt(1.0 - speed**2)
    xi = gamma * (x - origin - speed * instant)
    tau = gamma * (instant - speed * (x - origin))
    numerator = inverse_width * np.sin(frequency * tau)
    denominator = frequency * np.cosh(inverse_width * xi)
    norm = np.square(numerator) + np.square(denominator)
    if np.any(norm == 0.0) or not np.all(np.isfinite(norm)):
        raise NumericalFailure("moving-breather evaluation lost finite normalization")

    field = 4.0 * np.arctan2(numerator, denominator)
    dxi_dt = -gamma * speed
    dtau_dt = gamma
    dnum_dt = inverse_width * frequency * np.cos(frequency * tau) * dtau_dt
    dden_dt = frequency * inverse_width * np.sinh(inverse_width * xi) * dxi_dt
    velocity_samples = 4.0 * (
        denominator * dnum_dt - numerator * dden_dt
    ) / norm

    dxi_dx = gamma
    dtau_dx = -gamma * speed
    dnum_dx = inverse_width * frequency * np.cos(frequency * tau) * dtau_dx
    dden_dx = frequency * inverse_width * np.sinh(inverse_width * xi) * dxi_dx
    spatial_derivative = 4.0 * (
        denominator * dnum_dx - numerator * dden_dx
    ) / norm
    return (
        np.asarray(field, dtype=np.float64),
        np.asarray(velocity_samples, dtype=np.float64),
        np.asarray(spatial_derivative, dtype=np.float64),
    )


def endpoint_charge_coordinate(field_lower: float, field_upper: float) -> float:
    """Return ``(phi_upper-phi_lower)/(2*pi)`` without asserting vacuum data."""

    lower = _finite_scalar(field_lower, "field_lower")
    upper = _finite_scalar(field_upper, "field_upper")
    return float((upper - lower) / (2.0 * np.pi))


def sampled_boundary_sign_correlation(
    time: ArrayLike,
    boundary_velocity: ArrayLike,
    boundary_coordinate_derivative: ArrayLike,
) -> float:
    """Integrate ``sign(phi_t)*phi_x`` on a sampled boundary trace."""

    samples = _strict_time(time)
    velocity = _same_shape_vector(boundary_velocity, samples, "boundary_velocity")
    derivative = _same_shape_vector(
        boundary_coordinate_derivative,
        samples,
        "boundary_coordinate_derivative",
    )
    return trapezoid_integral(np.sign(velocity) * derivative, samples)


def periodic_sine_gordon_energy(
    field: ArrayLike,
    velocity: ArrayLike,
    spatial_step: float,
) -> float:
    """Return the centered-grid periodic Hamiltonian energy."""

    phi = _real_vector(field, "field")
    pi = _same_shape_vector(velocity, phi, "velocity")
    dx = _positive_scalar(spatial_step, "spatial_step")
    # The forward-edge difference is the discrete gradient whose variational
    # derivative is the centered Laplacian used by both periodic evolvers.
    # A centered nodal derivative would be a useful continuum estimator but
    # would not be the semi-discrete invariant of this operator.
    derivative = (np.roll(phi, -1) - phi) / dx
    density = (
        0.5 * np.square(pi)
        + 0.5 * np.square(derivative)
        + 1.0
        - np.cos(phi)
    )
    result = float(dx * np.sum(density))
    if not np.isfinite(result):
        raise NumericalFailure("periodic sine-Gordon energy became non-finite")
    return result


def homogeneous_dirichlet_sine_gordon_energy(
    field: ArrayLike,
    velocity: ArrayLike,
    spatial_step: float,
) -> float:
    """Return the centered-lattice homogeneous-Dirichlet Hamiltonian.

    The edge-gradient term is the potential whose variation gives the centered
    interior Laplacian used by the bulk-driven evolvers.  Fixed endpoint
    velocities leave no boundary-work term.
    """

    phi = _real_vector(field, "field")
    pi = _same_shape_vector(velocity, phi, "velocity")
    dx = _positive_scalar(spatial_step, "spatial_step")
    if phi.size < 5:
        raise ValueError("homogeneous-Dirichlet field requires at least five points")
    if not np.allclose(phi[[0, -1]], 0.0, atol=1.0e-12, rtol=0.0):
        raise ValueError("field endpoints must satisfy homogeneous Dirichlet data")
    if not np.allclose(pi[[0, -1]], 0.0, atol=1.0e-12, rtol=0.0):
        raise ValueError("velocity endpoints must satisfy fixed Dirichlet data")
    edge_energy = 0.5 * np.sum(np.square(np.diff(phi))) / dx
    nodal_energy = dx * np.sum(0.5 * np.square(pi) + 1.0 - np.cos(phi))
    result = float(edge_energy + nodal_energy)
    if not np.isfinite(result):
        raise NumericalFailure("homogeneous-Dirichlet energy became non-finite")
    return result


def fit_rest_breather_center_trace(
    time: ArrayLike,
    center_field: ArrayLike,
    *,
    frequency_bounds: tuple[float, float] = (0.05, 0.99),
) -> BreatherTraceFit:
    """Fit a trace to the exact C-SG-001 center-field waveform.

    The amplitude remains constrained by ``eta/omega`` and is not fitted
    independently.  This is one classifier diagnostic, not by itself proof
    that the full spatial state is a breather.
    """

    samples = _strict_time(time)
    trace = _same_shape_vector(center_field, samples, "center_field")
    lower = _open_unit_scalar(frequency_bounds[0], "frequency_bounds[0]")
    upper = _open_unit_scalar(frequency_bounds[1], "frequency_bounds[1]")
    if upper <= lower:
        raise ValueError("frequency bounds must be strictly increasing")
    centered = trace - np.mean(trace)
    scale = float(np.sqrt(np.mean(np.square(centered))))
    if scale == 0.0:
        raise ValueError("center_field has no oscillatory component")
    relative_time = samples - samples[0]

    def model(parameters: NDArray[np.float64]) -> FloatArray:
        frequency, phase = parameters
        inverse_width = np.sqrt(1.0 - frequency**2)
        return np.asarray(
            4.0
            * np.arctan(
                inverse_width
                / frequency
                * np.sin(frequency * relative_time + phase)
            ),
            dtype=np.float64,
        )

    def residual(parameters: NDArray[np.float64]) -> FloatArray:
        return np.asarray((model(parameters) - trace) / scale, dtype=np.float64)

    spacing = float(np.mean(np.diff(samples)))
    frequencies = 2.0 * np.pi * np.fft.rfftfreq(samples.size, spacing)
    spectrum = np.abs(np.fft.rfft(centered))
    eligible = np.flatnonzero((frequencies >= lower) & (frequencies <= upper))
    ranked = (
        eligible[np.argsort(spectrum[eligible])[-4:]]
        if eligible.size
        else np.array([], dtype=int)
    )
    initial_frequencies = [float(frequencies[index]) for index in ranked]
    initial_frequencies.extend(np.linspace(lower, upper, 7).tolist())
    best: tuple[float, NDArray[np.float64]] | None = None
    for frequency in initial_frequencies:
        clipped = min(max(frequency, lower + 1.0e-10), upper - 1.0e-10)
        for phase in (-np.pi, -0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi):
            result = least_squares(
                residual,
                np.array([clipped, phase], dtype=np.float64),
                bounds=(np.array([lower, -np.pi]), np.array([upper, np.pi])),
                max_nfev=2000,
            )
            score = float(np.mean(np.square(result.fun)))
            if best is None or score < best[0]:
                best = (score, np.asarray(result.x, dtype=np.float64))
    if best is None:
        raise NumericalFailure("breather center-trace fit found no admissible start")
    frequency, phase = (float(value) for value in best[1])
    fitted = model(best[1])
    difference = fitted - trace
    return BreatherTraceFit(
        angular_frequency=frequency,
        phase_at_window_start=phase,
        relative_rms_error=float(np.sqrt(np.mean(np.square(difference))) / scale),
        max_abs_error=float(np.max(np.abs(difference))),
        fitted_energy=float(16.0 * np.sqrt(1.0 - frequency**2)),
        cycles_in_window=float(
            frequency * (samples[-1] - samples[0]) / (2.0 * np.pi)
        ),
    )


def fit_rest_breather_snapshot(
    coordinate: ArrayLike,
    field: ArrayLike,
    velocity: ArrayLike,
    *,
    fit_radius: float,
    frequency_bounds: tuple[float, float] = (0.05, 0.99),
) -> BreatherSnapshotFit:
    """Fit one phase-space snapshot to the exact C-SG-001 family.

    The breather remains centered at the declared source origin, and its
    spatial width, amplitude, velocity, and energy are all fixed by one
    frequency and one phase.  No independent amplitude or width parameter is
    available to absorb a nonbreather state.
    """

    x, _dx = _uniform_inclusive_coordinate(coordinate)
    phi = _same_shape_vector(field, x, "field")
    pi = _same_shape_vector(velocity, x, "velocity")
    radius = _positive_scalar(fit_radius, "fit_radius")
    if not x[0] < -radius < radius < x[-1]:
        raise ValueError("fit_radius must lie strictly inside a domain spanning zero")
    mask = np.abs(x) < radius
    if np.count_nonzero(mask) < 5:
        raise ValueError("fit_radius must contain at least five coordinate points")
    x_fit = x[mask]
    phi_fit = phi[mask]
    pi_fit = pi[mask]
    lower = _open_unit_scalar(frequency_bounds[0], "frequency_bounds[0]")
    upper = _open_unit_scalar(frequency_bounds[1], "frequency_bounds[1]")
    if upper <= lower:
        raise ValueError("frequency bounds must be strictly increasing")
    scale = float(np.sqrt(np.mean(np.square(phi_fit) + np.square(pi_fit))))
    if scale == 0.0:
        raise ValueError("snapshot has no phase-space signal inside fit_radius")

    def model(
        parameters: NDArray[np.float64],
    ) -> tuple[FloatArray, FloatArray, FloatArray]:
        frequency, phase = parameters
        return moving_breather_samples(
            x_fit,
            float(phase / frequency),
            float(frequency),
        )

    def residual(parameters: NDArray[np.float64]) -> FloatArray:
        model_field, model_velocity, _model_gradient = model(parameters)
        return np.concatenate(
            (model_field - phi_fit, model_velocity - pi_fit)
        ) / scale

    best: tuple[float, NDArray[np.float64]] | None = None
    for frequency in np.linspace(lower, upper, 9):
        for phase in (-np.pi, -0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi):
            result = least_squares(
                residual,
                np.array([frequency, phase], dtype=np.float64),
                bounds=(np.array([lower, -np.pi]), np.array([upper, np.pi])),
                max_nfev=2000,
            )
            score = float(np.mean(np.square(result.fun)))
            if best is None or score < best[0]:
                best = (score, np.asarray(result.x, dtype=np.float64))
    if best is None:
        raise NumericalFailure("breather snapshot fit found no admissible start")
    frequency, phase = (float(value) for value in best[1])
    model_field, model_velocity, model_gradient = model(best[1])
    measured_gradient = np.gradient(phi, x)[mask]

    def relative_l2(difference: FloatArray, reference: FloatArray) -> float:
        numerator = trapezoid_integral(np.square(difference), x_fit)
        denominator = trapezoid_integral(np.square(reference), x_fit)
        if denominator <= np.finfo(np.float64).tiny:
            return float(np.sqrt(numerator / np.finfo(np.float64).tiny))
        return float(np.sqrt(numerator / denominator))

    joint_numerator = trapezoid_integral(
        np.square(model_field - phi_fit) + np.square(model_velocity - pi_fit),
        x_fit,
    )
    joint_denominator = trapezoid_integral(
        np.square(phi_fit) + np.square(pi_fit),
        x_fit,
    )
    return BreatherSnapshotFit(
        angular_frequency=frequency,
        phase=phase,
        joint_relative_l2_error=float(np.sqrt(joint_numerator / joint_denominator)),
        field_relative_l2_error=relative_l2(model_field - phi_fit, phi_fit),
        velocity_relative_l2_error=relative_l2(model_velocity - pi_fit, pi_fit),
        gradient_relative_l2_error=relative_l2(
            model_gradient - measured_gradient,
            measured_gradient,
        ),
        fitted_energy=float(16.0 * np.sqrt(1.0 - frequency**2)),
        fit_radius=radius,
    )


def evolve_periodic_sine_gordon_leapfrog(
    coordinate: ArrayLike,
    initial_field: ArrayLike,
    initial_velocity: ArrayLike,
    final_time: float,
    requested_time_step: float,
    *,
    sample_stride: int = 1,
) -> PeriodicSineGordonEvolution:
    """Evolve the periodic centered semi-discretization by velocity Verlet."""

    x, dx = _uniform_periodic_coordinate(coordinate)
    phi = _same_shape_vector(initial_field, x, "initial_field").copy()
    velocity = _same_shape_vector(initial_velocity, x, "initial_velocity").copy()
    duration = _positive_scalar(final_time, "final_time")
    proposed_step = _positive_scalar(requested_time_step, "requested_time_step")
    stride = _positive_integer(sample_stride, "sample_stride")
    steps = max(1, int(np.ceil(duration / proposed_step)))
    dt = duration / steps
    if dt / dx > 1.0:
        raise ValueError("periodic leapfrog requires time_step/spatial_step <= 1")

    sample_steps = set(range(0, steps + 1, stride))
    sample_steps.add(steps)
    times: list[float] = [0.0]
    fields: list[FloatArray] = [phi.copy()]
    velocities: list[FloatArray] = [velocity.copy()]
    energies: list[float] = [periodic_sine_gordon_energy(phi, velocity, dx)]

    acceleration = _periodic_acceleration(phi, dx)
    half_velocity = velocity + 0.5 * dt * acceleration
    for step in range(1, steps + 1):
        phi = phi + dt * half_velocity
        acceleration = _periodic_acceleration(phi, dx)
        velocity = half_velocity + 0.5 * dt * acceleration
        half_velocity = half_velocity + dt * acceleration
        if not np.all(np.isfinite(phi)) or not np.all(np.isfinite(velocity)):
            raise NumericalFailure(f"periodic leapfrog failed at step {step}")
        if step in sample_steps:
            times.append(step * dt)
            fields.append(phi.copy())
            velocities.append(velocity.copy())
            energies.append(periodic_sine_gordon_energy(phi, velocity, dx))

    return PeriodicSineGordonEvolution(
        coordinate=x,
        time=np.asarray(times, dtype=np.float64),
        field=np.asarray(fields, dtype=np.float64),
        velocity=np.asarray(velocities, dtype=np.float64),
        energy=np.asarray(energies, dtype=np.float64),
        method="periodic-centered-velocity-verlet",
        spatial_step=dx,
        time_step=dt,
        function_evaluations=None,
    )


def evolve_periodic_sine_gordon_mol(
    coordinate: ArrayLike,
    initial_field: ArrayLike,
    initial_velocity: ArrayLike,
    sample_times: ArrayLike,
    *,
    tolerances: SolverTolerances = SolverTolerances(),
    method: str = "DOP853",
) -> PeriodicSineGordonEvolution:
    """Evolve the same periodic spatial operator with an adaptive IVP method."""

    x, dx = _uniform_periodic_coordinate(coordinate)
    phi0 = _same_shape_vector(initial_field, x, "initial_field")
    velocity0 = _same_shape_vector(initial_velocity, x, "initial_velocity")
    times = _strict_time(sample_times)
    if times[0] != 0.0:
        raise ValueError("sample_times must begin at zero")
    point_count = x.size

    def rhs(_time: float, state: FloatArray) -> FloatArray:
        phi = state[:point_count]
        velocity = state[point_count:]
        return np.concatenate((velocity, _periodic_acceleration(phi, dx)))

    evidence = solve_method_of_lines(
        rhs,
        (float(times[0]), float(times[-1])),
        np.concatenate((phi0, velocity0)),
        sample_times=times,
        tolerances=tolerances,
        method=method,
    )
    fields = evidence.state[:point_count].T
    velocities = evidence.state[point_count:].T
    energies = np.asarray(
        [
            periodic_sine_gordon_energy(field, velocity, dx)
            for field, velocity in zip(fields, velocities)
        ],
        dtype=np.float64,
    )
    return PeriodicSineGordonEvolution(
        coordinate=x,
        time=evidence.time,
        field=fields,
        velocity=velocities,
        energy=energies,
        method=f"periodic-centered-{method}",
        spatial_step=dx,
        time_step=None,
        function_evaluations=evidence.function_evaluations,
    )


def evolve_bulk_driven_sine_gordon_leapfrog(
    coordinate: ArrayLike,
    initial_field: ArrayLike,
    initial_velocity: ArrayLike,
    bulk_source: BulkSource,
    damping: ArrayLike,
    final_time: float,
    requested_time_step: float,
    *,
    core_radius: float,
    sample_interval: float,
) -> BulkDrivenSineGordonEvolution:
    """Evolve a bulk-forced homogeneous-Dirichlet problem by leapfrog.

    The semi-discrete equation is
    ``phi_tt=Delta_h(phi)-sin(phi)+J(t)-gamma*phi_t``. Damping is centered in
    time. Source-work and damping-loss integrals use every time step rather
    than the coarser diagnostic interval.
    """

    x, dx = _uniform_inclusive_coordinate(coordinate)
    phi_previous = _same_shape_vector(initial_field, x, "initial_field").copy()
    velocity0 = _same_shape_vector(initial_velocity, x, "initial_velocity").copy()
    _require_homogeneous_dirichlet_data(phi_previous, velocity0)
    gamma = _nonnegative_shape_vector(damping, x, "damping")
    duration = _positive_scalar(final_time, "final_time")
    proposed_step = _positive_scalar(requested_time_step, "requested_time_step")
    radius = _positive_scalar(core_radius, "core_radius")
    sample = _positive_scalar(sample_interval, "sample_interval")
    if not x[0] < -radius < radius < x[-1]:
        raise ValueError("core_radius must lie strictly inside a domain spanning zero")
    steps = max(3, int(np.ceil(duration / proposed_step)))
    dt = duration / steps
    if dt / dx > 1.0:
        raise ValueError("bulk-driven leapfrog requires time_step/spatial_step <= 1")
    record_stride = max(1, int(round(sample / dt)))
    core_mask = np.abs(x) < radius
    center_index = int(np.argmin(np.abs(x)))

    source0 = _bulk_source_samples(bulk_source, 0.0, x)
    acceleration0 = _dirichlet_bulk_acceleration(phi_previous, dx, source0)
    phi_current = (
        phi_previous
        + dt * velocity0
        + 0.5 * dt**2 * (acceleration0 - gamma * velocity0)
    )
    phi_current[[0, -1]] = 0.0
    _require_finite(phi_current, "first bulk-driven leapfrog step")

    initial_energy = homogeneous_dirichlet_sine_gordon_energy(
        phi_previous, velocity0, dx
    )
    source_power_previous = float(dx * np.sum(source0[1:-1] * velocity0[1:-1]))
    damping_power_previous = float(
        dx * np.sum(gamma[1:-1] * np.square(velocity0[1:-1]))
    )
    source_work = 0.0
    damping_loss = 0.0
    times: list[float] = [0.0]
    centers: list[float] = [float(phi_previous[center_index])]
    core_energies: list[float] = [
        _sampled_core_energy(phi_previous, velocity0, x, core_mask)
    ]
    total_energies: list[float] = [initial_energy]

    for step in range(1, steps):
        time = step * dt
        source = _bulk_source_samples(bulk_source, time, x)
        acceleration = _dirichlet_bulk_acceleration(phi_current, dx, source)
        phi_next = (
            2.0 * phi_current
            - (1.0 - 0.5 * gamma * dt) * phi_previous
            + dt**2 * acceleration
        ) / (1.0 + 0.5 * gamma * dt)
        phi_next[[0, -1]] = 0.0
        _require_finite(phi_next, f"bulk-driven leapfrog step {step + 1}")
        velocity_current = (phi_next - phi_previous) / (2.0 * dt)
        _require_finite(velocity_current, f"bulk-driven velocity step {step}")
        source_power = float(dx * np.sum(source[1:-1] * velocity_current[1:-1]))
        damping_power = float(
            dx * np.sum(gamma[1:-1] * np.square(velocity_current[1:-1]))
        )
        source_work += 0.5 * dt * (source_power_previous + source_power)
        damping_loss += 0.5 * dt * (damping_power_previous + damping_power)
        source_power_previous = source_power
        damping_power_previous = damping_power
        if step % record_stride == 0:
            times.append(time)
            centers.append(float(phi_current[center_index]))
            core_energies.append(
                _sampled_core_energy(phi_current, velocity_current, x, core_mask)
            )
            total_energies.append(
                homogeneous_dirichlet_sine_gordon_energy(
                    phi_current, velocity_current, dx
                )
            )
        phi_previous, phi_current = phi_current, phi_next

    # Recover the velocity at ``duration`` from the backward field difference
    # and the endpoint acceleration.  The raw backward difference is only
    # first-order accurate and would manufacture an O(dt) energy-balance
    # residual even though the centered recurrence is second order.
    backward_velocity = (phi_current - phi_previous) / dt
    source_final = _bulk_source_samples(bulk_source, duration, x)
    final_acceleration = _dirichlet_bulk_acceleration(
        phi_current,
        dx,
        source_final,
    )
    final_velocity = (
        backward_velocity + 0.5 * dt * final_acceleration
    ) / (1.0 + 0.5 * gamma * dt)
    final_velocity[[0, -1]] = 0.0
    source_power_final = float(
        dx * np.sum(source_final[1:-1] * final_velocity[1:-1])
    )
    damping_power_final = float(
        dx * np.sum(gamma[1:-1] * np.square(final_velocity[1:-1]))
    )
    source_work += 0.5 * dt * (source_power_previous + source_power_final)
    damping_loss += 0.5 * dt * (damping_power_previous + damping_power_final)
    final_energy = homogeneous_dirichlet_sine_gordon_energy(
        phi_current, final_velocity, dx
    )
    if times[-1] != duration:
        times.append(duration)
        centers.append(float(phi_current[center_index]))
        core_energies.append(
            _sampled_core_energy(phi_current, final_velocity, x, core_mask)
        )
        total_energies.append(final_energy)
    return BulkDrivenSineGordonEvolution(
        coordinate=x,
        time=np.asarray(times, dtype=np.float64),
        center_field=np.asarray(centers, dtype=np.float64),
        core_energy=np.asarray(core_energies, dtype=np.float64),
        total_energy=np.asarray(total_energies, dtype=np.float64),
        final_field=np.asarray(phi_current, dtype=np.float64),
        final_velocity=np.asarray(final_velocity, dtype=np.float64),
        initial_energy=initial_energy,
        cumulative_source_work=float(source_work),
        cumulative_damping_loss=float(damping_loss),
        energy_balance_residual=float(
            final_energy - initial_energy - source_work + damping_loss
        ),
        core_radius=radius,
        spatial_step=dx,
        time_step=dt,
        method="bulk-dirichlet-centered-leapfrog",
        function_evaluations=None,
    )


def evolve_bulk_driven_sine_gordon_mol(
    coordinate: ArrayLike,
    initial_field: ArrayLike,
    initial_velocity: ArrayLike,
    bulk_source: BulkSource,
    damping: ArrayLike,
    sample_times: ArrayLike,
    *,
    core_radius: float,
    tolerances: SolverTolerances = SolverTolerances(),
    method: str = "DOP853",
) -> BulkDrivenSineGordonEvolution:
    """Integrate the same bulk-driven spatial operator by adaptive MOL."""

    x, dx = _uniform_inclusive_coordinate(coordinate)
    phi0 = _same_shape_vector(initial_field, x, "initial_field")
    velocity0 = _same_shape_vector(initial_velocity, x, "initial_velocity")
    _require_homogeneous_dirichlet_data(phi0, velocity0)
    gamma = _nonnegative_shape_vector(damping, x, "damping")
    times = _strict_time(sample_times)
    if times[0] != 0.0:
        raise ValueError("sample_times must begin at zero")
    radius = _positive_scalar(core_radius, "core_radius")
    if not x[0] < -radius < radius < x[-1]:
        raise ValueError("core_radius must lie strictly inside a domain spanning zero")
    point_count = x.size
    core_mask = np.abs(x) < radius

    def rhs(time: float, state: FloatArray) -> FloatArray:
        phi = state[:point_count]
        velocity = state[point_count:]
        source = _bulk_source_samples(bulk_source, time, x)
        acceleration = (
            _dirichlet_bulk_acceleration(phi, dx, source) - gamma * velocity
        )
        derivative = np.concatenate((velocity, acceleration))
        derivative[[0, point_count - 1, point_count, 2 * point_count - 1]] = 0.0
        return derivative

    evidence = solve_method_of_lines(
        rhs,
        (float(times[0]), float(times[-1])),
        np.concatenate((phi0, velocity0)),
        sample_times=times,
        tolerances=tolerances,
        method=method,
    )
    fields = evidence.state[:point_count]
    velocities = evidence.state[point_count:]
    source_powers: list[float] = []
    damping_powers: list[float] = []
    core_energies: list[float] = []
    total_energies: list[float] = []
    for index, time in enumerate(times):
        source = _bulk_source_samples(bulk_source, float(time), x)
        velocity = velocities[:, index]
        source_powers.append(float(dx * np.sum(source[1:-1] * velocity[1:-1])))
        damping_powers.append(
            float(dx * np.sum(gamma[1:-1] * np.square(velocity[1:-1])))
        )
        core_energies.append(
            _sampled_core_energy(fields[:, index], velocity, x, core_mask)
        )
        total_energies.append(
            homogeneous_dirichlet_sine_gordon_energy(
                fields[:, index], velocity, dx
            )
        )
    source_work = trapezoid_integral(source_powers, times)
    damping_loss = trapezoid_integral(damping_powers, times)
    initial_energy = total_energies[0]
    final_energy = total_energies[-1]
    return BulkDrivenSineGordonEvolution(
        coordinate=x,
        time=np.asarray(times, dtype=np.float64),
        center_field=np.asarray(fields[np.argmin(np.abs(x))], dtype=np.float64),
        core_energy=np.asarray(core_energies, dtype=np.float64),
        total_energy=np.asarray(total_energies, dtype=np.float64),
        final_field=np.asarray(fields[:, -1], dtype=np.float64),
        final_velocity=np.asarray(velocities[:, -1], dtype=np.float64),
        initial_energy=float(initial_energy),
        cumulative_source_work=float(source_work),
        cumulative_damping_loss=float(damping_loss),
        energy_balance_residual=float(
            final_energy - initial_energy - source_work + damping_loss
        ),
        core_radius=radius,
        spatial_step=dx,
        time_step=None,
        method=f"bulk-dirichlet-{method}",
        function_evaluations=evidence.function_evaluations,
    )


def evolve_driven_sine_gordon_leapfrog(
    coordinate: ArrayLike,
    initial_field: ArrayLike,
    initial_velocity: ArrayLike,
    left_neumann: ScalarTrace,
    final_time: float,
    requested_time_step: float,
    *,
    bulk_start: float,
) -> DrivenSineGordonEvolution:
    """Evolve a left-Neumann/right-outgoing finite-interval experiment.

    The grid includes both endpoints and its measured spacing is used by every
    stencil.  The left ghost value enforces ``phi_x=left_neumann(t)`` to second
    order.  The right endpoint uses the declared first-order Sommerfeld update
    ``phi_t+phi_x=0``; no extra cells are overwritten.
    """

    x, dx = _uniform_inclusive_coordinate(coordinate)
    phi_previous = _same_shape_vector(initial_field, x, "initial_field").copy()
    velocity0 = _same_shape_vector(initial_velocity, x, "initial_velocity").copy()
    duration = _positive_scalar(final_time, "final_time")
    proposed_step = _positive_scalar(requested_time_step, "requested_time_step")
    lower_bulk = _finite_scalar(bulk_start, "bulk_start")
    if not x[0] < lower_bulk < x[-1]:
        raise ValueError("bulk_start must lie strictly inside the coordinate domain")
    steps = max(3, int(np.ceil(duration / proposed_step)))
    dt = duration / steps
    courant = dt / dx
    if courant > 1.0:
        raise ValueError("driven leapfrog requires time_step/spatial_step <= 1")

    drive0 = _finite_scalar(left_neumann(0.0), "left_neumann(0)")
    acceleration0 = _neumann_left_acceleration(phi_previous, dx, drive0)
    phi_current = (
        phi_previous + dt * velocity0 + 0.5 * dt**2 * acceleration0
    )
    phi_current[-1] = phi_previous[-1] - courant * (
        phi_previous[-1] - phi_previous[-2]
    )
    _require_finite(phi_current, "first driven leapfrog step")

    initial_energy = _interval_energy(phi_previous, velocity0, x, drive0)
    trace_time: list[float] = []
    trace_field: list[float] = []
    trace_velocity: list[float] = []
    trace_derivative: list[float] = []
    left_power: list[float] = []
    right_power: list[float] = []

    for step in range(1, steps):
        time = step * dt
        drive = _finite_scalar(left_neumann(time), "left_neumann(time)")
        acceleration = _neumann_left_acceleration(phi_current, dx, drive)
        phi_next = 2.0 * phi_current - phi_previous + dt**2 * acceleration
        phi_next[-1] = phi_current[-1] - courant * (
            phi_current[-1] - phi_current[-2]
        )
        velocity_current = (phi_next - phi_previous) / (2.0 * dt)
        right_derivative = (phi_current[-1] - phi_current[-2]) / dx
        _require_finite(phi_next, f"driven leapfrog step {step + 1}")
        _require_finite(velocity_current, f"driven velocity step {step}")

        trace_time.append(time)
        trace_field.append(float(phi_current[0]))
        trace_velocity.append(float(velocity_current[0]))
        trace_derivative.append(drive)
        left_power.append(float(-velocity_current[0] * drive))
        right_power.append(float(velocity_current[-1] * right_derivative))
        phi_previous = phi_current
        phi_current = phi_next

    final_velocity = (phi_current - phi_previous) / dt
    final_drive = _finite_scalar(left_neumann(duration), "left_neumann(final_time)")
    final_energy = _interval_energy(phi_current, final_velocity, x, final_drive)
    times = np.asarray(trace_time, dtype=np.float64)
    velocities = np.asarray(trace_velocity, dtype=np.float64)
    derivatives = np.asarray(trace_derivative, dtype=np.float64)
    left_work = trapezoid_integral(left_power, times)
    outgoing_flux = trapezoid_integral(right_power, times)
    energy_residual = final_energy - initial_energy - left_work - outgoing_flux

    bulk_index = int(np.searchsorted(x, lower_bulk, side="left"))
    endpoint_difference = float(phi_current[-1] - phi_current[0])
    bulk_difference = float(phi_current[-1] - phi_current[bulk_index])
    return DrivenSineGordonEvolution(
        coordinate=x,
        final_field=np.asarray(phi_current, dtype=np.float64),
        final_velocity=np.asarray(final_velocity, dtype=np.float64),
        boundary_time=times,
        boundary_field=np.asarray(trace_field, dtype=np.float64),
        boundary_velocity=velocities,
        boundary_coordinate_derivative=derivatives,
        initial_energy=initial_energy,
        final_energy=final_energy,
        cumulative_left_work=left_work,
        cumulative_right_flux=outgoing_flux,
        energy_balance_residual=float(energy_residual),
        endpoint_field_difference=endpoint_difference,
        endpoint_charge_coordinate=endpoint_charge_coordinate(
            phi_current[0], phi_current[-1]
        ),
        bulk_start=float(x[bulk_index]),
        bulk_endpoint_field_difference=bulk_difference,
        bulk_endpoint_charge_coordinate=endpoint_charge_coordinate(
            phi_current[bulk_index], phi_current[-1]
        ),
        boundary_sign_correlation=sampled_boundary_sign_correlation(
            times,
            velocities,
            derivatives,
        ),
        spatial_step=dx,
        time_step=dt,
        method="left-neumann-right-sommerfeld-leapfrog",
        function_evaluations=None,
    )


def evolve_driven_sine_gordon_mol(
    coordinate: ArrayLike,
    initial_field: ArrayLike,
    initial_velocity: ArrayLike,
    left_neumann: ScalarTrace,
    sample_times: ArrayLike,
    *,
    bulk_start: float,
    tolerances: SolverTolerances = SolverTolerances(),
    method: str = "DOP853",
) -> DrivenSineGordonEvolution:
    """Evolve the driven problem with an adaptive method-of-lines route.

    The centered interior and left ghost-point operators match the leapfrog
    route.  The right Sommerfeld condition is evolved as
    ``pi_t=-pi_x`` at the boundary, a distinct time-integration formulation
    suitable for a method cross-check.  Flux quadrature uses the caller's
    declared sample times, which must therefore be refined independently when
    an energy-balance number is load-bearing.
    """

    x, dx = _uniform_inclusive_coordinate(coordinate)
    phi0 = _same_shape_vector(initial_field, x, "initial_field")
    velocity0 = _same_shape_vector(initial_velocity, x, "initial_velocity")
    times = _strict_time(sample_times)
    if times[0] != 0.0:
        raise ValueError("sample_times must begin at zero")
    lower_bulk = _finite_scalar(bulk_start, "bulk_start")
    if not x[0] < lower_bulk < x[-1]:
        raise ValueError("bulk_start must lie strictly inside the coordinate domain")
    point_count = x.size

    def rhs(time: float, state: FloatArray) -> FloatArray:
        phi = state[:point_count]
        velocity = state[point_count:]
        drive = _finite_scalar(left_neumann(time), "left_neumann(time)")
        acceleration = _neumann_left_acceleration(phi, dx, drive)
        acceleration[-1] = -(velocity[-1] - velocity[-2]) / dx
        return np.concatenate((velocity, acceleration))

    evidence = solve_method_of_lines(
        rhs,
        (float(times[0]), float(times[-1])),
        np.concatenate((phi0, velocity0)),
        sample_times=times,
        tolerances=tolerances,
        method=method,
    )
    fields = evidence.state[:point_count]
    velocities = evidence.state[point_count:]
    final_field = fields[:, -1]
    final_velocity = velocities[:, -1]
    derivative_trace = np.asarray([left_neumann(time) for time in times], dtype=np.float64)
    _require_finite(derivative_trace, "left Neumann trace")
    left_power = -velocities[0] * derivative_trace
    right_derivative = (fields[-1] - fields[-2]) / dx
    right_power = velocities[-1] * right_derivative
    left_work = trapezoid_integral(left_power, times)
    outgoing_flux = trapezoid_integral(right_power, times)
    initial_energy = _interval_energy(phi0, velocity0, x, derivative_trace[0])
    final_energy = _interval_energy(
        final_field,
        final_velocity,
        x,
        derivative_trace[-1],
    )
    bulk_index = int(np.searchsorted(x, lower_bulk, side="left"))
    endpoint_difference = float(final_field[-1] - final_field[0])
    bulk_difference = float(final_field[-1] - final_field[bulk_index])
    return DrivenSineGordonEvolution(
        coordinate=x,
        final_field=np.asarray(final_field, dtype=np.float64),
        final_velocity=np.asarray(final_velocity, dtype=np.float64),
        boundary_time=times,
        boundary_field=np.asarray(fields[0], dtype=np.float64),
        boundary_velocity=np.asarray(velocities[0], dtype=np.float64),
        boundary_coordinate_derivative=derivative_trace,
        initial_energy=initial_energy,
        final_energy=final_energy,
        cumulative_left_work=left_work,
        cumulative_right_flux=outgoing_flux,
        energy_balance_residual=float(
            final_energy - initial_energy - left_work - outgoing_flux
        ),
        endpoint_field_difference=endpoint_difference,
        endpoint_charge_coordinate=endpoint_charge_coordinate(
            final_field[0], final_field[-1]
        ),
        bulk_start=float(x[bulk_index]),
        bulk_endpoint_field_difference=bulk_difference,
        bulk_endpoint_charge_coordinate=endpoint_charge_coordinate(
            final_field[bulk_index], final_field[-1]
        ),
        boundary_sign_correlation=sampled_boundary_sign_correlation(
            times,
            velocities[0],
            derivative_trace,
        ),
        spatial_step=dx,
        time_step=None,
        method=f"left-neumann-right-sommerfeld-{method}",
        function_evaluations=evidence.function_evaluations,
    )


def gaussian_sine_neumann_drive(
    amplitude: float,
    angular_frequency: float,
    center_time: float,
    width: float,
    phase: float,
) -> ScalarTrace:
    """Return a declared Gaussian-enveloped sinusoidal Neumann trace."""

    scale = _finite_scalar(amplitude, "amplitude")
    frequency = _positive_scalar(angular_frequency, "angular_frequency")
    center = _finite_scalar(center_time, "center_time")
    sigma = _positive_scalar(width, "width")
    phase_value = _finite_scalar(phase, "phase")

    def drive(time: float) -> float:
        instant = _finite_scalar(time, "time")
        envelope = np.exp(-0.5 * ((instant - center) / sigma) ** 2)
        return float(
            scale * envelope * np.sin(frequency * (instant - center) + phase_value)
        )

    return drive


def gaussian_sine_trace(
    amplitude: float,
    angular_frequency: float,
    center_time: float,
    width: float,
    phase: float = 0.0,
) -> ScalarTrace:
    """Return a declared Gaussian-enveloped sinusoidal scalar trace."""

    return gaussian_sine_neumann_drive(
        amplitude,
        angular_frequency,
        center_time,
        width,
        phase,
    )


def localized_sech_bulk_source(
    coordinate: ArrayLike,
    inverse_width: float,
    temporal_trace: ScalarTrace,
) -> BulkSource:
    """Return ``J(t,x)=sech(inverse_width*x)*temporal_trace(t)``."""

    x = _real_vector(coordinate, "coordinate")
    width = _positive_scalar(inverse_width, "inverse_width")
    profile = 1.0 / np.cosh(width * x)

    def source(time: float) -> FloatArray:
        value = _finite_scalar(temporal_trace(time), "temporal_trace(time)")
        return np.asarray(value * profile, dtype=np.float64)

    return source


def gaussian_sine_full_line_l2(
    amplitude: float,
    angular_frequency: float,
    width: float,
    phase: float = 0.0,
) -> float:
    """Return the exact full-line integral of a Gaussian-sine trace squared."""

    scale = _finite_scalar(amplitude, "amplitude")
    frequency = _positive_scalar(angular_frequency, "angular_frequency")
    sigma = _positive_scalar(width, "width")
    phase_value = _finite_scalar(phase, "phase")
    return float(
        scale**2
        * np.sqrt(np.pi)
        * sigma
        / 2.0
        * (1.0 - np.exp(-(frequency * sigma) ** 2) * np.cos(2.0 * phase_value))
    )


def _dirichlet_bulk_acceleration(
    field: FloatArray,
    spatial_step: float,
    source: FloatArray,
) -> FloatArray:
    acceleration = np.zeros_like(field)
    acceleration[1:-1] = (
        (field[2:] - 2.0 * field[1:-1] + field[:-2]) / spatial_step**2
        - np.sin(field[1:-1])
        + source[1:-1]
    )
    return acceleration


def _bulk_source_samples(
    source: BulkSource,
    time: float,
    coordinate: FloatArray,
) -> FloatArray:
    instant = _finite_scalar(time, "source time")
    values = np.asarray(source(instant), dtype=np.float64)
    if values.shape != coordinate.shape:
        raise ValueError("bulk_source must return one value per coordinate point")
    _require_finite(values, "bulk_source")
    return values


def _nonnegative_shape_vector(
    values: ArrayLike,
    reference: FloatArray,
    name: str,
) -> FloatArray:
    result = _same_shape_vector(values, reference, name)
    if np.any(result < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    return result


def _require_homogeneous_dirichlet_data(
    field: FloatArray,
    velocity: FloatArray,
) -> None:
    if not np.allclose(field[[0, -1]], 0.0, atol=1.0e-12, rtol=0.0):
        raise ValueError("initial_field endpoints must be zero")
    if not np.allclose(velocity[[0, -1]], 0.0, atol=1.0e-12, rtol=0.0):
        raise ValueError("initial_velocity endpoints must be zero")


def _sampled_core_energy(
    field: FloatArray,
    velocity: FloatArray,
    coordinate: FloatArray,
    core_mask: NDArray[np.bool_],
) -> float:
    derivative = np.gradient(field, coordinate)
    density = (
        0.5 * np.square(velocity)
        + 0.5 * np.square(derivative)
        + 1.0
        - np.cos(field)
    )
    return trapezoid_integral(density[core_mask], coordinate[core_mask])


def _periodic_acceleration(field: FloatArray, spatial_step: float) -> FloatArray:
    laplacian = (
        np.roll(field, -1) - 2.0 * field + np.roll(field, 1)
    ) / spatial_step**2
    return np.asarray(laplacian - np.sin(field), dtype=np.float64)


def _neumann_left_acceleration(
    field: FloatArray,
    spatial_step: float,
    left_derivative: float,
) -> FloatArray:
    acceleration = np.empty_like(field)
    acceleration[1:-1] = (
        field[2:] - 2.0 * field[1:-1] + field[:-2]
    ) / spatial_step**2 - np.sin(field[1:-1])
    acceleration[0] = (
        2.0 * (field[1] - field[0]) / spatial_step**2
        - 2.0 * left_derivative / spatial_step
        - np.sin(field[0])
    )
    acceleration[-1] = 0.0
    return acceleration


def _interval_energy(
    field: FloatArray,
    velocity: FloatArray,
    coordinate: FloatArray,
    left_derivative: float,
) -> float:
    spatial_derivative = np.empty_like(field)
    spatial_derivative[0] = left_derivative
    spatial_derivative[1:-1] = (
        field[2:] - field[:-2]
    ) / (2.0 * (coordinate[1] - coordinate[0]))
    spatial_derivative[-1] = (
        field[-1] - field[-2]
    ) / (coordinate[1] - coordinate[0])
    density = (
        0.5 * np.square(velocity)
        + 0.5 * np.square(spatial_derivative)
        + 1.0
        - np.cos(field)
    )
    return trapezoid_integral(density, coordinate)


def _uniform_periodic_coordinate(coordinate: ArrayLike) -> tuple[FloatArray, float]:
    x = _real_vector(coordinate, "coordinate")
    if x.size < 5:
        raise ValueError("periodic coordinate requires at least five points")
    return x, _uniform_spacing(x)


def _uniform_inclusive_coordinate(coordinate: ArrayLike) -> tuple[FloatArray, float]:
    x = _real_vector(coordinate, "coordinate")
    if x.size < 5:
        raise ValueError("coordinate requires at least five points")
    return x, _uniform_spacing(x)


def _uniform_spacing(coordinate: FloatArray) -> float:
    differences = np.diff(coordinate)
    if np.any(differences <= 0.0):
        raise ValueError("coordinate must be strictly increasing")
    spacing = float(np.mean(differences))
    if not np.allclose(differences, spacing, rtol=1.0e-11, atol=1.0e-13):
        raise ValueError("coordinate must be uniformly spaced")
    return spacing


def _strict_time(values: ArrayLike) -> FloatArray:
    time = _real_vector(values, "time")
    if time.size < 2 or np.any(np.diff(time) <= 0.0):
        raise ValueError("time must contain at least two increasing samples")
    return time


def _same_shape_vector(values: ArrayLike, reference: FloatArray, name: str) -> FloatArray:
    vector = _real_vector(values, name)
    if vector.shape != reference.shape:
        raise ValueError(f"{name} must match the reference shape")
    return vector


def _real_vector(values: ArrayLike, name: str) -> FloatArray:
    vector = np.asarray(values, dtype=np.float64)
    if vector.ndim != 1 or vector.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional array")
    _require_finite(vector, name)
    return vector


def _require_finite(values: FloatArray, name: str) -> None:
    if not np.all(np.isfinite(values)):
        raise NumericalFailure(f"{name} contains non-finite values")


def _finite_scalar(value: float, name: str) -> float:
    scalar = float(value)
    if not np.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    return scalar


def _positive_scalar(value: float, name: str) -> float:
    scalar = _finite_scalar(value, name)
    if scalar <= 0.0:
        raise ValueError(f"{name} must be positive")
    return scalar


def _open_unit_scalar(value: float, name: str) -> float:
    scalar = _finite_scalar(value, name)
    if not 0.0 < scalar < 1.0:
        raise ValueError(f"{name} must satisfy 0 < {name} < 1")
    return scalar


def _open_signed_unit_scalar(value: float, name: str) -> float:
    scalar = _finite_scalar(value, name)
    if not abs(scalar) < 1.0:
        raise ValueError(f"{name} must satisfy abs({name}) < 1")
    return scalar


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or int(value) != value or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)
