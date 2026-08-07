"""Full-retardation radiation from the exact sine-Gordon breather stress.

The source is the longitudinal component of the same conserved stress whose
rest energy is used by ``sine_gordon_breather_gravity``.  The calculation does
not make a long-wavelength or slow-source approximation: every retained time
harmonic is spatially transformed at its actual retarded wavenumber.  The
Newton coefficient is derived by ``sine_gordon_teleparallel`` before this
numeric source transform is evaluated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_breather_observables,
)
from .linearized_einstein import retarded_line_source_radiation_ledger
from .numerics import trapezoid_integral, trapezoid_integral_axis
from .sine_gordon import breather_inverse_width
from .sine_gordon_teleparallel import (
    SineGordonTeleparallelGravity,
    sine_gordon_teleparallel_gravity,
)


def _positive_integer(value: int, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}; got {value!r}")
    return value


def _odd_grid_size(value: int, name: str, minimum: int) -> int:
    points = _positive_integer(value, name, minimum)
    if points % 2 == 0:
        raise ValueError(f"{name} must be odd for a symmetry-resolving grid; got {points}")
    return points


def _numeric_frequency(value: Any) -> tuple[sp.Expr, float, float]:
    exact = sp.sympify(value)
    if exact.has(sp.Float):
        raise ValueError("frequency must be exact rather than floating")
    if exact.is_real is not True or exact.is_positive is not True:
        raise ValueError("frequency must be provably positive and real")
    if exact.is_number and not bool(exact < 1):
        raise ValueError("frequency must satisfy 0 < frequency < 1")
    numeric = float(sp.N(exact, 17))
    inverse_width = float(sp.N(breather_inverse_width(exact), 17))
    if not np.isfinite(numeric) or not np.isfinite(inverse_width):
        raise ValueError(f"frequency produced non-finite values; got {exact!s}")
    return exact, numeric, inverse_width


def _finite_positive_float(value: Any, name: str) -> float:
    try:
        numeric = float(sp.N(sp.sympify(value), 17))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be numerically evaluable; got {value!r}") from error
    if not np.isfinite(numeric) or numeric <= 0.0:
        raise ValueError(f"{name} must be positive and finite; got {value!r}")
    return numeric


@dataclass(frozen=True)
class TeleparallelBreatherRadiationEvidence:
    """Resolution-stamped full-retardation power of one exact SG stress."""

    frequency: sp.Expr
    inverse_width: float
    spatial_extent_in_inverse_widths: float
    spatial_points: int
    time_points: int
    direction_points: int
    harmonic_count: int
    source_size_frequency_ratio: float
    source_energy: sp.Expr
    gravity: SineGordonTeleparallelGravity
    power_scale: float
    power: float
    power_over_onsite_speed: float
    leading_small_amplitude_power_ratio: float
    small_width_asymptotic_power_ratio: float
    relative_leading_difference: float
    field_cycle_energy_loss_fraction: float
    harmonic_pressure_maxima: tuple[float, ...]


def _normalized_breather_longitudinal_stress(
    inverse_width_coordinate: np.ndarray,
    phase: np.ndarray,
    frequency: float,
    inverse_width: float,
) -> np.ndarray:
    """Evaluate ``T_xx/mu`` from the exact accepted breather on a tensor grid."""

    sech = 1.0 / np.cosh(inverse_width_coordinate)[:, None]
    sine = np.sin(phase)[None, :]
    cosine = np.cos(phase)[None, :]
    tangent = np.tanh(inverse_width_coordinate)[:, None]
    tangent_argument = (inverse_width / frequency) * sech * sine
    denominator = 1.0 + tangent_argument**2
    field_time = 4.0 * inverse_width * sech * cosine / denominator
    field_space = -4.0 * inverse_width * tangent_argument * tangent / denominator
    potential = 8.0 * tangent_argument**2 / denominator**2
    return 0.5 * (field_time**2 + field_space**2) - potential


def leading_small_amplitude_teleparallel_breather_power_ratio(
    frequency: Any,
    *,
    direction_points: int = 4001,
) -> float:
    r"""Return the leading finite-wavelength ``P/(mu*c)`` exactly in width.

    For inverse width ``eta``, the leading pressure is
    ``8*eta**2*sech(eta*X)**2*cos(2*omega*tau)``.  Its exact Fourier transform
    and the derived coupling ``G*mu/c**4=1/(8*pi)`` give

    ``32*pi*omega**4 * integral_-1^1
    zeta**2*(1-zeta**2)**2/sinh(pi*omega*zeta/eta)**2 dzeta``.
    """

    _, numeric_frequency, inverse_width = _numeric_frequency(frequency)
    points = _odd_grid_size(direction_points, "direction_points", 129)
    direction = np.linspace(0.0, 1.0, points, dtype=np.float64)
    argument = np.pi * numeric_frequency * direction / inverse_width
    ratio = np.empty_like(direction)
    ratio[0] = (inverse_width / (np.pi * numeric_frequency)) ** 2
    with np.errstate(over="ignore", under="ignore"):
        ratio[1:] = (direction[1:] / np.sinh(argument[1:])) ** 2
    integrand = ratio * (1.0 - direction**2) ** 2
    full_integral = 2.0 * trapezoid_integral(integrand, direction)
    return float(32.0 * np.pi * numeric_frequency**4 * full_integral)


def _validated_radiation_grid(
    spatial_extent_in_inverse_widths: float,
    spatial_points: int,
    time_points: int,
    direction_points: int,
    harmonic_count: int,
) -> tuple[float, int, int, int, int]:
    """Validate and return all resolution controls together."""

    extent = float(spatial_extent_in_inverse_widths)
    if not np.isfinite(extent) or extent <= 0.0:
        raise ValueError(
            "spatial_extent_in_inverse_widths must be positive and finite; "
            f"got {spatial_extent_in_inverse_widths!r}"
        )
    space_count = _odd_grid_size(spatial_points, "spatial_points", 257)
    time_count = _odd_grid_size(time_points, "time_points", 33)
    direction_count = _odd_grid_size(direction_points, "direction_points", 129)
    modes = _positive_integer(harmonic_count, "harmonic_count", 1)
    if time_count < 4 * modes + 1:
        raise ValueError(
            "time_points must provide at least four samples per retained harmonic "
            f"plus the periodic endpoint; got {time_count} points and {modes} harmonics"
        )
    return extent, space_count, time_count, direction_count, modes


def _breather_pressure_harmonics(
    frequency: float,
    inverse_width: float,
    extent: float,
    space_count: int,
    time_count: int,
    modes: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project the exact normalized pressure onto temporal cosine modes."""

    inverse_width_coordinate = np.linspace(-extent, extent, space_count)
    longitudinal_coordinate = inverse_width_coordinate / inverse_width
    phase = np.linspace(0.0, np.pi, time_count)
    pressure = _normalized_breather_longitudinal_stress(
        inverse_width_coordinate,
        phase,
        frequency,
        inverse_width,
    )
    mode_numbers = np.arange(1, modes + 1, dtype=np.float64)
    temporal_basis = np.cos(2.0 * mode_numbers[:, None] * phase[None, :])
    pressure_modes = np.empty((modes, space_count), dtype=np.float64)
    for index in range(modes):
        pressure_modes[index] = (
            2.0
            / np.pi
            * trapezoid_integral_axis(
                pressure * temporal_basis[index][None, :],
                phase,
                axis=1,
            )
        )
    return longitudinal_coordinate, mode_numbers, pressure_modes


def _retarded_derivative_spectrum(
    frequency: float,
    longitudinal_coordinate: np.ndarray,
    mode_numbers: np.ndarray,
    pressure_modes: np.ndarray,
    direction_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``<dot(H_zeta)**2>`` after full spatial retardation."""

    direction = np.linspace(0.0, 1.0, direction_count)
    mean_derivative_squared = np.zeros_like(direction)
    chunk_size = 128
    for index, mode in enumerate(mode_numbers):
        spatial_amplitude = np.empty_like(direction)
        wavenumber = 2.0 * mode * frequency
        for start in range(0, direction_count, chunk_size):
            stop = min(start + chunk_size, direction_count)
            phase_matrix = (
                wavenumber
                * direction[start:stop, None]
                * longitudinal_coordinate[None, :]
            )
            spatial_amplitude[start:stop] = trapezoid_integral_axis(
                pressure_modes[index][None, :] * np.cos(phase_matrix),
                longitudinal_coordinate,
                axis=1,
            )
        mean_derivative_squared += 0.5 * (wavenumber * spatial_amplitude) ** 2
    return direction, mean_derivative_squared


def teleparallel_breather_radiation_evidence(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    *,
    spatial_extent_in_inverse_widths: float = 15.0,
    spatial_points: int = 4097,
    time_points: int = 257,
    direction_points: int = 2001,
    harmonic_count: int = 6,
) -> TeleparallelBreatherRadiationEvidence:
    r"""Compute exact-stress full-retardation breather radiation.

    In normalized variables ``X=x/ell`` and ``tau=c*t/ell``, write
    ``H_zeta=int dX (T_xx/mu)(tau+zeta*X,X)``.  The accepted TT projector,
    no-incoming Green function, Isaacson flux, and the independently derived
    Newton coefficient reduce the total power to

    ``P/(mu*c)=1/(16*pi) * integral_-1^1 (1-zeta**2)**2
    <(d_tau H_zeta)**2> dzeta``.

    The returned loss per source cycle is the one-way backreaction gate.
    """

    exact_frequency, numeric_frequency, inverse_width = _numeric_frequency(frequency)
    extent, space_count, time_count, direction_count, modes = _validated_radiation_grid(
        spatial_extent_in_inverse_widths,
        spatial_points,
        time_points,
        direction_points,
        harmonic_count,
    )
    longitudinal_coordinate, mode_numbers, pressure_modes = (
        _breather_pressure_harmonics(
            numeric_frequency,
            inverse_width,
            extent,
            space_count,
            time_count,
            modes,
        )
    )
    direction, mean_retarded_derivative_squared = _retarded_derivative_spectrum(
        numeric_frequency,
        longitudinal_coordinate,
        mode_numbers,
        pressure_modes,
        direction_count,
    )

    angular_integrand = (1.0 - direction**2) ** 2 * mean_retarded_derivative_squared
    leading_power = leading_small_amplitude_teleparallel_breather_power_ratio(
        exact_frequency,
        direction_points=direction_count,
    )
    asymptotic_power = 32.0 * inverse_width**3 * numeric_frequency / 3.0

    gravity = sine_gordon_teleparallel_gravity(coefficients)
    observables = dimensional_breather_observables(exact_frequency, coefficients)
    onsite = _finite_positive_float(coefficients.onsite, "coefficients.onsite")
    speed = _finite_positive_float(gravity.scales.signal_speed, "signal_speed")
    power_scale = onsite * speed
    center_direction = retarded_line_source_radiation_ledger(
        gravity.newton_constant,
        gravity.scales.signal_speed,
        0,
    )
    center_power_coefficient = _finite_positive_float(
        center_direction.azimuth_power_coefficient,
        "line radiation coefficient",
    )
    physical_stress_rate_scale = power_scale
    half_direction_integral = trapezoid_integral(angular_integrand, direction)
    power = float(
        2.0
        * center_power_coefficient
        * physical_stress_rate_scale**2
        * half_direction_integral
    )
    normalized_power = float(power / power_scale)
    independently_reduced_power = float(half_direction_integral / (8.0 * np.pi))
    if not np.isclose(
        normalized_power,
        independently_reduced_power,
        rtol=2.0e-14,
        atol=0.0,
    ):
        raise AssertionError("retarded Einstein coefficient did not reduce upstream")
    period = _finite_positive_float(observables.period, "breather period")
    source_energy = _finite_positive_float(observables.energy, "breather energy")
    loss_fraction = power * period / source_energy
    harmonic_maxima = tuple(
        float(np.max(np.abs(pressure_modes[index]))) for index in range(modes)
    )
    return TeleparallelBreatherRadiationEvidence(
        frequency=exact_frequency,
        inverse_width=inverse_width,
        spatial_extent_in_inverse_widths=extent,
        spatial_points=space_count,
        time_points=time_count,
        direction_points=direction_count,
        harmonic_count=modes,
        source_size_frequency_ratio=numeric_frequency / inverse_width,
        source_energy=observables.energy,
        gravity=gravity,
        power_scale=power_scale,
        power=power,
        power_over_onsite_speed=normalized_power,
        leading_small_amplitude_power_ratio=float(leading_power),
        small_width_asymptotic_power_ratio=float(asymptotic_power),
        relative_leading_difference=float(
            abs(normalized_power - leading_power) / leading_power
        ),
        field_cycle_energy_loss_fraction=float(loss_fraction),
        harmonic_pressure_maxima=harmonic_maxima,
    )
