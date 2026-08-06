"""Full-retardation Einstein radiation from one sine-Gordon breather stress.

The calculation uses the exact longitudinal pressure and does not assume that
the breather is slow compared with its profile size.
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
from .linearized_einstein import longitudinal_line_radiation_ledger
from .numerics import trapezoid_integral, trapezoid_integral_axis
from .sine_gordon import breather_inverse_width
from .sine_gordon_weave import SineGordonWeaveGravity, sine_gordon_weave_gravity


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
class BreatherRadiationEvidence:
    """Resolution-stamped full-retardation power of the exact breather stress."""

    frequency: sp.Expr
    inverse_width: float
    spatial_extent_in_inverse_widths: float
    spatial_points: int
    time_points: int
    direction_points: int
    harmonic_count: int
    source_size_frequency_ratio: float
    source_energy: sp.Expr
    newton_constant: sp.Expr
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
    """Evaluate the exact normalized breather ``T_xx`` on a tensor grid."""

    sech = 1.0 / np.cosh(inverse_width_coordinate)[:, None]
    sine = np.sin(phase)[None, :]
    cosine = np.cos(phase)[None, :]
    tangent = np.tanh(inverse_width_coordinate)[:, None]
    tangent_argument = (inverse_width / frequency) * sech * sine
    denominator = 1.0 + tangent_argument**2
    field_time = 4.0 * inverse_width * sech * cosine / denominator
    field_space = (
        -4.0
        * inverse_width
        * tangent_argument
        * tangent
        / denominator
    )
    potential = 8.0 * tangent_argument**2 / denominator**2
    return 0.5 * (field_time**2 + field_space**2) - potential


def leading_small_amplitude_breather_power_ratio(
    frequency: Any,
    *,
    direction_points: int = 4001,
) -> float:
    r"""Return the finite-wavelength leading power ``P/(mu*c)``.

    At small inverse width ``eta``, the longitudinal pressure is
    ``8*eta**2*sech(eta*X)**2*cos(2*omega*tau)``.  Its exact spatial Fourier
    transform gives

    ``64*pi**2*omega**4/log(2) * integral_-1^1
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
    integral = 2.0 * trapezoid_integral(integrand, direction)
    return float(
        64.0
        * np.pi**2
        * numeric_frequency**4
        * integral
        / np.log(2.0)
    )


def breather_radiation_evidence(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
    *,
    spatial_extent_in_inverse_widths: float = 15.0,
    spatial_points: int = 4097,
    time_points: int = 257,
    direction_points: int = 2001,
    harmonic_count: int = 6,
) -> BreatherRadiationEvidence:
    r"""Compute the exact-stress, full-retardation breather power.

    The exact pressure is Fourier decomposed over its period.  For every
    direction cosine ``zeta``, each temporal harmonic is spatially transformed
    at its actual retarded wavenumber.  Orthogonality performs the time
    average, and the accepted TT projection plus Isaacson flux reduce the
    result to

    ``P/(mu*c)=integral_-1^1 (1-zeta**2)**2
    <(d_tau integral dX T_xx(tau+zeta*X,X))**2> dzeta/(8*log(2))``.

    This is a one-way leading radiation prediction: backreaction is controlled
    a posteriori by the reported energy loss per field cycle.
    """

    exact_frequency, numeric_frequency, inverse_width = _numeric_frequency(frequency)
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

    inverse_width_coordinate = np.linspace(-extent, extent, space_count)
    longitudinal_coordinate = inverse_width_coordinate / inverse_width
    phase = np.linspace(0.0, np.pi, time_count)
    pressure = _normalized_breather_longitudinal_stress(
        inverse_width_coordinate,
        phase,
        numeric_frequency,
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

    direction = np.linspace(0.0, 1.0, direction_count)
    mean_retarded_derivative_squared = np.zeros_like(direction)
    chunk_size = 128
    for index, mode in enumerate(mode_numbers):
        spatial_amplitude = np.empty_like(direction)
        wavenumber = 2.0 * mode * numeric_frequency
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
        derivative_amplitude = wavenumber * spatial_amplitude
        mean_retarded_derivative_squared += 0.5 * derivative_amplitude**2

    angular_integrand = (1.0 - direction**2) ** 2 * mean_retarded_derivative_squared
    leading_power = leading_small_amplitude_breather_power_ratio(
        exact_frequency,
        direction_points=direction_count,
    )
    asymptotic_power = (
        64.0
        * np.pi
        * inverse_width**3
        * numeric_frequency
        / (3.0 * np.log(2.0))
    )

    weave = sine_gordon_weave_gravity(coefficients)
    observables = dimensional_breather_observables(exact_frequency, coefficients)
    onsite = _finite_positive_float(coefficients.onsite, "coefficients.onsite")
    speed = _finite_positive_float(weave.scales.signal_speed, "signal_speed")
    power_scale = onsite * speed
    center_direction = longitudinal_line_radiation_ledger(
        weave.horizon_equilibrium.newton_constant,
        weave.scales.signal_speed,
        0,
    )
    center_power_coefficient = _finite_positive_float(
        center_direction.azimuth_power_coefficient,
        "line radiation coefficient",
    )
    physical_stress_rate_scale = power_scale
    power = float(
        2.0
        * center_power_coefficient
        * physical_stress_rate_scale**2
        * trapezoid_integral(angular_integrand, direction)
    )
    normalized_power = float(power / power_scale)
    reduced_power = float(
        2.0
        * trapezoid_integral(angular_integrand, direction)
        / (8.0 * np.log(2.0))
    )
    if not np.isclose(normalized_power, reduced_power, rtol=2.0e-14, atol=0.0):
        raise AssertionError("retarded Einstein coefficient did not reduce upstream")
    period = _finite_positive_float(observables.period, "breather period")
    source_energy = _finite_positive_float(observables.energy, "breather energy")
    loss_fraction = power * period / source_energy
    harmonic_maxima = tuple(
        float(np.max(np.abs(pressure_modes[index]))) for index in range(modes)
    )
    return BreatherRadiationEvidence(
        frequency=exact_frequency,
        inverse_width=inverse_width,
        spatial_extent_in_inverse_widths=extent,
        spatial_points=space_count,
        time_points=time_count,
        direction_points=direction_count,
        harmonic_count=modes,
        source_size_frequency_ratio=numeric_frequency / inverse_width,
        source_energy=observables.energy,
        newton_constant=weave.horizon_equilibrium.newton_constant,
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
