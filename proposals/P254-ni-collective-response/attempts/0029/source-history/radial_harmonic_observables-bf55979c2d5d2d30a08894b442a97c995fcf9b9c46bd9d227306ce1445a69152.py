"""Energy observables for finite radial sine-Gordon harmonic expansions.

The functions in this module operate on declared harmonic coefficients. They
do not turn a finite-harmonic or finite-box approximation into an exact PDE
solution. In particular, a temporal line in a spherically symmetric scalar
energy moment is not a nonzero STF quadrupole or a physical radiation result.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .numerics import trapezoid_integral

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class RadialHarmonicKinematics:
    """Field, time derivative, and radial derivative on a phase-radius grid."""

    phase: FloatArray
    field: FloatArray
    time_derivative: FloatArray
    radial_derivative: FloatArray


@dataclass(frozen=True)
class PeriodicFourierCoefficients:
    """Real Fourier coefficients for endpoint-excluded uniform samples.

    The represented series is ``a_0 + sum_k(a_k*cos(k*tau)+b_k*sin(k*tau))``.
    Hence the DC row is the sample mean and positive-frequency rows use the
    normalization ``2/N``.
    """

    cosine: FloatArray
    sine: FloatArray
    sample_count: int

    @property
    def amplitude(self) -> FloatArray:
        """Return ``sqrt(a_k**2+b_k**2)`` for every stored harmonic."""

        return np.asarray(
            np.sqrt(np.square(self.cosine) + np.square(self.sine)),
            dtype=np.float64,
        )


def radial_harmonic_kinematics(
    amplitudes: ArrayLike,
    radial_derivatives: ArrayLike,
    harmonics: tuple[int, ...],
    frequency: float,
    phase: ArrayLike,
) -> RadialHarmonicKinematics:
    """Reconstruct harmonic field kinematics on a phase-radius product grid.

    ``phase`` is the dimensionless angle ``tau=frequency*t``. The returned
    time derivative is the physical derivative with respect to ``t`` and
    therefore includes both the harmonic index and ``frequency``.
    """

    modes = _harmonics(harmonics)
    coefficients = _coefficient_matrix(amplitudes, len(modes), "amplitudes")
    derivatives = _coefficient_matrix(
        radial_derivatives,
        len(modes),
        "radial_derivatives",
        radial_size=coefficients.shape[1],
    )
    omega = _positive_finite(frequency, "frequency")
    phases = _finite_vector(phase, "phase")
    mode_array = np.asarray(modes, dtype=np.float64)
    angles = np.outer(mode_array, phases)
    cosine = np.cos(angles)
    sine = np.sin(angles)
    field = cosine.T @ coefficients
    radial = cosine.T @ derivatives
    time = -(omega * mode_array[:, None] * sine).T @ coefficients
    return RadialHarmonicKinematics(
        phase=np.asarray(phases, dtype=np.float64),
        field=np.asarray(field, dtype=np.float64),
        time_derivative=np.asarray(time, dtype=np.float64),
        radial_derivative=np.asarray(radial, dtype=np.float64),
    )


def radial_harmonic_energy_density(
    amplitudes: ArrayLike,
    radial_derivatives: ArrayLike,
    harmonics: tuple[int, ...],
    frequency: float,
    phase: ArrayLike,
) -> FloatArray:
    """Return canonical radial sine-Gordon ``T00`` on a phase-radius grid."""

    kinematics = radial_harmonic_kinematics(
        amplitudes,
        radial_derivatives,
        harmonics,
        frequency,
        phase,
    )
    return np.asarray(
        0.5 * np.square(kinematics.time_derivative)
        + 0.5 * np.square(kinematics.radial_derivative)
        + 1.0
        - np.cos(kinematics.field),
        dtype=np.float64,
    )


def periodic_fourier_coefficients(
    samples: ArrayLike,
    *,
    max_harmonic: int,
) -> PeriodicFourierCoefficients:
    """Return direct real Fourier coefficients along the leading sample axis.

    Samples must cover one full period on a uniform endpoint-excluded phase
    grid. The maximum stored harmonic is restricted below the Nyquist index,
    where the generic positive-frequency ``2/N`` normalization remains valid.
    """

    values = np.asarray(samples, dtype=np.float64)
    if values.ndim == 0 or values.shape[0] < 4 or not np.all(np.isfinite(values)):
        raise ValueError(
            "samples must have at least four finite values along the leading axis"
        )
    maximum = int(max_harmonic)
    if maximum != max_harmonic or maximum < 0:
        raise ValueError("max_harmonic must be a nonnegative integer")
    count = values.shape[0]
    if maximum >= count // 2:
        raise ValueError("max_harmonic must be strictly below the Nyquist index")
    phase = 2.0 * np.pi * np.arange(count, dtype=np.float64) / count
    harmonic = np.arange(maximum + 1, dtype=np.float64)
    cosine_basis = np.cos(np.outer(harmonic, phase))
    sine_basis = np.sin(np.outer(harmonic, phase))
    flat = values.reshape(count, -1)
    cosine = (2.0 / count) * cosine_basis @ flat
    sine = (2.0 / count) * sine_basis @ flat
    cosine[0] *= 0.5
    sine[0] = 0.0
    output_shape = (maximum + 1, *values.shape[1:])
    return PeriodicFourierCoefficients(
        cosine=np.asarray(cosine.reshape(output_shape), dtype=np.float64),
        sine=np.asarray(sine.reshape(output_shape), dtype=np.float64),
        sample_count=count,
    )


def integrate_spherical_radial_density(
    radius: ArrayLike,
    density: ArrayLike,
    *,
    radial_power: int = 0,
) -> float | FloatArray:
    """Integrate ``4*pi*integral r**(2+radial_power)*density dr``.

    The final density axis is radial; any leading axes are retained. Power
    zero is total spherical energy and power two is its scalar radial second
    moment. Quadrature dispatch is centralized in the shared NumPy-version-
    independent helper.
    """

    coordinate = _radius(radius)
    values = np.asarray(density, dtype=np.float64)
    if values.ndim == 0 or values.shape[-1] != coordinate.size:
        raise ValueError("density must have one final-axis value per radius")
    if not np.all(np.isfinite(values)):
        raise ValueError("density must contain only finite values")
    power = int(radial_power)
    if power != radial_power or power < 0:
        raise ValueError("radial_power must be a nonnegative integer")
    weight = 4.0 * np.pi * np.power(coordinate, 2 + power)
    flat = values.reshape(-1, coordinate.size)
    result = np.asarray(
        [trapezoid_integral(row * weight, coordinate) for row in flat],
        dtype=np.float64,
    ).reshape(values.shape[:-1])
    if result.ndim == 0:
        return float(result)
    return result


def spherical_radial_second_moment_tensor(
    radius: ArrayLike,
    density: ArrayLike,
) -> FloatArray:
    """Return the isotropic Cartesian second moment of a radial density.

    If ``M2=integral rho*r**2*d^3x``, exact angular integration gives
    ``I_ij=(M2/3)*delta_ij``. The returned array has any leading density axes
    followed by the two tensor axes.
    """

    scalar = np.asarray(
        integrate_spherical_radial_density(radius, density, radial_power=2),
        dtype=np.float64,
    )
    return np.asarray(
        scalar[..., None, None] * np.eye(3, dtype=np.float64) / 3.0,
        dtype=np.float64,
    )


def time_averaged_per_axis_energy_variance(
    radius: ArrayLike,
    energy_density: ArrayLike,
) -> float:
    """Return ``<r**2>/3`` for a positive time-averaged radial energy density."""

    values = np.asarray(energy_density, dtype=np.float64)
    coordinate = _radius(radius)
    if values.ndim != 2 or values.shape[1] != coordinate.size:
        raise ValueError("energy_density must have phase rows and radius columns")
    if not np.all(np.isfinite(values)):
        raise ValueError("energy_density must contain only finite values")
    averaged = np.mean(values, axis=0)
    energy = float(integrate_spherical_radial_density(coordinate, averaged))
    if energy <= 0.0:
        raise ValueError("time-averaged total energy must be positive")
    second = float(
        integrate_spherical_radial_density(
            coordinate,
            averaged,
            radial_power=2,
        )
    )
    return second / (3.0 * energy)


def _harmonics(harmonics: tuple[int, ...]) -> tuple[int, ...]:
    modes = tuple(int(value) for value in harmonics)
    if not modes or modes != harmonics:
        raise ValueError("harmonics must be a nonempty tuple of integers")
    if any(mode <= 0 for mode in modes) or tuple(sorted(set(modes))) != modes:
        raise ValueError("harmonics must be positive, increasing, and unique")
    return modes


def _coefficient_matrix(
    values: ArrayLike,
    harmonic_count: int,
    name: str,
    *,
    radial_size: int | None = None,
) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1:
        array = array[:, None]
    if array.ndim != 2 or array.shape[0] != harmonic_count:
        raise ValueError(f"{name} must have one row per harmonic")
    if radial_size is not None and array.shape[1] != radial_size:
        raise ValueError(f"{name} must match the radial sample count")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _finite_vector(values: ArrayLike, name: str) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite one-dimensional array")
    return array


def _radius(values: ArrayLike) -> FloatArray:
    radius = _finite_vector(values, "radius")
    if radius.size < 2 or radius[0] < 0.0 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("radius must be strictly increasing and nonnegative")
    return radius


def _positive_finite(value: float, name: str) -> float:
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number
