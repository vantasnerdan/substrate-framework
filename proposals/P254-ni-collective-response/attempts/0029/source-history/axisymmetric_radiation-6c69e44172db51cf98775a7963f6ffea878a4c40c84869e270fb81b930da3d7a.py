"""Convention-safe conditional readout of an axisymmetric STF time trace.

This module evaluates the external linearized-gravity convention used by the
accepted conditional GW claims. It does not derive a gravitational field
equation, coupling, radiation zone, backreaction law, or physical scale.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ConditionalAxisymmetricRadiationCoefficients:
    """Dimensionless coefficients after removing the declared ``G`` and ``R``.

    ``conventional_plus_R_over_G`` means ``h_plus*R/G`` and
    ``power_over_G`` means ``P/G``. The cross coefficient is provided
    explicitly so consumers cannot silently infer a second polarization.
    """

    conventional_plus_R_over_G: FloatArray
    conventional_cross_R_over_G: FloatArray
    power_over_G: FloatArray
    inclination: float
    quadrupole_scale: float


def conditional_axisymmetric_radiation_coefficients(
    axial_second_derivative: ArrayLike,
    axial_third_derivative: ArrayLike,
    *,
    inclination: float,
    quadrupole_scale: float,
) -> ConditionalAxisymmetricRadiationCoefficients:
    """Map an axisymmetric STF eigenvalue trace to conditional coefficients.

    The supplied tensor convention has symmetry-axis eigenvalue ``lambda`` and
    matrix ``diag(-lambda/2,-lambda/2,lambda)`` in its principal frame.
    ``quadrupole_scale=1`` denotes the normalized STF moment and scale three
    denotes ``Q=3*I_STF``. Under the declared normalized waveform/flux inputs,

    ``h_plus*R/G = 3*lambda_ddot*sin(inclination)^2/(2*scale)``,

    ``h_cross=0``, and ``P/G = 3*lambda'''^2/(10*scale^2)``.

    Thus a triple-normalized ``Qzz=q`` gives ``h_plus*R/G=q_ddot*sin(i)^2/2``
    and ``P/G=q'''^2/30``. Passing that same triple tensor with scale one is
    the factor-three waveform and factor-nine power convention error.
    """

    second = _finite_array(axial_second_derivative, "axial_second_derivative")
    third = _finite_array(axial_third_derivative, "axial_third_derivative")
    if second.shape != third.shape:
        raise ValueError("axial derivative arrays must have the same shape")
    angle = float(inclination)
    scale = float(quadrupole_scale)
    if not np.isfinite(angle):
        raise ValueError("inclination must be finite")
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("quadrupole_scale must be finite and nonzero")
    sine_squared = float(np.sin(angle) ** 2)
    plus = 3.0 * second * sine_squared / (2.0 * scale)
    cross = np.zeros_like(plus)
    power = 3.0 * np.square(third) / (10.0 * scale**2)
    return ConditionalAxisymmetricRadiationCoefficients(
        conventional_plus_R_over_G=np.asarray(plus, dtype=np.float64),
        conventional_cross_R_over_G=cross,
        power_over_G=np.asarray(power, dtype=np.float64),
        inclination=angle,
        quadrupole_scale=scale,
    )


def _finite_array(values: ArrayLike, name: str) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0 or array.size == 0:
        raise ValueError(f"{name} must be a nonempty array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array
