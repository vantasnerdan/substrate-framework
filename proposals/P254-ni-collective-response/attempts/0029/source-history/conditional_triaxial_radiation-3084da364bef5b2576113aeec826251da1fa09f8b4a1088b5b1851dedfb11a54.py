"""Convention-safe conditional radiation algebra for real-m2 STF moments.

Every waveform and power result in this module is conditional on the wave
equation, far-zone Green function, and flux premises recorded by C-GW-001.
The functions do not derive gravity, a source deformation, or a graviton
interpretation from the scalar framework.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .triaxial_l2 import real_l2_tt_readout
from .tt_angular import conditional_tt_power


@dataclass(frozen=True)
class ConditionalScaledSTFWaveform:
    """Conditional TT waveform for ``Q_s = s*I_STF``.

    Coordinates labeled ``normalized`` contract against unit-Frobenius TT
    basis tensors. Conventional readouts are smaller by ``sqrt(2)`` and match
    ``h_plus=(h_11-h_22)/2`` and ``h_cross=h_12`` in the transverse frame.
    """

    quadrupole_scale: sp.Expr
    gravitational_coupling: sp.Expr
    distance: sp.Expr
    waveform_prefactor: sp.Expr
    quadrupole_second_derivative: sp.Matrix
    projected_quadrupole_second_derivative: sp.Matrix
    waveform_tensor: sp.Matrix
    normalized_plus: sp.Expr
    normalized_cross: sp.Expr
    conventional_plus: sp.Expr
    conventional_cross: sp.Expr


def _stf_three(tensor: Any, name: str) -> sp.Matrix:
    value = sp.Matrix(tensor)
    if value.shape != (3, 3):
        raise ValueError(f"{name} must be 3 by 3")
    if sp.simplify(value - value.T) != sp.zeros(3):
        raise ValueError(f"{name} must be symmetric")
    if sp.simplify(sp.trace(value)) != 0:
        raise ValueError(f"{name} must be trace free")
    return value


def _nonzero(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if sp.simplify(expression) == 0:
        raise ValueError(f"{name} must be nonzero")
    return expression


def conditional_scaled_stf_waveform(
    quadrupole_second_derivative: Any,
    direction: Any,
    gravitational_coupling: Any,
    distance: Any,
    quadrupole_scale: Any = 1,
    reference: Any | None = None,
) -> ConditionalScaledSTFWaveform:
    """Return the conditional waveform for a scaled STF quadrupole.

    The declared normalized-convention premise is
    ``h_TT=(2*G/R)*TT(I_STF'')``. If the supplied tensor is
    ``Q_s=s*I_STF``, the coefficient multiplying it is ``2*G/(s*R)``.
    """

    derivative = _stf_three(
        quadrupole_second_derivative, "quadrupole_second_derivative"
    )
    coupling = sp.sympify(gravitational_coupling)
    radius = _nonzero(distance, "distance")
    scale = _nonzero(quadrupole_scale, "quadrupole_scale")
    readout = real_l2_tt_readout(derivative, direction, reference)
    prefactor = sp.simplify(2 * coupling / (scale * radius))
    return ConditionalScaledSTFWaveform(
        quadrupole_scale=scale,
        gravitational_coupling=coupling,
        distance=radius,
        waveform_prefactor=prefactor,
        quadrupole_second_derivative=derivative,
        projected_quadrupole_second_derivative=readout.projected_tensor,
        waveform_tensor=sp.simplify(prefactor * readout.projected_tensor),
        normalized_plus=sp.simplify(
            prefactor * readout.normalized_plus_coordinate
        ),
        normalized_cross=sp.simplify(
            prefactor * readout.normalized_cross_coordinate
        ),
        conventional_plus=sp.simplify(
            prefactor * readout.conventional_plus_readout
        ),
        conventional_cross=sp.simplify(
            prefactor * readout.conventional_cross_readout
        ),
    )


def conditional_scaled_stf_power(
    quadrupole_third_derivative: Any,
    gravitational_coupling: Any,
    quadrupole_scale: Any = 1,
) -> sp.Expr:
    """Return ``G*|Q_s'''|_F^2/(5*s^2)`` under C-GW-001 premises."""

    derivative = _stf_three(
        quadrupole_third_derivative, "quadrupole_third_derivative"
    )
    coupling = _nonzero(gravitational_coupling, "gravitational_coupling")
    scale = _nonzero(quadrupole_scale, "quadrupole_scale")
    return conditional_tt_power(
        derivative,
        2 * coupling / scale,
        1 / (32 * sp.pi * coupling),
    )


def real_m2_triple_stf_tensor(
    cosine_component: Any,
    sine_component: Any = 0,
) -> sp.Matrix:
    """Return a triple-STF real-m2 tensor from its Cartesian components.

    ``cosine_component`` is ``Q_xx=-Q_yy`` and ``sine_component`` is
    ``Q_xy=Q_yx``. In C-GW-007's radial-coefficient notation these equal
    ``2*H_2c/5`` and ``2*H_2s/5`` respectively.
    """

    cosine = sp.sympify(cosine_component)
    sine = sp.sympify(sine_component)
    return sp.Matrix(
        [[cosine, sine, 0], [sine, -cosine, 0], [0, 0, 0]]
    )


def conditional_real_m2_natural_axis_waveform(
    cosine_second_derivative: Any,
    sine_second_derivative: Any,
    gravitational_coupling: Any,
    distance: Any,
) -> ConditionalScaledSTFWaveform:
    """Return the triple-STF real-m2 waveform viewed along its z axis."""

    return conditional_scaled_stf_waveform(
        real_m2_triple_stf_tensor(
            cosine_second_derivative, sine_second_derivative
        ),
        [0, 0, 1],
        gravitational_coupling,
        distance,
        quadrupole_scale=3,
        reference=[1, 0, 0],
    )


def conditional_real_m2_power(
    cosine_third_derivative: Any,
    sine_third_derivative: Any,
    gravitational_coupling: Any,
) -> sp.Expr:
    """Return ``2*G*(q_c'''**2+q_s'''**2)/45`` conditionally."""

    return conditional_scaled_stf_power(
        real_m2_triple_stf_tensor(
            cosine_third_derivative, sine_third_derivative
        ),
        gravitational_coupling,
        quadrupole_scale=3,
    )
