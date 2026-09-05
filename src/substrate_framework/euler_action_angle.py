"""Exact mode algebra for proposed Euler action-angle preparations.

The caller supplies an actual oriented volume action-angle chart, a constant
curl factor and a nonresonant frequency interval. This module constructs the
cohomological coefficient; it does not construct that chart or an Euler field.
"""

from __future__ import annotations

import sympy as sp


def cohomological_kelvin_mode(
    vorticity_mode, angular_velocity, angular_mode, action, curl_factor
):
    """Solve ``curl_factor * [u, xi] = w`` for one angular Fourier mode.

    Components are contravariant in ``(action, theta1, theta2)`` and
    ``u=(0, angular_velocity[0], angular_velocity[1])``. The supplied
    ``vorticity_mode`` contains its three coefficient functions of action.
    The returned immutable column omits the common angular exponential.

    The angular shear term is derived from the supplied velocity. A
    divergence-free curl input with zero mean gives a divergence-free lift
    under the actual cohomological hypotheses; arbitrary input coefficients
    are not asserted to satisfy those physical conditions. Symbolic nonzero
    frequency and curl-factor assumptions remain the caller's responsibility.
    """
    w = sp.ImmutableMatrix(vorticity_mode)
    omega = sp.ImmutableMatrix(angular_velocity)
    mode = sp.ImmutableMatrix(angular_mode)
    action = sp.sympify(action)
    factor = sp.sympify(curl_factor)
    if w.shape != (3, 1) or omega.shape != (2, 1) or mode.shape != (2, 1):
        raise ValueError("expected a three-component curl and two angular components")
    if not isinstance(action, sp.Symbol):
        raise ValueError("action must be a symbolic coordinate")
    if any(entry.is_integer is False for entry in mode):
        raise ValueError("angular Fourier indices must be integers")
    if factor.has(action):
        raise ValueError("this lift requires a constant curl factor")
    frequency = sp.simplify(mode.dot(omega))
    if factor.is_zero or frequency.is_zero:
        raise ValueError("curl factor and mode frequency must be nonzero")
    normal = w[0] / (sp.I * factor * frequency)
    return sp.ImmutableMatrix([
        normal,
        w[1] / (sp.I * factor * frequency)
        + normal * sp.diff(omega[0], action) / (sp.I * frequency),
        w[2] / (sp.I * factor * frequency)
        + normal * sp.diff(omega[1], action) / (sp.I * frequency),
    ])
