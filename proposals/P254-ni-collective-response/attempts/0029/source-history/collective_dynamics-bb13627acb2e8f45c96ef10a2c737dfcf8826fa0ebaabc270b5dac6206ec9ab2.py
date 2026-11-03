"""Exact one-coordinate dynamics for the declared static optical action."""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _nonnegative(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return expression


def optical_collective_lagrangian(
    coordinate: sp.Expr,
    parameter: sp.Symbol,
    index: Any,
    signal_speed: Any,
    energy_scale: Any,
) -> sp.Expr:
    """Return the declared timelike optical action integrand.

    The result is ``-E0/sqrt(n) * sqrt(1-n**2*qdot**2/c0**2)``.
    """

    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    e0 = _positive(energy_scale, "energy_scale")
    velocity = sp.diff(coordinate, parameter)
    return -e0 / sp.sqrt(n_value) * sp.sqrt(
        1 - n_value**2 * velocity**2 / c0**2
    )


def optical_collective_acceleration(
    coordinate: sp.Expr,
    parameter: sp.Symbol,
    index: Any,
    signal_speed: Any,
) -> sp.Expr:
    """Return the exact coordinate-time acceleration of the optical action."""

    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    velocity = sp.diff(coordinate, parameter)
    index_gradient = sp.diff(n_value, coordinate)
    return sp.simplify(
        (c0**2 - 3 * n_value**2 * velocity**2)
        * index_gradient
        / (2 * n_value**3)
    )


def slow_optical_collective_acceleration(
    coordinate: sp.Expr,
    index: Any,
    signal_speed: Any,
) -> sp.Expr:
    """Return the zero-velocity limit ``c0**2*n_q/(2*n**3)``."""

    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    return sp.simplify(
        c0**2 * sp.diff(n_value, coordinate) / (2 * n_value**3)
    )


def slow_optical_profile_width_correction(
    coordinate: sp.Expr,
    index: Any,
    signal_speed: Any,
    spatial_variance: Any,
) -> sp.Expr:
    """Return the centered profile-average term ``variance*a''/2``.

    Here ``a`` is :func:`slow_optical_collective_acceleration`.  The result is
    the second-order term for a declared normalized centered scaled-profile
    average.  It is not a Mathisson-Papapetrou-Dixon force law.
    """

    variance = _nonnegative(spatial_variance, "spatial_variance")
    point_acceleration = slow_optical_collective_acceleration(
        coordinate, index, signal_speed
    )
    return sp.simplify(
        variance * sp.diff(point_acceleration, coordinate, 2) / 2
    )


def virial_scaling_exponents(
    quadratic_exponent: Any,
    quartic_exponent: Any,
) -> tuple[sp.Expr, sp.Expr]:
    """Return conditional ``(width_slope, energy_slope)`` formulas."""

    a = sp.sympify(quadratic_exponent)
    b = sp.sympify(quartic_exponent)
    return sp.simplify((a - b) / 2), sp.simplify(-(a + b) / 2)
