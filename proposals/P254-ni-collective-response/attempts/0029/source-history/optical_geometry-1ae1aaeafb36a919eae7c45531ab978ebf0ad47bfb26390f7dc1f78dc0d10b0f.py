"""Exact geometry for the static 1+1 optical metric.

The declared model metric is ``diag(-1/n, n/c0**2)`` for a positive static
index ``n(x)`` and positive medium signal speed ``c0``.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) > 0.0:
            raise ValueError(f"{name} must be real and positive")
    return expression


def optical_metric_1d(index: Any, signal_speed: Any) -> sp.Matrix:
    """Return ``diag(-1/n, n/c0**2)`` in coordinates ``(t, x)``."""

    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    return sp.diag(-1 / n_value, n_value / c0**2)


def optical_dilaton(index: Any) -> sp.Expr:
    """Return the canonical 1+1 optical dilaton ``log(n)``."""

    return sp.log(_positive(index, "index"))


def optical_ricci_scalar_1d(
    index: Any, coordinate: sp.Symbol, signal_speed: Any
) -> sp.Expr:
    """Return the exact Ricci scalar of the static 1+1 optical metric."""

    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    first = sp.diff(n_value, coordinate)
    second = sp.diff(n_value, coordinate, 2)
    return sp.simplify(c0**2 * (n_value * second - 2 * first**2) / n_value**3)


def optical_box_static_1d(
    scalar: Any, index: Any, coordinate: sp.Symbol, signal_speed: Any
) -> sp.Expr:
    """Return the scalar wave operator for a static scalar in the metric.

    The metric volume density is constant, so this is
    ``c0**2 * d_x[(d_x scalar)/n]``.
    """

    expression = sp.sympify(scalar)
    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    return sp.simplify(c0**2 * sp.diff(sp.diff(expression, coordinate) / n_value, coordinate))


def index_from_potential(potential: Any, signal_speed: Any) -> sp.Expr:
    """Return the conditional TF map ``n = 1/(1 + 2*Phi/c0**2)``."""

    phi = sp.sympify(potential)
    c0 = _positive(signal_speed, "signal_speed")
    denominator = sp.simplify(1 + 2 * phi / c0**2)
    if denominator.is_number:
        if denominator.is_real is not True or not float(denominator) > 0.0:
            raise ValueError("potential must produce a positive optical index")
    return 1 / denominator


def slow_geodesic_acceleration_1d(
    index: Any, coordinate: sp.Symbol, signal_speed: Any
) -> sp.Expr:
    """Return the static slow-geodesic acceleration ``-Gamma^x_tt``."""

    n_value = _positive(index, "index")
    c0 = _positive(signal_speed, "signal_speed")
    return sp.simplify(c0**2 * sp.diff(n_value, coordinate) / (2 * n_value**3))


def slow_geodesic_acceleration_from_potential(
    potential: Any, coordinate: sp.Symbol, signal_speed: Any
) -> sp.Expr:
    """Return the geodesic acceleration under the conditional TF index map."""

    n_value = index_from_potential(potential, signal_speed)
    return slow_geodesic_acceleration_1d(n_value, coordinate, signal_speed)


def optical_dilaton_source_operator_1d(
    potential: Any, coordinate: sp.Symbol, signal_speed: Any
) -> sp.Expr:
    """Return ``-Box_g(log(n(Phi)))`` under the conditional TF index map.

    This is the geometric source-side operator, not a matter field equation.
    In the accepted 1+1 conventions it simplifies exactly to ``2*Phi_xx``.
    """

    n_value = index_from_potential(potential, signal_speed)
    dilaton = optical_dilaton(n_value)
    return sp.simplify(
        -optical_box_static_1d(dilaton, n_value, coordinate, signal_speed)
    )
