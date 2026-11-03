"""Conditional exact consequences of declared rectangular Wilson-loop laws."""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def rectangular_area_law(
    separation: Any, euclidean_time: Any, string_tension: Any
) -> sp.Expr:
    """Return the declared ansatz ``exp(-string_tension*R*T)``."""

    distance = _positive(separation, "separation")
    duration = _positive(euclidean_time, "euclidean_time")
    tension = _positive(string_tension, "string_tension")
    return sp.exp(-tension * distance * duration)


def rectangular_perimeter_law(
    separation: Any, euclidean_time: Any, perimeter_coefficient: Any
) -> sp.Expr:
    """Return the declared ansatz ``exp(-2*coefficient*(R+T))``."""

    distance = _positive(separation, "separation")
    duration = _positive(euclidean_time, "euclidean_time")
    coefficient = _positive(perimeter_coefficient, "perimeter_coefficient")
    return sp.exp(-2 * coefficient * (distance + duration))


def static_potential_from_loop(loop: Any, euclidean_time: sp.Symbol) -> sp.Expr:
    """Extract ``-lim(T->infinity) log(loop)/T`` from a declared loop law."""

    if not isinstance(euclidean_time, sp.Symbol):
        raise TypeError("euclidean_time must be a SymPy Symbol")
    expression = sp.sympify(loop)
    return sp.simplify(
        -sp.limit(sp.log(expression) / euclidean_time, euclidean_time, sp.oo)
    )
