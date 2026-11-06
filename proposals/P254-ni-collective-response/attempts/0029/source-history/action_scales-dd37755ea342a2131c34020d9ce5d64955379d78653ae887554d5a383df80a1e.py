"""Generic exact relations among energy, frequency, and normalized action."""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) > 0.0:
            raise ValueError(f"{name} must be real and positive")
    return expression


def secant_action_scale(energy: Any, frequency: Any) -> sp.Expr:
    """Return the action-dimensioned secant ``energy/frequency``."""

    return sp.simplify(
        _positive(energy, "energy") / _positive(frequency, "frequency")
    )


def rigid_rotor_energy(inertia: Any, angular_frequency: Any) -> sp.Expr:
    """Return ``I*omega**2/2`` for a positive rigid rotor."""

    moment = _positive(inertia, "inertia")
    frequency = _positive(angular_frequency, "angular_frequency")
    return sp.simplify(moment * frequency**2 / 2)


def rigid_rotor_normalized_action(
    inertia: Any, angular_frequency: Any
) -> sp.Expr:
    """Return ``(1/(2*pi))*closed_integral(p dtheta) = I*omega``."""

    moment = _positive(inertia, "inertia")
    frequency = _positive(angular_frequency, "angular_frequency")
    return sp.simplify(moment * frequency)
