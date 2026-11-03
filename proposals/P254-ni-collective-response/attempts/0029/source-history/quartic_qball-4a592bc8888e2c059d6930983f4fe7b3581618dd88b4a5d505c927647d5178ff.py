"""Exact conditional profile relations for the declared quartic Q-ball ODE.

The accepted equation is a dimensionless 1+1 stationary-profile model,

``f_xx = (1/2 - omega**2 - f**2/12) f``.

This module does not derive that equation from an action and does not assert
spectral, orbital, or nonlinear stability.  A sign of ``dQ/domega`` is exposed
as calculus only; applying a Vakhitov-Kolokolov theorem requires separately
verified hypotheses and a fluctuation operator.
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from .u1_charge import sech_profile_u1_charge


def _frequency(omega: Any) -> sp.Expr:
    value = sp.sympify(omega)
    upper = sp.sqrt(sp.Rational(1, 2))
    if value.is_number:
        if (
            value.is_real is not True
            or not 0.0 < float(value) < float(upper)
        ):
            raise ValueError(
                "omega must be real and satisfy 0 < omega < 1/sqrt(2)"
            )
    return value


def quartic_qball_inverse_width(omega: Any) -> sp.Expr:
    """Return ``kappa=sqrt(1/2-omega**2)`` on the localized branch."""

    frequency = _frequency(omega)
    return sp.sqrt(sp.Rational(1, 2) - frequency**2)


def quartic_qball_amplitude(omega: Any) -> sp.Expr:
    """Return the positive-profile amplitude ``sqrt(24)*kappa``."""

    return sp.simplify(
        sp.sqrt(24) * quartic_qball_inverse_width(omega)
    )


def quartic_qball_profile(
    coordinate: Any, omega: Any, center: Any = 0
) -> sp.Expr:
    """Return the positive translated sech profile of the declared ODE."""

    x = sp.sympify(coordinate)
    origin = sp.sympify(center)
    kappa = quartic_qball_inverse_width(omega)
    amplitude = sp.sqrt(24) * kappa
    return amplitude * sp.sech(kappa * (x - origin))


def quartic_qball_residual(
    profile: Any, coordinate: sp.Symbol, omega: Any
) -> sp.Expr:
    """Return the fixed quartic ODE residual ``f_xx-(1/2-w^2-f^2/12)f``."""

    field = sp.sympify(profile)
    frequency = _frequency(omega)
    return sp.diff(field, coordinate, 2) - (
        sp.Rational(1, 2) - frequency**2 - field**2 / 12
    ) * field


def quartic_qball_charge(omega: Any) -> sp.Expr:
    """Return ``Q=96*omega*sqrt(1/2-omega**2)`` in C-U1-001's convention."""

    frequency = _frequency(omega)
    kappa = quartic_qball_inverse_width(frequency)
    amplitude = quartic_qball_amplitude(frequency)
    return sp.simplify(
        sech_profile_u1_charge(frequency, kappa, amplitude)
    )


def quartic_qball_charge_derivative(omega: Any) -> sp.Expr:
    """Return the exact charge-family slope, without a stability interpretation."""

    frequency = _frequency(omega)
    symbol = sp.Dummy("omega", positive=True)
    derivative = sp.diff(
        96 * symbol * sp.sqrt(sp.Rational(1, 2) - symbol**2),
        symbol,
    )
    return sp.simplify(derivative.subs(symbol, frequency))
