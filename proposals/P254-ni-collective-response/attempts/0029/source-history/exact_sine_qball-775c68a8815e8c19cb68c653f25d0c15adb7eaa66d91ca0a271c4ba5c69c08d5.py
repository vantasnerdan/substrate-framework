"""Conditional exact-sine stationary Q-ball relations.

The declared dimensionless 1+1 profile equation is

``f_xx = sin(f)/2 - omega**2*f``.

The localized profile is represented by its exact first integral and inverse
quadrature.  No elementary closed form, stability theorem, physical charge
map, or identity with the real sine-Gordon breather is asserted.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq


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


def exact_sine_qball_residual(
    profile: Any, coordinate: sp.Symbol, omega: Any
) -> sp.Expr:
    """Return ``f_xx-sin(f)/2+omega**2*f`` for the declared equation."""

    field = sp.sympify(profile)
    frequency = _frequency(omega)
    return (
        sp.diff(field, coordinate, 2)
        - sp.sin(field) / 2
        + frequency**2 * field
    )


def exact_sine_qball_effective_square(
    field_value: Any, omega: Any
) -> sp.Expr:
    """Return ``f_x**2=1-cos(f)-omega**2*f**2`` on the localized orbit."""

    field = sp.sympify(field_value)
    frequency = _frequency(omega)
    return 1 - sp.cos(field) - frequency**2 * field**2


def exact_sine_qball_first_integral_residual(
    profile: Any, coordinate: sp.Symbol, omega: Any
) -> sp.Expr:
    """Return the localized first-integral residual ``f_x**2-G(f)``."""

    field = sp.sympify(profile)
    return (
        sp.diff(field, coordinate) ** 2
        - exact_sine_qball_effective_square(field, omega)
    )


def exact_sine_qball_coordinate_quadrature(
    field_value: Any, peak: Any, omega: Any
) -> sp.Integral:
    """Return the inverse positive-half profile ``x(f)`` as an exact integral."""

    field = sp.sympify(field_value)
    amplitude = sp.sympify(peak)
    frequency = _frequency(omega)
    integration_field = sp.Dummy("u", real=True, positive=True)
    square = exact_sine_qball_effective_square(
        integration_field, frequency
    )
    return sp.Integral(
        1 / sp.sqrt(square),
        (integration_field, field, amplitude),
    )


def exact_sine_qball_charge_quadrature(
    peak: Any, omega: Any
) -> sp.Expr:
    """Return ``4*omega*integral_0^peak f^2/sqrt(G(f)) df`` unevaluated."""

    amplitude = sp.sympify(peak)
    frequency = _frequency(omega)
    integration_field = sp.Dummy("u", real=True, positive=True)
    square = exact_sine_qball_effective_square(
        integration_field, frequency
    )
    return 4 * frequency * sp.Integral(
        integration_field**2 / sp.sqrt(square),
        (integration_field, 0, amplitude),
    )


def exact_sine_qball_scaled_rhs(
    scaled_profile: Any, inverse_width: Any
) -> sp.Expr:
    """Return the exact RHS for ``f=kappa*F(kappa*x)`` divided by ``kappa**3``."""

    profile = sp.sympify(scaled_profile)
    kappa = sp.sympify(inverse_width)
    if kappa.is_number and (
        kappa.is_real is not True or not float(kappa) > 0.0
    ):
        raise ValueError("inverse_width must be real and positive")
    frequency_squared = sp.Rational(1, 2) - kappa**2
    return sp.simplify(
        (
            sp.sin(kappa * profile) / 2
            - frequency_squared * kappa * profile
        )
        / kappa**3
    )


def _effective_square_float(field_value: float, frequency: float) -> float:
    if field_value == 0.0:
        return 0.0
    half = field_value / 2.0
    sinc = math.sin(half) / half
    return field_value**2 * (0.5 * sinc**2 - frequency**2)


def exact_sine_qball_peak_amplitude(omega: Any) -> float:
    """Return the unique nonzero peak in ``(0,2*pi)`` by bracketed root solve."""

    frequency = _frequency(omega)
    if frequency.is_number is not True:
        raise ValueError("omega must be numeric for peak root solving")
    value = float(frequency)
    inverse_width = math.sqrt(0.5 - value**2)
    lower = min(1.0e-6, inverse_width * 1.0e-3)
    return float(
        brentq(
            lambda field: _effective_square_float(field, value),
            lower,
            2.0 * math.pi,
            xtol=5.0e-15,
            rtol=1.0e-14,
        )
    )


@dataclass(frozen=True)
class ExactSineChargeEvidence:
    """Numerical evaluation of the exact charge quadrature."""

    frequency: float
    peak: float
    charge: float
    absolute_error: float
    epsabs: float
    epsrel: float


def evaluate_exact_sine_qball_charge(
    omega: Any,
    *,
    epsabs: float = 1.0e-10,
    epsrel: float = 1.0e-10,
) -> ExactSineChargeEvidence:
    """Evaluate the exact charge quadrature with an endpoint-regularizing map."""

    frequency_expression = _frequency(omega)
    if frequency_expression.is_number is not True:
        raise ValueError("omega must be numeric for charge quadrature")
    frequency = float(frequency_expression)
    if epsabs <= 0.0 or epsrel <= 0.0:
        raise ValueError("quadrature tolerances must be positive")
    peak = exact_sine_qball_peak_amplitude(frequency)
    peak_derivative = math.sin(peak) - 2.0 * frequency**2 * peak
    endpoint_limit = peak**3 / math.sqrt(-peak_derivative * peak / 2.0)

    def transformed_integrand(angle: float) -> float:
        cosine = math.cos(angle)
        if abs(cosine) < 1.0e-7:
            return endpoint_limit
        field = peak * math.sin(angle)
        square = _effective_square_float(field, frequency)
        if square <= 0.0:
            return endpoint_limit
        return field**2 * peak * cosine / math.sqrt(square)

    integral, error = quad(
        transformed_integrand,
        0.0,
        math.pi / 2.0,
        epsabs=epsabs,
        epsrel=epsrel,
        limit=300,
    )
    return ExactSineChargeEvidence(
        frequency=frequency,
        peak=peak,
        charge=4.0 * frequency * integral,
        absolute_error=4.0 * frequency * error,
        epsabs=epsabs,
        epsrel=epsrel,
    )
