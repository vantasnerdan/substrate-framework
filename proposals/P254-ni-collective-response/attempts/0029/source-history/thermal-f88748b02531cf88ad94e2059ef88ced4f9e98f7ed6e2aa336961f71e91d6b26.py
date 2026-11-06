"""Exact two-level gates and conditional coth-activation utilities.

The coth scale and response family below are declared mathematical premises.
They do not identify a physical bath, mode, stochastic escape process, or
operating objective.  A dimensionful conditional rate additionally requires a
supplied positive attempt frequency.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _real_splitting(splitting: Any) -> sp.Expr:
    value = sp.sympify(splitting)
    if value.is_number and value.is_real is not True:
        raise ValueError("splitting must be real and dimensionless")
    return value


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _exact_nonnegative(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be provably nonnegative")
    return expression


def _open_unit_coordinate(value: Any, name: str) -> sp.Expr:
    expression = _exact_positive(value, name)
    if expression.is_number and bool(expression >= 1):
        raise ValueError(f"{name} must lie in the open unit interval")
    return expression


def two_level_upper_occupation(splitting: Any) -> sp.Expr:
    """Return ``P = exp(-x)/(1+exp(-x)) = 1/(1+exp(x))``."""

    x = _real_splitting(splitting)
    return 1 / (1 + sp.exp(x))


def two_level_occupation_variance(splitting: Any) -> sp.Expr:
    """Return the Bernoulli upper-state variance ``P*(1-P)``."""

    probability = two_level_upper_occupation(splitting)
    return sp.simplify(probability * (1 - probability))


def symmetric_two_level_gate(splitting: Any) -> sp.Expr:
    """Return ``W = 2*P*(1-P) = sech(x/2)**2/2``."""

    return sp.simplify(2 * two_level_occupation_variance(splitting))


def declared_coth_effective_scale(
    quantum_energy: Any,
    thermal_energy: Any,
) -> sp.Expr:
    r"""Return the declared scale ``q*coth(q/(2*vartheta))/2``.

    Both inputs are exact positive energies.  The expression tends to ``q/2``
    as ``vartheta`` tends to zero and is asymptotic to ``vartheta`` at high
    thermal energy.  Those identities do not make it a physical noise or bath
    temperature without an independently approved model.
    """

    quantum = _exact_positive(quantum_energy, "quantum_energy")
    thermal = _exact_positive(thermal_energy, "thermal_energy")
    return quantum * sp.coth(quantum / (2 * thermal)) / 2


def coth_gated_response_shape(
    barrier: Any,
    quantum_energy: Any,
    thermal_energy: Any,
    *,
    prefactor_exponent: Any = sp.Rational(1, 2),
) -> sp.Expr:
    r"""Return a declared dimensionless coth-gated activation family.

    With ``Theta=q*coth(q/(2*vartheta))/2`` this is

    ``(E/Theta)**a * exp(-E/Theta) * W(q/vartheta)``.

    The exponent ``a`` is exact and nonnegative.  ``a=1/2`` is the BD2
    capillary-prefactor shape after eliminating line tension and area drive;
    ``a=0`` is an equally dimensionally closed constant-prefactor comparison.
    """

    energy_barrier = _exact_positive(barrier, "barrier")
    exponent = _exact_nonnegative(prefactor_exponent, "prefactor_exponent")
    quantum = _exact_positive(quantum_energy, "quantum_energy")
    thermal = _exact_positive(thermal_energy, "thermal_energy")
    scale = declared_coth_effective_scale(quantum, thermal)
    ratio = energy_barrier / scale
    gate = symmetric_two_level_gate(quantum / thermal)
    return ratio**exponent * sp.exp(-ratio) * gate


def conditional_coth_gated_capillary_rate(
    attempt_frequency: Any,
    barrier: Any,
    quantum_energy: Any,
    thermal_energy: Any,
) -> sp.Expr:
    r"""Return BD2's capillary-reduced conditional rate family.

    For ``E=pi*tau**2/p``, the declared source prefactor satisfies
    ``tau/sqrt(p*Theta)=sqrt(E/(pi*Theta))``.  The returned expression is
    therefore

    ``nu/sqrt(pi) * (E/Theta)**(1/2) * exp(-E/Theta) * W(q/vartheta)``.

    A supplied frequency gives units of inverse time.  The function composes
    premises; it derives no kinetic or stochastic escape law.
    """

    frequency = _exact_positive(attempt_frequency, "attempt_frequency")
    return (
        frequency
        * coth_gated_response_shape(
            barrier,
            quantum_energy,
            thermal_energy,
        )
        / sp.sqrt(sp.pi)
    )


def coth_gated_reduced_shape(
    reduced_coordinate: Any,
    barrier_ratio: Any,
    *,
    prefactor_exponent: Any = sp.Rational(1, 2),
) -> sp.Expr:
    r"""Return the temperature-dependent reduced shape.

    Put ``u=tanh(q/(2*vartheta))`` and ``b=E/q``.  Up to a positive
    temperature-independent factor, the response is
    ``u**a*(1-u**2)*exp(-2*b*u)``.  Numeric ``u`` must lie in ``(0,1)``;
    symbolic callers are responsible for carrying the upper-domain premise.
    """

    coordinate = _open_unit_coordinate(reduced_coordinate, "reduced_coordinate")
    ratio = _exact_positive(barrier_ratio, "barrier_ratio")
    exponent = _exact_nonnegative(prefactor_exponent, "prefactor_exponent")
    return coordinate**exponent * (1 - coordinate**2) * sp.exp(
        -2 * ratio * coordinate
    )


def coth_gated_log_stationarity_residual(
    reduced_coordinate: Any,
    barrier_ratio: Any,
    *,
    prefactor_exponent: Any = sp.Rational(1, 2),
) -> sp.Expr:
    r"""Return ``d log(shape)/du`` for the reduced response.

    The exact residual is ``a/u - 2*u/(1-u**2) - 2*b``.  For ``a>0`` it is
    strictly decreasing on ``0<u<1`` and has one root.  For ``a=0`` it is
    strictly negative, so the response has no finite-temperature maximum.
    """

    coordinate = _open_unit_coordinate(reduced_coordinate, "reduced_coordinate")
    ratio = _exact_positive(barrier_ratio, "barrier_ratio")
    exponent = _exact_nonnegative(prefactor_exponent, "prefactor_exponent")
    return exponent / coordinate - 2 * coordinate / (1 - coordinate**2) - 2 * ratio


def coth_gated_stationary_coordinate_upper_bound(
    prefactor_exponent: Any = sp.Rational(1, 2),
) -> sp.Expr:
    r"""Return ``sqrt(a/(a+2))``, the strict root bound for ``a>0``.

    A positive barrier ratio shifts the unique stationary coordinate below
    the zero-barrier root.  At ``a=1/2`` the bound is ``1/sqrt(5)``.
    """

    exponent = _exact_positive(prefactor_exponent, "prefactor_exponent")
    return sp.sqrt(exponent / (exponent + 2))


def activated_barrier_log_elasticity(
    barrier: Any,
    activation_scale: Any,
    *,
    prefactor_exponent: Any = sp.Rational(1, 2),
) -> sp.Expr:
    r"""Return the log elasticity ``a-E/Theta`` of the response in ``E``."""

    energy_barrier = _exact_positive(barrier, "barrier")
    scale = _exact_positive(activation_scale, "activation_scale")
    exponent = _exact_nonnegative(prefactor_exponent, "prefactor_exponent")
    return exponent - energy_barrier / scale


def inverse_power_input_log_elasticity(
    barrier: Any,
    activation_scale: Any,
    inverse_power: Any,
    *,
    prefactor_exponent: Any = sp.Rational(1, 2),
) -> sp.Expr:
    r"""Return the response elasticity for ``E proportional to z**(-m)``.

    It is ``m*(E/Theta-a)``.  BD2's amplitude and wavenumber each have
    ``m=2`` only while the quantum, gate, scale, attempt frequency, and all
    other inputs are held fixed.
    """

    power = _exact_positive(inverse_power, "inverse_power")
    return -power * activated_barrier_log_elasticity(
        barrier,
        activation_scale,
        prefactor_exponent=prefactor_exponent,
    )
