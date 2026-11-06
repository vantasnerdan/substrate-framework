"""Conditional global-U(1) current and declared-profile charge relations.

This module does not complexify the accepted real sine-Gordon field. Its
field-theory functions apply to an independently declared complex scalar on
1+1 Minkowski spacetime with signature ``(+,-)``. The breather-parameterized
helpers additionally declare the profile width to be
``sqrt(1-omega**2)`` so they can be composed with accepted sine-Gordon
quantities without claiming that the complex profile solves sine-Gordon.
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from .sine_gordon import (
    breather_energy,
    breather_inverse_width,
    breather_secant_action_scale,
)


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) > 0.0:
            raise ValueError(f"{name} must be real and positive")
    return expression


def minkowski_dalembertian(
    field: Any, coordinate: sp.Symbol, time: sp.Symbol
) -> sp.Expr:
    """Return ``Box(field)=field_tt-field_xx`` for signature ``(+,-)``."""

    expression = sp.sympify(field)
    return sp.diff(expression, time, 2) - sp.diff(expression, coordinate, 2)


def u1_current_components(
    field: Any,
    conjugate_field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the raised components ``(j^0,j^1)`` of the declared U(1) current.

    The convention is
    ``j^mu=i*(Psi_conj*d^mu(Psi)-Psi*d^mu(Psi_conj))``. Consequently the
    spatial component contains the minus sign from ``d^1=-d_x``.
    ``conjugate_field`` is explicit so symbolic callers can declare the
    conjugacy relation without relying on SymPy assumptions.
    """

    psi = sp.sympify(field)
    psi_conjugate = sp.sympify(conjugate_field)
    density = sp.I * (
        psi_conjugate * sp.diff(psi, time)
        - psi * sp.diff(psi_conjugate, time)
    )
    flux = -sp.I * (
        psi_conjugate * sp.diff(psi, coordinate)
        - psi * sp.diff(psi_conjugate, coordinate)
    )
    return sp.simplify(density), sp.simplify(flux)


def u1_current_divergence(
    field: Any,
    conjugate_field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
) -> sp.Expr:
    """Return ``d_t(j^0)+d_x(j^1)`` for the declared current convention."""

    density, flux = u1_current_components(
        field, conjugate_field, coordinate, time
    )
    return sp.simplify(
        sp.diff(density, time) + sp.diff(flux, coordinate)
    )


def stationary_phase_field(profile: Any, time: Any, frequency: Any) -> sp.Expr:
    """Return the declared stationary-phase ansatz ``profile*exp(-i*omega*t)``."""

    amplitude_profile = sp.sympify(profile)
    temporal_coordinate = sp.sympify(time)
    angular_frequency = _positive(frequency, "frequency")
    return amplitude_profile * sp.exp(-sp.I * angular_frequency * temporal_coordinate)


def stationary_u1_charge_density(profile: Any, frequency: Any) -> sp.Expr:
    """Return ``2*omega*profile**2`` for a declared real stationary profile."""

    amplitude_profile = sp.sympify(profile)
    angular_frequency = _positive(frequency, "frequency")
    return sp.simplify(2 * angular_frequency * amplitude_profile**2)


def sech_profile_u1_charge(
    frequency: Any, inverse_width: Any, amplitude: Any = 1
) -> sp.Expr:
    """Return the exact charge of ``A*sech(eta*x)*exp(-i*omega*t)``.

    The result ``4*A**2*omega/eta`` is conditional on the declared profile;
    this function does not assert an equation of motion that produces it.
    """

    angular_frequency = _positive(frequency, "frequency")
    eta = _positive(inverse_width, "inverse_width")
    scale = _positive(amplitude, "amplitude")
    return sp.simplify(4 * scale**2 * angular_frequency / eta)


def breather_parameterized_u1_charge(
    frequency: Any, amplitude: Any = 1
) -> sp.Expr:
    """Return the declared sech charge with ``eta=sqrt(1-omega**2)``.

    The shared parameterization permits exact composition with accepted
    breather quantities but is not an ontology map between the complex field
    and the accepted real sine-Gordon breather.
    """

    eta = breather_inverse_width(frequency)
    return sech_profile_u1_charge(frequency, eta, amplitude)


def breather_charge_energy_product(
    frequency: Any, amplitude: Any = 1
) -> sp.Expr:
    """Return ``Q*E=64*A**2*omega`` under the declared shared parameterization."""

    return sp.simplify(
        breather_parameterized_u1_charge(frequency, amplitude)
        * breather_energy(frequency)
    )


def breather_charge_secant_product(
    frequency: Any, amplitude: Any = 1
) -> sp.Expr:
    """Return ``Q*(E/omega)=64*A**2`` under the declared parameterization."""

    return sp.simplify(
        breather_parameterized_u1_charge(frequency, amplitude)
        * breather_secant_action_scale(frequency)
    )


def charge_scale_exponent_matrix() -> sp.Matrix:
    """Return frequency/width exponents for primitives ``(Q,H,E,omega)``.

    Rows correspond to powers of ``omega`` and ``eta``. Under the declared
    forms ``Q~omega/eta``, ``H~eta/omega``, and ``E~eta``, its kernel is the
    complete monomial family independent of both ``omega`` and ``eta``.
    """

    return sp.Matrix([[1, -1, 0, 1], [-1, 1, 1, 0]])
