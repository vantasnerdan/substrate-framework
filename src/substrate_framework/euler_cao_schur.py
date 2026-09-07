"""Exact leading parameter and Schur algebra for a thin Cao vortex ring.

The functions in this module only evaluate the source-derived leading jet.
They do not prove the Green-operator Fredholm theorem or construct a charged
Euler--Maxwell branch.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class CaoThinRingSchurJet:
    """Leading thin-core parameter, moment, and Jacobian expressions."""

    chemical_potential: sp.Expr
    translation_speed: sp.Expr
    physical_axial_impulse: sp.Expr
    parameter_jacobian: sp.Expr
    moment_jacobian: sp.Expr
    circulation_response_at_fixed_speed: sp.Expr
    impulse_response_at_fixed_circulation: sp.Expr
    physical_schur_determinant: sp.Expr


def cao_thin_ring_schur_jet(circulation, radius, log_inverse_core, density):
    """Return the exact algebraic leading jet used by P253/0080.

    ``circulation``, ``radius``, ``log_inverse_core``, and ``density`` are
    assumed positive.  The returned formulas use physical axial impulse
    ``I_z = pi*rho*kappa*R**2`` and order the parameter columns as
    ``(mu, c)`` and constraint rows as ``(kappa, I_z)``.
    """

    kappa = sp.sympify(circulation)
    ring_radius = sp.sympify(radius)
    core_log = sp.sympify(log_inverse_core)
    mass_density = sp.sympify(density)
    for name, value in (
        ("circulation", kappa),
        ("radius", ring_radius),
        ("log_inverse_core", core_log),
        ("density", mass_density),
    ):
        if value.is_positive is False or value.is_zero is True:
            raise ValueError(f"{name} must be positive")

    chemical = sp.Rational(3, 8) * kappa * ring_radius * core_log / sp.pi
    speed = kappa * core_log / (4 * sp.pi * ring_radius)
    impulse = sp.pi * mass_density * kappa * ring_radius**2

    # Differentiate with temporary symbols so callers may pass numbers.
    k, r = sp.symbols("k r", positive=True)
    mu_kr = sp.Rational(3, 8) * k * r * core_log / sp.pi
    c_kr = k * core_log / (4 * sp.pi * r)
    impulse_kr = sp.pi * mass_density * k * r**2
    parameter_jacobian = sp.factor(
        sp.Matrix(
            [
                [sp.diff(mu_kr, k), sp.diff(mu_kr, r)],
                [sp.diff(c_kr, k), sp.diff(c_kr, r)],
            ]
        ).det()
    ).subs({k: kappa, r: ring_radius})
    moment_jacobian = sp.factor(
        sp.Matrix(
            [
                [1, 0],
                [sp.diff(impulse_kr, k), sp.diff(impulse_kr, r)],
            ]
        ).det()
    ).subs({k: kappa, r: ring_radius})
    schur = sp.factor(moment_jacobian / parameter_jacobian)
    circulation_response = sp.factor(
        sp.diff(c_kr, r) / parameter_jacobian
    ).subs({k: kappa, r: ring_radius})
    impulse_response = sp.factor(
        sp.diff(impulse_kr, r) / sp.diff(c_kr, r)
    ).subs({k: kappa, r: ring_radius})

    return CaoThinRingSchurJet(
        chemical_potential=sp.factor(chemical),
        translation_speed=sp.factor(speed),
        physical_axial_impulse=sp.factor(impulse),
        parameter_jacobian=sp.factor(parameter_jacobian),
        moment_jacobian=sp.factor(moment_jacobian),
        circulation_response_at_fixed_speed=sp.factor(circulation_response),
        impulse_response_at_fixed_circulation=sp.factor(impulse_response),
        physical_schur_determinant=schur,
    )
