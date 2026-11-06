"""Conditional algebraic relations used by candidate Skyrme mass formulas.

The functions encode premises for governed composition. They do not establish
the physical validity or empirical accuracy of either mass formula.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def conditional_topological_mass(
    hedgehog_coefficient: Any,
    electron_rest_energy: Any,
) -> sp.Expr:
    """Return the conditional premise ``48*pi**3*B1*m_e*c**2``."""

    coefficient = _positive(hedgehog_coefficient, "hedgehog_coefficient")
    rest_energy = _positive(electron_rest_energy, "electron_rest_energy")
    return 48 * sp.pi**3 * coefficient * rest_energy


def conditional_anw_mass(
    hedgehog_coefficient: Any,
    pion_scale: Any,
    skyrme_coupling: Any,
) -> sp.Expr:
    """Return the conditional premise ``3*pi**2*B1*F_pi/e``."""

    coefficient = _positive(hedgehog_coefficient, "hedgehog_coefficient")
    scale = _positive(pion_scale, "pion_scale")
    coupling = _positive(skyrme_coupling, "skyrme_coupling")
    return 3 * sp.pi**2 * coefficient * scale / coupling


def matched_pion_coupling_ratio(electron_rest_energy: Any) -> sp.Expr:
    """Return ``F_pi/e = 16*pi*m_e*c**2`` implied by matching the premises."""

    rest_energy = _positive(electron_rest_energy, "electron_rest_energy")
    return 16 * sp.pi * rest_energy
