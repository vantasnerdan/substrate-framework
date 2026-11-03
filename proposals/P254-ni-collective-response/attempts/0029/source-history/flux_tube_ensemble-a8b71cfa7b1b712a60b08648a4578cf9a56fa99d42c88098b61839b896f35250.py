"""Exact conditional stress atoms for an isotropic static-tube ensemble.

Authority status: conditional, unpromoted infrastructure linked to open goal
issue #28 (vantasnerdan/substrate-framework). No accepted claim backs these
symbols.

These functions establish only the declared inverse-quartic static-tube model:
transverse pressure equals the density, longitudinal pressure is its negative,
and an isotropic spatial orientation average has ``w = 1/3``. They make no
claim about a Lorentz-invariant ensemble or completed vacuum.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def tube_transverse_pressure(coefficient: Any) -> sp.Expr:
    """Return ``p_perp`` for ``rho(L) = -coefficient * L**-4`` and ``A=L**2``.

    The pressure is computed from ``-d(rho*A)/dA`` and equals ``rho`` in this
    declared model.
    """

    c = sp.sympify(coefficient)
    length = sp.Symbol("L", positive=True)
    density = -c * length**-4
    area = length**2
    energy_per_length = density * area
    return -sp.simplify(
        sp.diff(energy_per_length, length) / sp.diff(area, length)
    )


def tube_longitudinal_pressure(coefficient: Any) -> sp.Expr:
    """Return ``p_parallel = -rho`` for the declared static uniform tube."""

    c = sp.sympify(coefficient)
    length = sp.Symbol("L", positive=True)
    density = -c * length**-4
    return sp.simplify(-density)


def isotropic_average_equation_of_state(coefficient: Any) -> sp.Expr:
    """Return ``w = 1/3`` for a nonzero-density isotropic static-tube average."""

    c = sp.sympify(coefficient)
    if c.is_zero is True:
        raise ValueError("coefficient must be nonzero so p/rho is defined")
    length = sp.Symbol("L", positive=True)
    density = -c * length**-4
    transverse_pressure = tube_transverse_pressure(c)
    longitudinal_pressure = tube_longitudinal_pressure(c)
    average_pressure = sp.simplify(
        (2 * transverse_pressure + longitudinal_pressure) / 3
    )
    return sp.simplify(average_pressure / density)


def orientation_average_nn(i: int, j: int) -> sp.Expr:
    """Return the exact unit-sphere average ``<n_i n_j> = delta_ij/3``."""

    if i not in range(3) or j not in range(3):
        raise ValueError(f"orientation indices must be in range(3); got {(i, j)!r}")
    theta, phi = sp.symbols("theta phi", positive=True)
    direction = (
        sp.sin(theta) * sp.cos(phi),
        sp.sin(theta) * sp.sin(phi),
        sp.cos(theta),
    )
    measure = sp.sin(theta)
    sphere_volume = sp.integrate(
        measure,
        (theta, 0, sp.pi),
        (phi, 0, 2 * sp.pi),
    )
    moment = sp.integrate(
        direction[i] * direction[j] * measure,
        (theta, 0, sp.pi),
        (phi, 0, 2 * sp.pi),
    )
    return sp.simplify(moment / sphere_volume)
