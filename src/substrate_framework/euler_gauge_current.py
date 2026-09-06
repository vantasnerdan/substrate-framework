"""Conditional U(1) current and Coulomb helpers for an Euler material tag.

This module describes an explicit enlarged Euler--gauge model.  It does not
derive electromagnetism, charge normalization, or a particle from Euler.
"""

import sympy as sp


def transported_charge_current(coupling, tag, velocity):
    """Return ``(rho_q, J)=(g*chi, g*chi*u)`` for a transported tag."""

    vector = sp.ImmutableMatrix(velocity)
    if vector.shape != (3, 1):
        raise ValueError("velocity must have three components")
    return sp.sympify(coupling) * tag, sp.ImmutableMatrix(sp.sympify(coupling) * tag * vector)


def maxwell_speed(permittivity, permeability):
    """Return the principal gauge-field speed ``1/sqrt(epsilon*mu)``."""

    if permittivity <= 0 or permeability <= 0:
        raise ValueError("permittivity and permeability must be positive")
    return sp.sqrt(1 / (sp.sympify(permittivity) * sp.sympify(permeability)))


def coulomb_pair_energy(permittivity, charge_1, charge_2, displacement):
    """Return the leading static Gauss-field cross energy."""

    if permittivity <= 0:
        raise ValueError("permittivity must be positive")
    position = sp.ImmutableMatrix(displacement)
    if position.shape != (3, 1):
        raise ValueError("displacement must have three components")
    radius_squared = sp.simplify(position.dot(position))
    if radius_squared == 0:
        raise ValueError("displacement must be nonzero")
    return sp.simplify(charge_1 * charge_2 / (4 * sp.pi * permittivity * sp.sqrt(radius_squared)))


def coulomb_force(permittivity, charge_1, charge_2, displacement):
    """Return ``-grad_d E`` for the leading Coulomb cross energy."""

    if permittivity <= 0:
        raise ValueError("permittivity must be positive")
    position = sp.ImmutableMatrix(displacement)
    if position.shape != (3, 1):
        raise ValueError("displacement must have three components")
    radius_squared = sp.simplify(position.dot(position))
    if radius_squared == 0:
        raise ValueError("displacement must be nonzero")
    radius = sp.sqrt(radius_squared)
    return sp.ImmutableMatrix(charge_1 * charge_2 * position / (4 * sp.pi * permittivity * radius**3))
