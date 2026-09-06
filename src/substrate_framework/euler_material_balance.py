"""Unpromoted exact material balances and an exterior Euler pressure fixture.

P253/0001 derives the continuum transport identities. The finite-mass API
evaluates their algebra from supplied simultaneous positions, velocities and
PHYSICAL pressure gradients; it does not solve the pressure Poisson equation
or certify that arbitrary samples belong to one Euler field. Continuum use
requires the corresponding integrals and actual Euler pressure.

The radial-swirl fixture is exact outside its declared support at the initial
time. It supplies smooth finite-energy initial data when the radial profile is
a smooth bump, not a stationary carrier or an all-time particle interaction.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework.euler_observation import (
    MaterialTagMoments,
    material_tag_moments,
)


def _column(value, name):
    result = sp.Matrix(value)
    if result.shape == (1, 3):
        result = result.T
    if result.shape != (3, 1):
        raise ValueError(f"{name} must have three components")
    return result


def _positive(value, name):
    result = sp.sympify(value)
    if (result.is_positive is False or result.is_real is False
            or result.is_finite is False
            or (result.is_number and result.is_positive is not True)):
        raise ValueError(f"{name} must be positive finite real")
    return result


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.simplify))


@dataclass(frozen=True)
class MaterialTagBalance:
    """Unnormalized central tensors and their exact pressure-driven rates."""

    moments: MaterialTagMoments
    centroid_velocity: sp.ImmutableMatrix
    force: sp.ImmutableMatrix
    centroid_acceleration: sp.ImmutableMatrix
    kinetic_covariance: sp.ImmutableMatrix
    pressure_moment: sp.ImmutableMatrix
    first_momentum_rate: sp.ImmutableMatrix
    second_mass_acceleration: sp.ImmutableMatrix
    spin_rate: sp.ImmutableMatrix
    kinetic_covariance_rate: sp.ImmutableMatrix
    internal_energy: sp.Expr
    internal_energy_rate: sp.Expr
    total_energy_rate: sp.Expr


def material_tag_balance(masses, positions, velocities, pressure_gradients, *, density):
    """Evaluate exact central transport rates with acceleration -grad(p)/rho.

    ``pressure_moment[i,j] = sum mass*grad(p)[i]*r[j]/rho``.
    No zero surface force/torque or rigid-shape hypothesis is imposed.
    Positive symbolic density/mass signs undecidable here remain hypotheses.
    """
    rho = _positive(density, "density")
    if len(pressure_gradients) != len(masses):
        raise ValueError("pressure gradients and masses must have equal length")
    moments = material_tag_moments(masses, positions, velocities)
    weights = tuple(map(sp.sympify, masses))
    points = tuple(_column(x, "position") for x in positions)
    rates = tuple(_column(v, "velocity") for v in velocities)
    gradients = tuple(_column(g, "pressure gradient") for g in pressure_gradients)
    mean_v = moments.momentum / moments.mass
    radii = tuple(x - moments.centroid for x in points)
    relative_v = tuple(v - mean_v for v in rates)
    accelerations = tuple(-g/rho for g in gradients)
    force = sum((m*a for m, a in zip(weights, accelerations)), sp.zeros(3, 1))
    mean_a = force / moments.mass
    kinetic = sum((m*c*c.T for m, c in zip(weights, relative_v)), sp.zeros(3))
    pressure = sum((m*g*r.T/rho for m, g, r in zip(weights, gradients, radii)), sp.zeros(3))
    first_rate = kinetic-pressure
    spin_rate = sum((m*r.cross(a) for m, r, a in zip(weights, radii, accelerations)), sp.zeros(3, 1))
    kinetic_rate = sum((m*((a-mean_a)*c.T+c*(a-mean_a).T)
                        for m, a, c in zip(weights, accelerations, relative_v)), sp.zeros(3))
    internal_rate = sp.simplify(sp.trace(kinetic_rate)/2)
    return MaterialTagBalance(
        moments, *map(_matrix, (mean_v, force, mean_a, kinetic, pressure,
                                first_rate, first_rate+first_rate.T,
                                spin_rate, kinetic_rate)),
        sp.simplify(sp.trace(kinetic)/2), internal_rate,
        sp.simplify(mean_v.dot(force)+internal_rate),
    )


def radial_swirl_exterior_pressure(position, center, axis, radial_moment, *,
                                   density, support_radius):
    """Exact initial pressure outside u=f(s)*(axis cross (x-center)).

    The caller supplies J=integral_0^a s^4*f(s)^2 ds>0, a positive support
    radius a, and a unit axis. For a smooth compact radial f(s)=g(s^2), the
    velocity is smooth, divergence free and finite energy. This API does not
    construct f or verify J from it. Symbolically undecidable exterior-domain
    inequalities remain explicit caller hypotheses; known interior points
    are rejected. Pressure tends to zero at infinity.
    """
    rho = _positive(density, "density")
    moment = _positive(radial_moment, "radial moment")
    radius = _positive(support_radius, "support radius")
    r = _column(position, "position")-_column(center, "center")
    n = _column(axis, "axis")
    if sp.simplify(n.dot(n)-1) != 0 or any(v.is_real is False for v in n):
        raise ValueError("axis must be a real unit vector")
    distance2 = sp.simplify(r.dot(r))
    gap = sp.simplify(distance2-radius**2)
    if gap.is_nonpositive is True:
        raise ValueError("position must lie strictly outside the support radius")
    return sp.simplify(-rho*moment*(3*n.dot(r)**2-distance2)/(3*distance2**sp.Rational(5, 2)))
