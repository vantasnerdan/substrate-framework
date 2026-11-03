"""Exact cross energies for declared massive triplet dipoles in three dimensions.

The APIs in this module start from a static three-component linear field with
positive stiffness ``K`` and operator ``-Delta + m**2``. Each point source is
a spatial dipole in every field component. Isolated point-dipole self energies
are divergent; the returned interaction is only the finite cross term after
subtracting the two isolated self energies.

This conditional long-range model is useful for auditing a Skyrmion dipole
approximation. It does not derive the linear field or its sources from a
nonlinear Skyrme action, construct a two-Skyrmion solution, quantize a nucleon,
describe a short-range core, or establish binding.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _nonnegative(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return expression


def _matrix3(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(value)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3 by 3 matrix")
    return matrix


def _unit_vector3(value: Sequence[Any], name: str) -> sp.ImmutableMatrix:
    vector = sp.ImmutableMatrix(value)
    if vector.shape not in {(3, 1), (1, 3)}:
        raise ValueError(f"{name} must have three components")
    vector = sp.ImmutableMatrix(3, 1, list(vector))
    if sp.simplify(vector.dot(vector) - 1) != 0:
        raise ValueError(f"{name} must be exactly unit normalized")
    return vector


def _proper_rotation3(value: Any) -> sp.ImmutableMatrix:
    rotation = _matrix3(value, "relative_orientation")
    orthogonality = rotation.T * rotation - sp.eye(3)
    if any(sp.simplify(entry) != 0 for entry in orthogonality):
        raise ValueError("relative_orientation must be orthogonal")
    if sp.simplify(rotation.det() - 1) != 0:
        raise ValueError("relative_orientation must have determinant one")
    return rotation


@dataclass(frozen=True)
class YukawaRadialHessian:
    """Radial Green function and Hessian coefficients for ``R>0``."""

    radius: sp.Expr
    mass: sp.Expr
    green: sp.Expr
    transverse: sp.Expr
    longitudinal: sp.Expr
    anisotropic: sp.Expr
    third_radial_derivative: sp.Expr


def yukawa_radial_hessian(radius: Any, mass: Any) -> YukawaRadialHessian:
    r"""Return the exact Hessian data of ``exp(-m*R)/(4*pi*R)``.

    For a unit radial vector ``u``, the Cartesian Hessian is
    ``H = transverse*I + anisotropic*u*u.T``, where
    ``transverse=G'/R`` and ``anisotropic=G''-G'/R``. ``mass`` may be zero,
    giving the massless dipole limit, while ``radius`` must be positive.
    """

    separation = _positive(radius, "radius")
    field_mass = _nonnegative(mass, "mass")
    dummy = sp.Dummy("R", positive=True)
    dummy_green = sp.exp(-field_mass * dummy) / (4 * sp.pi * dummy)
    green = dummy_green.subs(dummy, separation)
    first = sp.diff(dummy_green, dummy).subs(dummy, separation)
    second = sp.diff(dummy_green, dummy, 2).subs(dummy, separation)
    third = sp.diff(dummy_green, dummy, 3).subs(dummy, separation)
    transverse = sp.factor(first / separation)
    longitudinal = sp.factor(second)
    return YukawaRadialHessian(
        radius=separation,
        mass=field_mass,
        green=sp.factor(green),
        transverse=transverse,
        longitudinal=longitudinal,
        anisotropic=sp.factor(longitudinal - transverse),
        third_radial_derivative=sp.factor(third),
    )


@dataclass(frozen=True)
class MassiveTripletDipoleInteraction:
    """One declared relative orientation and its finite cross energy."""

    radius: sp.Expr
    mass: sp.Expr
    stiffness: sp.Expr
    dipole_strength: sp.Expr
    radial_direction: sp.ImmutableMatrix
    relative_orientation: sp.ImmutableMatrix
    orientation_contraction: sp.Expr
    interaction_energy: sp.Expr


def massive_triplet_dipole_interaction(
    radius: Any,
    mass: Any,
    stiffness: Any,
    dipole_strength: Any,
    radial_direction: Sequence[Any],
    relative_orientation: Any,
) -> MassiveTripletDipoleInteraction:
    r"""Return the isolated-self-energy-subtracted dipole cross term.

    The declared static energy is
    ``E[phi;J] = K/2*int((grad phi_a)^2+m^2 phi_a^2) - int J_a phi_a``
    with ``K*(-Delta+m^2)phi_a=J_a``. Two equal-strength triplet dipoles have
    relative proper rotation ``D`` and separation direction ``u``. After the
    two isolated self energies are removed, their on-shell cross energy is
    ``P^2/K * (G'/R*tr(D) + (G''-G'/R)*u.T*D*u)``.

    The source sign and on-shell energy convention are part of this API. A
    different source coupling changes the interaction sign.
    """

    spring = _positive(stiffness, "stiffness")
    strength = _nonnegative(dipole_strength, "dipole_strength")
    direction = _unit_vector3(radial_direction, "radial_direction")
    rotation = _proper_rotation3(relative_orientation)
    hessian = yukawa_radial_hessian(radius, mass)
    radial_projection = (direction.T * rotation * direction)[0]
    contraction = sp.simplify(
        hessian.transverse * sp.trace(rotation)
        + hessian.anisotropic * radial_projection
    )
    energy = sp.factor(strength**2 * contraction / spring)
    return MassiveTripletDipoleInteraction(
        radius=hessian.radius,
        mass=hessian.mass,
        stiffness=spring,
        dipole_strength=strength,
        radial_direction=direction,
        relative_orientation=rotation,
        orientation_contraction=contraction,
        interaction_energy=energy,
    )


@dataclass(frozen=True)
class MassiveTripletDipoleExtrema:
    """Global ``SO(3)`` energy extrema at one positive separation."""

    radius: sp.Expr
    mass: sp.Expr
    stiffness: sp.Expr
    dipole_strength: sp.Expr
    most_attractive_orientation: sp.ImmutableMatrix
    most_repulsive_orientation: sp.ImmutableMatrix
    identity_orientation: sp.ImmutableMatrix
    most_attractive_energy: sp.Expr
    most_repulsive_energy: sp.Expr
    identity_energy: sp.Expr
    most_attractive_radial_force: sp.Expr


def massive_triplet_dipole_extrema(
    radius: Any,
    mass: Any,
    stiffness: Any,
    dipole_strength: Any,
) -> MassiveTripletDipoleExtrema:
    r"""Return exact global orientation extrema and the attractive force.

    Representatives use separation direction ``u=e3``. The global minimum is
    any rotation by pi about an axis perpendicular to ``u``; the representative
    is ``diag(1,-1,-1)``. The global maximum is the rotation by pi about
    ``u``, ``diag(-1,-1,1)``. The force is ``-dE_min/dR`` at fixed relative
    orientation and is strictly negative for positive dipole strength.
    """

    spring = _positive(stiffness, "stiffness")
    strength = _nonnegative(dipole_strength, "dipole_strength")
    hessian = yukawa_radial_hessian(radius, mass)
    scale = strength**2 / spring
    return MassiveTripletDipoleExtrema(
        radius=hessian.radius,
        mass=hessian.mass,
        stiffness=spring,
        dipole_strength=strength,
        most_attractive_orientation=sp.ImmutableMatrix(sp.diag(1, -1, -1)),
        most_repulsive_orientation=sp.ImmutableMatrix(sp.diag(-1, -1, 1)),
        identity_orientation=sp.ImmutableMatrix(sp.eye(3)),
        most_attractive_energy=sp.factor(-scale * hessian.longitudinal),
        most_repulsive_energy=sp.factor(
            scale * (hessian.longitudinal - 2 * hessian.transverse)
        ),
        identity_energy=sp.factor(
            scale * (hessian.longitudinal + 2 * hessian.transverse)
        ),
        most_attractive_radial_force=sp.factor(
            scale * hessian.third_radial_derivative
        ),
    )
