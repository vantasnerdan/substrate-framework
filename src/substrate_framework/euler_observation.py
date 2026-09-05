"""Unpromoted exact material-tag observations for P251's Euler continuation.

Mass is retained: the hybrid momentum replaces each tagged parcel by its
actual centroid momentum and keeps the complete untagged ambient fluid.
The Euler point-filter momentum differs by both an intrinsic-spin dipole
and a symmetric shape-rate dipole. These kinematic identities do not supply
an autonomous microrotation law, a stationary ensemble or PDE regularity.
Fourier convention: exp(-I*k.x). Positions/velocities are simultaneous
material data; velocities must be the actual derivatives of positions.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework.conserved_moments import discrete_mass_moments


def _column(value, name):
    result = sp.Matrix(value)
    if result.shape == (1, 3):
        result = result.T
    if result.shape != (3, 1):
        raise ValueError(f"{name} must have three components")
    return result


def _immutable(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.simplify))


@dataclass(frozen=True)
class MaterialTagMoments:
    """Actual central moments; ``first_momentum[i,j]=sum m*v_i*r_j``."""

    mass: sp.Expr
    centroid: sp.ImmutableMatrix
    momentum: sp.ImmutableMatrix
    second_mass: sp.ImmutableMatrix
    first_momentum: sp.ImmutableMatrix
    spin: sp.ImmutableMatrix
    shape_rate: sp.ImmutableMatrix


def material_tag_moments(masses, positions, velocities):
    """Compute exact finite-mass material moments, not a continuum quadrature.

    Positive finite real masses are required; undecidable symbolic signs
    remain caller hypotheses. Continuous parcels obey the same identities
    with integrals. This finite representation proves only its moment algebra.
    """
    if len(velocities) != len(masses):
        raise ValueError("masses, positions and velocities must have equal length")
    weights = tuple(sp.sympify(m) for m in masses)
    for mass in weights:
        if (mass.is_positive is False or mass.is_real is False
                or mass.is_finite is False
                or (mass.is_number and mass.is_positive is not True)):
            raise ValueError("masses must be positive finite real values")
    raw = discrete_mass_moments(weights, positions)
    if raw.monopole.is_zero is True:
        raise ValueError("total mass must be nonzero")
    points = tuple(_column(x, "position") for x in positions)
    rates = tuple(_column(v, "velocity") for v in velocities)
    centroid = raw.dipole / raw.monopole
    momentum = sum((m*v for m, v in zip(weights, rates)), sp.zeros(3, 1))
    radii = tuple(x-centroid for x in points)
    first = sum((m*v*r.T for m, v, r in zip(weights, rates, radii)), sp.zeros(3))
    spin = sum((m*r.cross(v) for m, r, v in zip(weights, radii, rates)), sp.zeros(3, 1))
    second = raw.second_moment - raw.monopole*centroid*centroid.T
    return MaterialTagMoments(raw.monopole, *map(_immutable, (
        centroid, momentum, second, first, spin, first+first.T)))


def material_tag_fourier_dipole(wavevector, spin, shape_rate):
    """Return the point-minus-centroid momentum first jet, without its phase.

    The result is I*k.cross(S)/2-I*dotI*k/2. Multiply by
    exp(-I*k.X) for a tag at X. Inputs S and dotI are actual simultaneous
    material moments; dotI must be symmetric. The omitted point-momentum
    remainder is bounded by |k|^2*integral rho*|v|*|r|^2/2 for real k.
    No compact-velocity or rigid-tag assumption removes the shape term.
    """
    k = _column(wavevector, "wavevector")
    angular = _column(spin, "spin")
    shape = sp.Matrix(shape_rate)
    if shape.shape != (3, 3) or sp.simplify(shape-shape.T) != sp.zeros(3):
        raise ValueError("shape rate must be a symmetric three by three matrix")
    return _immutable(sp.I*(k.cross(angular)-shape*k)/2)
