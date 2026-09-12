"""Volume-preserving column bending, before Hodge velocity recovery.

These are exact vorticity initial data. The API does not equate the supplied
curve with an Euler material trajectory or a persistent finite-core soliton.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class BentColumn:
    material_map: sp.ImmutableMatrix
    deformation_gradient: sp.ImmutableMatrix
    vorticity: sp.ImmutableMatrix
    vorticity_difference: sp.ImmutableMatrix


def bent_column(coordinates, center, radial_profile, radial_square):
    """Push w(x²+y²)e_z through (x,y,z)->(x+a(z),y+b(z),z).

    radial_profile is a supplied function of radial_square. Smooth compact
    transverse support, decay of a,b and their derivatives, and convergence
    of all desired integrals are caller hypotheses. Returned vorticity is
    Eulerian; the map and its gradient use material coordinates with the
    same symbolic names. A decaying Hodge inverse of the difference supplies
    the velocity perturbation; no pointwise velocity pushforward is asserted.
    """
    x, y, z = coordinates
    if len({x, y, z, radial_square}) != 4 or not all(
        isinstance(q, sp.Symbol) for q in (x, y, z, radial_square)
    ):
        raise ValueError("coordinates and radial square must be distinct symbols")
    a, b = map(sp.sympify, center)
    profile = sp.sympify(radial_profile)
    if a.has(x, y) or b.has(x, y) or profile.has(x, y, z):
        raise ValueError("center depends only on z; profile only on radial-square variable")
    material_map = sp.ImmutableMatrix([x+a, y+b, z])
    gradient = material_map.jacobian((x, y, z))
    w = profile.subs(radial_square, (x-a)**2+(y-b)**2)
    vorticity = sp.ImmutableMatrix([w*sp.diff(a, z), w*sp.diff(b, z), w])
    baseline = sp.ImmutableMatrix([0, 0, profile.subs(radial_square, x*x+y*y)])
    return BentColumn(material_map, sp.ImmutableMatrix(gradient), vorticity,
                      sp.ImmutableMatrix(vorticity-baseline))


def bent_column_relative_impulse_density(center, axial_coordinate, circulation):
    """Transversely integrated (x cross delta_omega)/2 per unit axial length.

    Requires a radial integrable profile of circulation Gamma, with zero
    transverse first moments. Axial integration and decay are not inferred.
    No mass density is included; physical impulse is rho times the result.
    """
    a, b = map(sp.sympify, center)
    z, g = axial_coordinate, sp.sympify(circulation)
    return sp.ImmutableMatrix([
        g*(b-z*sp.diff(b, z))/2,
        g*(z*sp.diff(a, z)-a)/2,
        g*(a*sp.diff(b, z)-b*sp.diff(a, z))/2,
    ])
