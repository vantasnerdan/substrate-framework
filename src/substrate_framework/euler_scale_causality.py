"""Euler similarity weights and the compact-source pressure quadrupole.

These identities expose two boundaries relevant to particle models.  They do
not select an action quantum and they do not turn an effective carrier band
into a relativistic causal cone.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class EulerSimilarityWeights:
    """Multiplicative weights for ``u_AB(x,t)=A*u(B*x,A*B*t)``."""

    length: sp.Expr
    time: sp.Expr
    velocity: sp.Expr
    vorticity: sp.Expr
    pressure: sp.Expr
    mass: sp.Expr
    energy: sp.Expr
    linear_momentum: sp.Expr
    vorticity_impulse: sp.Expr
    angular_momentum: sp.Expr
    kks_action: sp.Expr
    helicity: sp.Expr
    circulation: sp.Expr
    tag_inertia: sp.Expr
    topological_charge: sp.Expr


def euler_similarity_weights(velocity_scale, inverse_length_scale):
    """Return exact fixed-density similarity weights.

    ``inverse_length_scale`` is ``B``: lengths shrink by ``B**-1``.  The
    topological entry records invariance under the orientation-preserving
    dilation; it is not a dimensionful observable.
    """

    a = sp.sympify(velocity_scale)
    b = sp.sympify(inverse_length_scale)
    if a.is_positive is not True or b.is_positive is not True:
        raise ValueError("Euler scales must be positive")
    return EulerSimilarityWeights(
        length=1 / b,
        time=1 / (a * b),
        velocity=a,
        vorticity=a * b,
        pressure=a**2,
        mass=b**-3,
        energy=a**2 / b**3,
        linear_momentum=a / b**3,
        vorticity_impulse=a / b**3,
        angular_momentum=a / b**4,
        kks_action=a / b**4,
        helicity=a**2 / b**2,
        circulation=a / b,
        tag_inertia=b**-5,
        topological_charge=sp.Integer(1),
    )


def pressure_quadrupole(second_moment, position, *, density):
    """Return the leading pressure of a compact quadratic Euler source.

    For ``M_ij=integral u_i*u_j dx`` and

    ``p=rho*(-Delta)^-1 partial_i partial_j(u_i*u_j)``,

    the returned term is ``rho*M_ij*partial_i partial_j(1/(4*pi*|x|))``.
    It is the far-field leading term, not the complete pressure near the
    source.
    """

    moment = sp.ImmutableMatrix(second_moment)
    point = sp.ImmutableMatrix(position)
    rho = sp.sympify(density)
    if moment.shape != (3, 3) or moment != moment.T:
        raise ValueError("second_moment must be a symmetric 3 by 3 matrix")
    if point.shape != (3, 1):
        raise ValueError("position must have three components")
    if rho.is_positive is not True:
        raise ValueError("density must be positive")
    radius_squared = sp.expand((point.T * point)[0])
    if radius_squared.is_zero is True:
        raise ValueError("far-field position must be nonzero")
    numerator = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            numerator += moment[i, j] * (
                3 * point[i] * point[j]
                - (1 if i == j else 0) * radius_squared
            )
    return sp.simplify(
        rho * numerator / (4 * sp.pi * radius_squared ** sp.Rational(5, 2))
    )


def axisymmetric_swirl_pressure_quadrupole(moment, cylindrical_radius, axial, *, density):
    """Specialize the pressure quadrupole to ``M=diag(moment,moment,0)``."""

    m = sp.sympify(moment)
    s = sp.sympify(cylindrical_radius)
    z = sp.sympify(axial)
    return pressure_quadrupole(
        sp.diag(m, m, 0), sp.Matrix([s, 0, z]), density=density
    )
