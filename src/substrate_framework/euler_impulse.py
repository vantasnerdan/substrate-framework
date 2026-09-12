"""Conditional, unpromoted compact-vorticity far-field coefficients (P253/0007).

For smooth compact divergence-free vorticity, I=(1/2) integral x cross omega.
The decaying Hodge velocity has the dipole coefficient returned here. The
energy is the leading cross term of the full Euler kinetic energy, not an
assumed particle potential, Coulomb law or autonomous modulation equation.
The caller supplies actual impulses and checks support separation and moment
remainder bounds; algebraic evaluation alone supplies no carrier or stability.
"""

import sympy as sp


def _vector(value):
    vector = sp.Matrix(value)
    if vector.shape == (1, 3):
        vector = vector.T
    if vector.shape != (3, 1) or any(x.is_real is False for x in vector):
        raise ValueError("vectors must have three real components")
    return vector


def impulse_dipole_velocity(impulse, separation):
    """Return [3*r*(r.I)-|r|^2*I]/(4*pi*|r|^5).

    This is the leading velocity for nonzero impulse, with absolute remainder
    O(integral |y|^2|omega|/|r|^4) when |r| exceeds twice the core radius.
    A zero impulse gives zero leading coefficient, not zero exact velocity.
    Symbolically undecidable nonzero separation remains a caller hypothesis.
    """
    moment, r = _vector(impulse), _vector(separation)
    distance2 = sp.simplify(r.dot(r))
    if distance2.is_zero is True:
        raise ValueError("separation must be nonzero")
    result = (3*r*r.dot(moment)-distance2*moment)/(4*sp.pi*distance2**sp.Rational(5, 2))
    return sp.ImmutableMatrix(result.applyfunc(sp.simplify))


def impulse_dipole_cross_energy(impulse_a, impulse_b, separation, *, density):
    """Leading rho*integral u_a.u_b, with positive kinetic-energy convention.

    Equals rho*I_b.u_a(d) at dipole order. A mechanical force additionally
    requires a derived collective action; negative separation differentiation
    by itself does not establish it. Known nonpositive/nonfinite densities
    are rejected; unknown symbolic positivity remains a hypothesis.
    """
    rho = sp.sympify(density)
    if (rho.is_positive is False or rho.is_real is False or rho.is_finite is False
            or (rho.is_number and rho.is_positive is not True)):
        raise ValueError("density must be positive finite real")
    return sp.simplify(rho*_vector(impulse_b).dot(impulse_dipole_velocity(impulse_a, separation)))
