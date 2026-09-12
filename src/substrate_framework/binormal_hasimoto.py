"""Exact binormal-flow/NLS data; no finite-core Euler transfer is presumed.

Geometric time is normalized by gamma_t = gamma_s cross gamma_ss. The NLS
Hamiltonian below is a curvature functional, not the Euler kinetic energy.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class BinormalSoliton:
    curve: sp.ImmutableMatrix
    hasimoto_field: sp.Expr


def binormal_soliton(arclength, time, *, eta, xi):
    """Return the exact traveling curve and its focusing Hasimoto field.

    eta must be positive, xi real, both constant. For xi**2 > eta**2 the
    axial coordinate is strictly increasing and the curve is embedded.
    Undecidable symbolic signs remain caller hypotheses. The curve has an
    infinite straight background; its localized displacement has finite
    excess arclength, not finite total length.
    """
    s, t = arclength, time
    if not isinstance(s, sp.Symbol) or not isinstance(t, sp.Symbol) or s == t:
        raise ValueError("arclength and time must be distinct symbols")
    e, x = sp.sympify(eta), sp.sympify(xi)
    if e.is_positive is False or e.is_finite is False:
        raise ValueError("eta must be positive and finite")
    if x.is_real is False or x.is_finite is False or e.has(s, t) or x.has(s, t):
        raise ValueError("eta and xi must be finite real constants")
    y = e*(s-2*x*t)
    theta = x*s+(e*e-x*x)*t
    amplitude = 2*e/(e*e+x*x)
    curve = sp.ImmutableMatrix([
        s-amplitude*sp.tanh(y),
        amplitude*sp.sech(y)*sp.cos(theta),
        amplitude*sp.sech(y)*sp.sin(theta),
    ])
    return BinormalSoliton(curve, 2*e*sp.sech(y)*sp.exp(sp.I*theta))


def nls_densities(field, coordinate):
    """Return mass, momentum and H densities for i psi_t+psi_ss+|psi|²psi/2=0."""
    q = sp.sympify(field)
    derivative = sp.diff(q, coordinate)
    mass = sp.conjugate(q)*q
    momentum = (sp.conjugate(q)*derivative-q*sp.conjugate(derivative))/(2*sp.I)
    energy = sp.conjugate(derivative)*derivative-mass**2/4
    return mass, momentum, energy


def nls_deficit_square(field, coordinate, *, cumulative_mass, total_mass):
    """Return |q_s+(F-M/2)q/4|² for a field with zero total momentum.

    With F_s=|q|², F(-infinity)=0 and F(infinity)=M, its integral is
    H+M**3/192. For a general field first apply q=exp(-i P s/M)psi.
    This function computes the expression and does not certify the integral
    hypotheses, positivity of a density chart, or Euler conservation.
    """
    q = sp.sympify(field)
    residual = sp.diff(q, coordinate)+(sp.sympify(cumulative_mass)
                                           -sp.sympify(total_mass)/2)*q/4
    return sp.conjugate(residual)*residual
