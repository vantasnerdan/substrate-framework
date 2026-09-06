"""Exact column linearization and same-fluid exterior for Euler wave analysis.

The expressions are conditional operators, not a numerical wave constructor
or a certification of the operator/implicit-function proof in P253/0027.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ColumnWaveCoefficients:
    linear_potential: sp.Expr
    quadratic_coefficient: sp.Expr
    bernoulli_label_derivative: sp.Expr


def column_wave_coefficients(angular_momentum, axial_velocity, radius, speed):
    """Return Q and J in Delta_* f+Q f+J f²/2+...=0.

    The base streamfunction has psi0'=r*(W-c). This function requires no
    stagnation of W-c where the labels are inverted. Unresolved symbolic
    nonvanishing remains a hypothesis. Angular momentum is r*u_theta.
    """
    r = radius
    if not isinstance(r, sp.Symbol) or r.is_nonpositive is True:
        raise ValueError("radius must be a symbol on r>0")
    L, W, c = map(sp.sympify, (angular_momentum, axial_velocity, speed))
    q = W-c
    if c.has(r) or q.is_zero is True:
        raise ValueError("speed is constant and W-c must be nonzero")
    def label_derivative(expr):
        return sp.diff(expr, r)/(r*q)
    FFprime = L*sp.diff(L, r)/(r*q)
    Bprime = sp.diff(W, r)/r+FFprime/r**2
    Q = label_derivative(FFprime)-r*r*label_derivative(Bprime)
    J = label_derivative(label_derivative(FFprime))-r*r*label_derivative(label_derivative(Bprime))
    return ColumnWaveCoefficients(sp.simplify(Q), sp.simplify(J), Bprime)


@dataclass(frozen=True)
class ColumnExterior:
    unit_trace_streamfunction: sp.Expr
    outward_radial_derivative: sp.Expr
    energy_trace_symbol: sp.Expr


def column_exterior(wavenumber, matching_radius, radius):
    """Decaying Delta_* solution on r>=R, normalized to boundary value one.

    The outward derivative is with respect to increasing r at R. Exterior
    energy per Fourier mode is energy_trace_symbol*|trace|²/2, using dr/r.
    k=0 denotes the constant-in-r limit, whose velocity energy is zero for
    that single mode; it is not an isolated square-integrable radial state.
    """
    k, R, r = map(sp.sympify, (wavenumber, matching_radius, radius))
    if k.is_real is False or R.is_nonpositive is True:
        raise ValueError("wavenumber must be real and matching radius positive")
    if k.has(r) or R.has(r):
        raise ValueError("wavenumber and matching radius are constant in r")
    if k.is_zero is True:
        return ColumnExterior(sp.S.One, sp.S.Zero, sp.S.Zero)
    q = sp.Abs(k)
    trace = r*sp.besselk(1, q*r)/(R*sp.besselk(1, q*R))
    derivative = -q*sp.besselk(0, q*R)/sp.besselk(1, q*R)
    return ColumnExterior(trace, derivative, -derivative/R)
