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


@dataclass(frozen=True)
class ColumnLinearDynamics:
    vorticity_rhs: sp.Expr
    angular_momentum_rhs: sp.Expr
    casimir_second_derivative: sp.Expr
    swirl_energy_weight: sp.Expr


def column_linear_dynamics(angular_momentum, poloidal_streamfunction,
                           swirl_perturbation, radius, axial_coordinate):
    """Axisymmetric Euler linearization at the pure-swirl column L(r)/r.

    eta=delta(omega_theta/r), chi=delta(r*u_theta), and psi=K eta.
    The returned RHSs retain the whole-space pressure through this Hodge
    relation. The normalized conserved quadratic density is
    (eta*psi+swirl_energy_weight*chi**2)/2, integrated with r dr dz;
    its physical multiplier is 2*pi*rho. The weight applies on L'>0;
    the weighted domain at a flat edge is part of the caller's hypotheses.
    Positivity additionally uses L>0. Nonlinear wave stability is not implied.
    """
    r, z = radius, axial_coordinate
    if not isinstance(r, sp.Symbol) or r.is_nonpositive is True:
        raise ValueError("radius must be a symbol on r>0")
    if not isinstance(z, sp.Symbol) or z == r:
        raise ValueError("axial coordinate must be a distinct symbol")
    L, psi, chi = map(sp.sympify, (angular_momentum,
                                 poloidal_streamfunction, swirl_perturbation))
    Lr = sp.diff(L, r)
    if L.has(z) or Lr.is_nonpositive is True:
        raise ValueError("column L depends only on r and needs L'>0")
    Csecond = -sp.diff(L/r**2, r)/Lr
    return ColumnLinearDynamics(
        2*L*sp.diff(chi, z)/r**4,
        Lr*sp.diff(psi, z)/r,
        sp.simplify(Csecond),
        sp.simplify(1/r**2+Csecond),
    )
