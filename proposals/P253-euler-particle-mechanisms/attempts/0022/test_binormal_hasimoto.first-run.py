"""Exact curve, dynamics, conserved functional and normalization checks."""

import pytest
import sympy as sp

from substrate_framework.binormal_hasimoto import (
    binormal_soliton, nls_deficit_square, nls_densities,
)


def reduce_hyperbolic(expr):
    return sp.factor(sp.cancel(sp.trigsimp(sp.expand_trig(expr)).rewrite(sp.exp)))


def test_full_curve_binormal_equation_metric_and_torsion():
    s, t = sp.symbols("s t", real=True)
    eta = sp.Symbol("eta", positive=True)
    xi = sp.Symbol("xi", real=True)
    y, theta = sp.symbols("y theta", real=True)
    data = binormal_soliton(s, t, eta=eta, xi=xi)
    # Independent jet coordinates keep the complete two-variable identity small.
    curve = data.curve.subs({eta*(s-2*xi*t): y,
                            xi*s+(eta**2-xi**2)*t: theta})
    def ds(v):
        return v.diff(s)+eta*v.diff(y)+xi*v.diff(theta)
    tangent = ds(curve)
    second = ds(tangent)
    velocity = -2*eta*xi*curve.diff(y)+(eta**2-xi**2)*curve.diff(theta)
    residual = velocity-tangent.cross(second)
    assert all(reduce_hyperbolic(x) == 0 for x in residual)
    assert reduce_hyperbolic(tangent.dot(tangent)-1) == 0
    assert reduce_hyperbolic(second.dot(second)-4*eta**2*sp.sech(y)**2) == 0
    cross = tangent.cross(second)
    assert reduce_hyperbolic(cross.dot(ds(second))-xi*cross.dot(cross)) == 0
    # Opposite geometric-flow orientation must be detected.
    assert velocity.subs({y: 0, theta: 0, eta: 1, xi: 2}) != sp.zeros(3, 1)


def test_actual_hasimoto_pde_and_independently_integrated_invariants():
    s, t = sp.symbols("s t", real=True)
    eta = sp.Symbol("eta", positive=True)
    xi = sp.Symbol("xi", real=True)
    psi = binormal_soliton(s, t, eta=eta, xi=xi).hasimoto_field
    residual = sp.I*psi.diff(t)+psi.diff(s, 2)+sp.conjugate(psi)*psi**2/2
    assert reduce_hyperbolic(residual/psi) == 0
    mass, momentum, energy = nls_densities(psi.subs(t, 0), s)
    # z=tanh(eta*s) integrates the actual densities, not supplied totals.
    z = sp.Symbol("z", real=True)
    jacobian = 1/(eta*(1-z*z))
    def integrate_density(density):
        density = sp.simplify(sp.expand_complex(density))
        rational = sp.simplify(density.subs(s, sp.atanh(z)/eta)*jacobian)
        return sp.simplify(sp.integrate(rational, (z, -1, 1)))
    assert integrate_density(mass) == 8*eta
    assert integrate_density(momentum) == 8*eta*xi
    assert integrate_density(energy) == 8*eta*xi**2-sp.Rational(8, 3)*eta**3


def test_square_identity_quantile_factor_and_wrong_coefficient():
    s = sp.Symbol("s", real=True)
    M = sp.Symbol("M", positive=True)
    amplitude = sp.Function("a", real=True)(s)
    phase = sp.Function("theta", real=True)(s)
    F = sp.Function("F", real=True)(s)
    q = amplitude*sp.exp(sp.I*phase)
    square = nls_deficit_square(q, s, cumulative_mass=F, total_mass=M)
    _, _, energy = nls_densities(q, s)
    W = F-M/2
    boundary = sp.diff(W*amplitude**2/4, s)
    # Equality as local densities modulo the explicit boundary derivative.
    expected = energy+W**2*amplitude**2/16+boundary
    assert sp.simplify((square-expected).subs(sp.diff(F, s), amplitude**2)) == 0
    x = sp.Symbol("F", real=True)
    a, a_F, theta_F = sp.symbols("a a_F theta_F", real=True, positive=True)
    quantile = ((a_F/2+(x-M/2)/4)**2+a**2*theta_F**2)
    ground = x*(M-x)/4
    assert sp.expand(quantile-(a_F-sp.diff(ground, x))**2/4-a**2*theta_F**2) == 0
    assert sp.integrate((x-M/2)**2/16, (x, 0, M)) == M**3/192
    # A changed square coefficient cannot preserve the quartic-energy row.
    wrong = sp.diff(q, s)+W*q/2
    error = sp.expand_complex(sp.conjugate(wrong)*wrong-square)
    assert sp.simplify(error.subs({amplitude: 1, phase: 0, F: M/2+1}).doit()) != 0


def test_ground_state_square_and_embedded_excess_length():
    s, t = sp.symbols("s t", real=True)
    eta = sp.Symbol("eta", positive=True)
    xi = sp.Symbol("xi", real=True)
    q = 2*eta*sp.sech(eta*s)
    F = 4*eta*(1+sp.tanh(eta*s))
    assert reduce_hyperbolic(F.diff(s)-q*q) == 0
    square = nls_deficit_square(q, s, cumulative_mass=F, total_mass=8*eta)
    assert reduce_hyperbolic(square) == 0
    axial = binormal_soliton(s, t, eta=eta, xi=xi).curve[0].subs(t, 0)
    primitive = s-axial
    excess = sp.limit(primitive, s, sp.oo)-sp.limit(primitive, s, -sp.oo)
    assert excess == 4*eta/(eta**2+xi**2)
    with pytest.raises(ValueError):
        binormal_soliton(s, t, eta=0, xi=2)
