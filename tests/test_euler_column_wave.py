"""Physical cylindrical residual and actual exterior/schur checks."""

import sympy as sp

from substrate_framework.euler_column_wave import (
    column_exterior, column_wave_coefficients, column_linear_dynamics,
)
from substrate_framework.euler_compact_ring import compact_ring_fields


def test_full_steady_euler_residual_is_actual_grad_shafranov_equation():
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    psi = sp.Function("psi")(r, z)
    F = sp.Function("F")(psi)
    B = sp.Function("B")(psi)
    pressure = B-(psi.diff(r)**2+psi.diff(z)**2+F*F)/(2*r*r)
    data = compact_ring_fields(psi, F, pressure, r, z)
    GS = psi.diff(r, 2)-psi.diff(r)/r+psi.diff(z, 2)+F*sp.diff(F, psi)-r*r*sp.diff(B, psi)
    expected = -GS*sp.Matrix([psi.diff(r), 0, psi.diff(z)])/(r*r)
    assert (data.residual-expected).applyfunc(sp.simplify) == sp.zeros(3, 1)


def test_general_column_and_pure_swirl_quadratic_coefficient():
    r = sp.Symbol("r", positive=True)
    c = sp.Symbol("c", nonzero=True, real=True)
    L, W = sp.Function("L")(r), sp.Function("W")(r)
    data = column_wave_coefficients(L, W, r, c)
    q = W-c
    expected = 2*L*L.diff(r)/(r**3*q*q)-(W.diff(r, 2)-W.diff(r)/r)/q
    assert sp.simplify(data.linear_potential-expected) == 0
    pure = column_wave_coefficients(L, 0, r, c)
    assert sp.simplify(pure.quadratic_coefficient+2*pure.linear_potential.diff(r)/(r*c)) == 0
    uniform = column_wave_coefficients(r*r, 0, r, c)
    assert uniform.quadratic_coefficient == 0
    assert uniform.linear_potential == 4/c**2


def test_exterior_dtn_and_actual_weighted_kinetic_energy_flux():
    r, R, k = sp.symbols("r R k", positive=True)
    data = column_exterior(k, R, r)
    f = data.unit_trace_streamfunction
    assert sp.simplify(sp.expand_func(f.diff(r, 2)-f.diff(r)/r-k*k*f)) == 0
    assert f.subs(r, R) == 1
    assert sp.simplify(sp.expand_func(f.diff(r).subs(r, R)-data.outward_radial_derivative)) == 0
    assert sp.simplify(sp.expand_func(sp.diff(f*f.diff(r)/r, r)-(f.diff(r)**2+k*k*f*f)/r)) == 0
    assert sp.simplify(data.energy_trace_symbol+data.outward_radial_derivative/R) == 0
    assert column_exterior(-k, R, r) == data
    assert column_exterior(0, R, r).unit_trace_streamfunction == 1


def test_exact_rank_one_inverse_and_limiting_homoclinic_zero_mode():
    mu, k, T, a, b = sp.symbols("mu k T a b", positive=True)
    h, q, p = sp.symbols("h q p", real=True)
    boundary = sp.Matrix([h, q, p])
    full = sp.diag(mu+k*k, a+k*k, b+k*k)+T*boundary*boundary.T
    Rq = sp.diag(0, 1/(a+k*k), 1/(b+k*k))
    v = Rq*boundary
    d = (boundary.T*v)[0]
    e = sp.Matrix([1, 0, 0])-T*h*v/(1+T*d)
    m = mu+k*k+T*h*h/(1+T*d)
    inverse = Rq-T*v*v.T/(1+T*d)+e*e.T/m
    assert (full*inverse-sp.eye(3)).applyfunc(sp.factor) == sp.zeros(3)
    X = sp.Symbol("X", real=True)
    beta = sp.Symbol("beta", positive=True)
    A = 3*sp.sech(X/2)**2/(2*beta)
    def reduce(expr):
        return sp.factor(sp.cancel(expr.rewrite(sp.exp)))
    assert reduce(A-A.diff(X, 2)-beta*A*A) == 0
    assert reduce(-A.diff(X, 3)+(1-2*beta*A)*A.diff(X)) == 0
    assert sp.simplify(A.diff(X).subs(X, -X)+A.diff(X)) == 0


def test_column_dynamics_from_full_cylindrical_euler_and_energy_flux():
    r = sp.Symbol("r", positive=True)
    z, t, amplitude = sp.symbols("z t amplitude", real=True)
    L = sp.Function("L")(r)
    psi, chi, p = (sp.Function(name)(r, z, t) for name in ("psi", "chi", "p"))
    ur, ut, uz = -amplitude*psi.diff(z)/r, (L+amplitude*chi)/r, amplitude*psi.diff(r)/r
    ar = ur.diff(t)+ur*ur.diff(r)+uz*ur.diff(z)-ut**2/r+amplitude*p.diff(r)
    at = ut.diff(t)+ur*ut.diff(r)+uz*ut.diff(z)+ur*ut/r
    az = uz.diff(t)+ur*uz.diff(r)+uz*uz.diff(z)+amplitude*p.diff(z)
    eta = -(psi.diff(r, 2)-psi.diff(r)/r+psi.diff(z, 2))/r**2
    curl_row = sp.diff(ar, z)-sp.diff(az, r)
    linear_curl = sp.diff(curl_row, amplitude).subs(amplitude, 0)/r
    linear_swirl = r*sp.diff(at, amplitude).subs(amplitude, 0)
    data = column_linear_dynamics(L, psi, chi, r, z)
    assert sp.simplify(linear_curl-eta.diff(t)+data.vorticity_rhs) == 0
    assert sp.simplify(linear_swirl-chi.diff(t)+data.angular_momentum_rhs) == 0
    energy_rate = psi*data.vorticity_rhs+data.swirl_energy_weight*chi*data.angular_momentum_rhs
    assert sp.simplify(energy_rate-sp.diff(2*L*psi*chi/r**4, z)) == 0
    # An omitted centrifugal factor cannot pass either curl or energy row.
    assert sp.simplify(psi*data.vorticity_rhs/2+
                       data.swirl_energy_weight*chi*data.angular_momentum_rhs-
                       sp.diff(2*L*psi*chi/r**4, z)) != 0


def test_column_metric_and_pressure_dispersion_are_the_same_physical_system():
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    k, omega = sp.symbols("k omega", nonzero=True, real=True)
    L, f = sp.Function("L")(r), sp.Function("f")(r)
    phase = sp.exp(sp.I*k*z)
    psi = f*phase
    chi = -k/omega*L.diff(r)/r*psi
    data = column_linear_dynamics(L, psi, chi, r, z)
    eta = -(f.diff(r, 2)-f.diff(r)/r-k*k*f)*phase/r**2
    Phi = 2*L*L.diff(r)/r**3
    dispersion = (omega/k)**2*(-f.diff(r, 2)+f.diff(r)/r+k*k*f)-Phi*f
    assert sp.simplify((-sp.I*omega*eta-data.vorticity_rhs)*
                       (omega*r*r/( -sp.I*k*k*phase))-dispersion) == 0
    assert sp.simplify(-sp.I*omega*chi-data.angular_momentum_rhs) == 0
    assert sp.simplify(data.swirl_energy_weight-2*L/(r**3*L.diff(r))) == 0
    solid = column_linear_dynamics(r*r, psi, chi, r, z)
    assert solid.casimir_second_derivative == 0
    assert solid.swirl_energy_weight == 1/r**2


def test_column_metric_rejects_flat_or_negative_circulation_slope():
    import pytest
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    for L in (1, -r, r*r+z):
        with pytest.raises(ValueError):
            column_linear_dynamics(L, r*r*z, r*z, r, z)
