"""Expose actual Euler source, gauge, measure and end-term identities."""

import sympy as sp

from substrate_framework.euler_axisymmetric_action import (
    axisymmetric_bracket, axisymmetric_canonical_euler,
    axisymmetric_linearized_euler,
)


def test_canonical_evolution_is_full_axisymmetric_euler_vorticity():
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    xi, beta, psi = (sp.Function(n)(r, z) for n in ("xi", "beta", "psi"))
    data = axisymmetric_canonical_euler(xi, beta, psi, r, z)
    def bracket(f, g):
        return axisymmetric_bracket(f, g, r, z)
    differentiated = bracket(data.swirl_rhs, beta)+bracket(xi, data.clebsch_rhs)
    assert sp.simplify(differentiated-data.euler_vorticity_rhs) == 0
    passive = -bracket(psi, data.poloidal_vorticity)
    assert sp.simplify(differentiated-passive) != 0
    # Reversing the centrifugal term produces the opposite swirl source.
    wrong = bracket(data.swirl_rhs, beta)+bracket(xi, data.clebsch_rhs-2*xi/r**2)
    assert sp.simplify(wrong-data.euler_vorticity_rhs) != 0


def test_clebsch_gauge_preserves_field_and_canonical_dynamics():
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    xi, beta, psi = (sp.Function(n)(r, z) for n in ("xi", "beta", "psi"))
    F = sp.Function("F")(xi)
    data = axisymmetric_canonical_euler(xi, beta, psi, r, z)
    gauged = axisymmetric_canonical_euler(xi, beta+F, psi, r, z)
    assert sp.simplify(data.poloidal_vorticity-gauged.poloidal_vorticity) == 0
    assert sp.simplify(gauged.clebsch_rhs-data.clebsch_rhs-
                       sp.diff(F, xi)*data.swirl_rhs) == 0


def test_impulse_surface_term_from_green_identity_and_nonzero_rectangle():
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    xi, beta = sp.Function("xi")(r, z), sp.Function("beta")(r, z)
    data = axisymmetric_canonical_euler(xi, beta, 0, r, z)
    br, bz = data.impulse_boundary_oneform
    boundary_density = sp.diff(bz, r)-sp.diff(br, z)
    bulk = r**3*data.poloidal_vorticity/2-r*data.normalized_translation_density
    assert sp.simplify(bulk-boundary_density) == 0
    sample = axisymmetric_canonical_euler(r*r+z, r*z, 0, r, z)
    br, bz = sample.impulse_boundary_oneform
    boundary = (sp.integrate(br.subs(z, -1), (r, 1, 2))+
                sp.integrate(bz.subs(r, 2), (z, -1, 1))+
                sp.integrate(br.subs(z, 1), (r, 2, 1))+
                sp.integrate(bz.subs(r, 1), (z, 1, -1)))
    impulse = sp.integrate(r**3*sample.poloidal_vorticity/2, (z, -1, 1), (r, 1, 2))
    momentum = sp.integrate(r*sample.normalized_translation_density, (z, -1, 1), (r, 1, 2))
    assert impulse == momentum+boundary
    assert boundary != 0 and impulse != momentum-boundary


def test_mixed_casimir_is_an_axial_end_jump_with_physical_measure():
    r = sp.Symbol("r", positive=True)
    z = sp.Symbol("z", real=True)
    xi, beta = r*r, r**4*z
    data = axisymmetric_canonical_euler(xi, beta, 0, r, z)
    a = sp.Symbol("a", positive=True)
    D = a**2
    direct = sp.integrate(data.poloidal_vorticity*D.subs(a, xi)*r,
                          (z, -1, 1), (r, 1, 2))
    jump = beta.subs(r, sp.sqrt(a)).subs(z, 1)-beta.subs(r, sp.sqrt(a)).subs(z, -1)
    assert direct == sp.integrate(D*jump, (a, 1, 4))
    rho = sp.Symbol("rho", positive=True)
    physical = 2*sp.pi*rho*direct
    assert physical != direct and jump != 0


def test_full_solitary_linearization_is_actual_euler_variation():
    r = sp.Symbol("r", positive=True)
    z, amplitude, c = sp.symbols("z amplitude c", real=True)
    psi, zeta, xi, dpsi, eta, chi = (
        sp.Function(n)(r, z) for n in ("psi", "zeta", "xi", "dpsi", "eta", "chi"))
    def bracket(f, g):
        return axisymmetric_bracket(f, g, r, z)
    ps, ze, xx = psi+amplitude*dpsi, zeta+amplitude*eta, xi+amplitude*chi
    nonlinear = sp.Matrix([
        c*sp.diff(ze, z)-bracket(ps, ze)+sp.diff(xx**2, z)/r**4,
        c*sp.diff(xx, z)-bracket(ps, xx),
    ])
    actual = axisymmetric_linearized_euler(psi, zeta, xi, dpsi, eta, chi, r, z, c)
    assert (actual-nonlinear.diff(amplitude).subs(amplitude, 0)).applyfunc(sp.simplify) == sp.zeros(2, 1)
    # Translation derivative commutes with the exact steady residual.
    neutral = axisymmetric_linearized_euler(psi, zeta, xi,
        psi.diff(z), zeta.diff(z), xi.diff(z), r, z, c)
    assert (neutral-nonlinear.subs(amplitude, 0).diff(z)).applyfunc(sp.simplify) == sp.zeros(2, 1)


def test_exact_solitary_label_jump_and_first_nonzero_coefficient():
    a, c = sp.symbols("a c", positive=True)
    z, amplitude = sp.symbols("z amplitude", real=True)
    L = sp.Function("L")(a)
    R = sp.Function("R")(a, z)
    # F(psi0(a))=L(a), psi0=-c*a^2/2; B'=F*F'/a^2.
    Fprime = L.diff(a)/(-c*a)
    Bprime = L*Fprime/a**2
    zeta = L*Fprime/R**2-Bprime
    xi_r = L.diff(a)/R.diff(a)
    density = zeta*R/xi_r
    exact = L/(c*a)*(R*R.diff(a)/a**2-R.diff(a)/R)
    assert sp.simplify(density-exact) == 0
    h = sp.Function("h")(a, z)
    perturbed = exact.subs(R, a+amplitude*h).doit()
    assert sp.simplify(perturbed.subs(amplitude, 0)) == 0
    assert sp.simplify(perturbed.diff(amplitude).subs(amplitude, 0)-2*L*h/(c*a**3)) == 0
    # The derivative of the radial displacement cancels at first order.
    assert not sp.simplify(perturbed.diff(amplitude).subs(amplitude, 0)).has(h.diff(a))
