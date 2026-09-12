"""Exact metric, Hodge, topology and steady-transport exposing calculations."""

import sympy as sp

from substrate_framework.euler_twisted_carrier import (
    radial_l1_poisson_inverse, radial_twisted_carrier,
)


def test_berry_form_from_actual_unit_spinor_and_spherical_metric():
    r = sp.Symbol("r", positive=True)
    th, ph = sp.symbols("theta phi", real=True)
    k, rho = sp.symbols("kappa rho", positive=True)
    F = sp.Function("F", real=True)(r)
    data = radial_twisted_carrier(F, r, th, density=rho, circulation_scale=k)
    z = sp.Matrix([sp.cos(F)+sp.I*sp.cos(th)*sp.sin(F),
                   (-sp.sin(ph)+sp.I*sp.cos(ph))*sp.sin(th)*sp.sin(F)])
    direct = sp.Matrix([-sp.I*k*(z.conjugate().T*z.diff(q))[0] for q in (r, th, ph)])
    assert (direct-data.berry_oneform).applyfunc(sp.simplify) == sp.zeros(3, 1)
    ar, at, ap = data.berry_oneform
    divergence = sp.diff(r*r*ar, r)/(r*r)+sp.diff(sp.sin(th)*at, th)/(r*r*sp.sin(th))
    assert sp.trigsimp(divergence-k*sp.cos(th)*data.radial_divergence_source) == 0
    # Missing spherical metric factors would change both divergence and j.
    angular_moment = sp.integrate(r*sp.sin(th)*(ap/(r*sp.sin(th)))*r*r*sp.sin(th),
                                 (th, 0, sp.pi))*2*sp.pi*rho
    assert sp.simplify(angular_moment-data.axial_moment_radial_integrand) == 0


def test_full_l1_poisson_inverse_and_physical_divergence():
    r = sp.Symbol("r", positive=True)
    source = sp.Function("S")(r)
    h = radial_l1_poisson_inverse(source, r)
    operator = sp.diff(h, r, 2)+2*sp.diff(h, r)/r-2*h/r**2
    assert sp.simplify(operator-source) == 0
    th = sp.Symbol("theta", real=True)
    data = radial_twisted_carrier(sp.Function("F", real=True)(r), r, th,
                                  density=1, circulation_scale=1)
    ur, ut, _ = data.physical_spherical_velocity
    divergence = sp.diff(r*r*ur, r)/(r*r)+sp.diff(sp.sin(th)*ut, th)/(r*sp.sin(th))
    assert sp.trigsimp(sp.simplify(divergence)) == 0


def test_helicity_from_actual_threeform_and_profile_orientation():
    r = sp.Symbol("r", positive=True)
    th = sp.Symbol("theta", real=True)
    k = sp.Symbol("kappa", positive=True)
    # r here is the profile-value integration variable after dF=F' dr.
    data = radial_twisted_carrier(r, r, th, density=1, circulation_scale=k)
    helicity = sp.integrate(data.helicity_threeform_coefficient,
                            (th, 0, sp.pi), (r, sp.pi, 0))*2*sp.pi
    assert helicity == -4*sp.pi**2*k*k


def test_common_translation_cannot_satisfy_flat_swirl_transport_edge():
    r, b, C, k = sp.symbols("r b C kappa", positive=True)
    U = sp.Symbol("U", real=True)
    # Exact exterior Green potential and its actual poloidal streamfunction.
    h = -C/(3*r*r)
    G = -r*r*sp.diff(h, r)/2
    Q = k*G-U*r*r/2
    matched_speed = sp.solve(Q.subs(r, b), U)[0]
    derivative = sp.simplify(sp.diff(Q, r).subs({r: b, U: matched_speed}))
    assert derivative == k*C/b**2
    assert derivative.is_positive
