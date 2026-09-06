"""Conditional spherical-profile/Hodge data for a twisted Euler initial state.

P253/0017 derives smooth localization when F is pi near zero and zero outside
a finite radius, with smooth flat joining. This API computes exact expressions
from supplied profiles; it does not certify these hypotheses, steady Euler,
nonlinear stability, quantum charge, or particle identity.
"""

from dataclasses import dataclass

import sympy as sp


def radial_l1_poisson_inverse(source, radius):
    """Return the regular-at-zero, decaying l=1 radial Poisson Green integral.

    The source must be regular and have convergent displayed integrals; smooth
    compact support away from zero suffices. No convergence is inferred for
    arbitrary symbolic input. The radial operator is h''+2h'/r-2h/r**2.
    """
    if not isinstance(radius, sp.Symbol) or radius.is_positive is False:
        raise ValueError("radius must be a symbol on the positive radial domain")
    r = radius
    s = sp.Dummy("s", positive=True)
    forcing = sp.sympify(source).subs(r, s)
    return -(sp.Integral(s**3*forcing, (s, 0, r))/r**2
             + r*sp.Integral(forcing, (s, r, sp.oo)))/3


@dataclass(frozen=True)
class RadialTwistedCarrier:
    """Coordinate one-form, metric divergence, Hodge field and actual moments."""

    berry_oneform: sp.ImmutableMatrix
    radial_divergence_source: sp.Expr
    hodge_radial_potential: sp.Expr
    physical_spherical_velocity: sp.ImmutableMatrix
    helicity_threeform_coefficient: sp.Expr
    axial_moment_radial_integrand: sp.Expr


def radial_twisted_carrier(profile, radius, polar_angle, *, density, circulation_scale):
    """Return data for z0=(cos F+i*nz*sin F,(-ny+i*nx)*sin F).

    The Berry entries multiply (dr,dtheta,dphi); they are coordinate one-form
    coefficients, not orthonormal vector components. The physical velocity
    entries instead are orthonormal spherical components after Hodge projection.
    The returned three-form coefficient multiplies dr wedge dtheta wedge dphi.
    The axial moment integrand already includes angular integration and rho.
    Positive finite constant rho/kappa, real radial F, smooth flat edges and
    Green-integral convergence remain explicit caller hypotheses when undecidable.
    """
    r, th = radius, polar_angle
    if not isinstance(r, sp.Symbol) or not isinstance(th, sp.Symbol) or r == th:
        raise ValueError("radius and polar angle must be distinct symbols")
    F = sp.sympify(profile)
    rho, k = sp.sympify(density), sp.sympify(circulation_scale)
    if F.has(th) or F.is_real is False:
        raise ValueError("profile must be real and independent of polar angle")
    for value in (rho, k):
        if (value.is_positive is False or value.is_finite is False
                or value.is_real is False or value.has(r, th)):
            raise ValueError("density and circulation scale must be positive finite constants")
    ar = k*sp.cos(th)*sp.diff(F, r)
    at = -k*sp.sin(F)*sp.cos(F)*sp.sin(th)
    ap = k*sp.sin(F)**2*sp.sin(th)**2
    source = sp.diff(F, r, 2)+2*sp.diff(F, r)/r-sp.sin(2*F)/r**2
    h = radial_l1_poisson_inverse(source, r)
    velocity = sp.ImmutableMatrix([
        k*(sp.diff(F, r)-sp.diff(h, r))*sp.cos(th),
        k*(h-sp.sin(F)*sp.cos(F))*sp.sin(th)/r,
        k*sp.sin(F)**2*sp.sin(th)/r,
    ])
    helicity = sp.trigsimp(ar*sp.diff(ap, th)-at*sp.diff(ap, r)
                          +ap*(sp.diff(at, r)-sp.diff(ar, th)))
    return RadialTwistedCarrier(sp.ImmutableMatrix([ar, at, ap]), source, h,
                                velocity, helicity,
                                8*sp.pi*rho*k*r**2*sp.sin(F)**2/3)
