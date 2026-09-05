"""Conditional linear Euler modes of an isolated Rankine vortex.

Unpromoted P251 infrastructure, issue #198. The core rotates with Omega>0;
the exterior velocity is Omega*a^2/r, density is constant and perturbations
are exp(i*(m*theta+k*z-omega*t)). Exterior perturbations are irrotational,
u'=grad(potential), with the decaying modified-Bessel radial profile.
Pressure and radial displacement match at
r=a. These APIs establish mode equations, not Cosserat elastic moduli.
The Bessel residual is dimensionless: x=|k|*a, s=omega/Omega-m. Its real
J_m branch assumes 0<|s|<2 and excludes pressure-amplitude poles.
"""
from __future__ import annotations

from typing import Any

import sympy as sp


def core_velocity(pressure: Any, radius: Any, axial_k: Any, azimuthal_m: Any,
                  doppler: Any, rotation: Any, density: Any) -> sp.ImmutableMatrix:
    """Solve the Cartesian-derived radial/azimuthal/axial Euler equations.

Nonzero density, radius, doppler and 4*rotation**2-doppler**2 are required.
Pressure is a symbolic radial function; no time integration is performed.
"""
    p, r, k, m, s, Om, rho = map(sp.sympify, (
        pressure, radius, axial_k, azimuthal_m, doppler, rotation, density))
    matrix = sp.Matrix([[-sp.I*s, -2*Om, 0], [2*Om, -sp.I*s, 0], [0, 0, -sp.I*s]])
    gradient = sp.Matrix([sp.diff(p, r), sp.I*m*p/r, sp.I*k*p]) / rho
    return sp.ImmutableMatrix(matrix.inv() * -gradient)


def boundary_determinant(doppler: Any, rotation: Any, azimuthal_m: Any,
                         interior_log_derivative: Any, exterior_log_derivative: Any) -> sp.Expr:
    """Pressure match using core a*P'/P and exterior a*potential'/potential.

    Exterior pressure includes a radially varying Doppler factor. Its own
    logarithmic derivative is not the supplied exterior argument.
    """
    s, Om, m, J, K = map(sp.sympify, (doppler, rotation, azimuthal_m,
                                      interior_log_derivative, exterior_log_derivative))
    return sp.expand((4*Om**2-s**2)*K + s**2*J - 2*Om*m*s)


def rankine_residual(context: Any, x: Any, azimuthal_m: int, doppler: Any) -> Any:
    """Exact Bessel residual in a caller-owned mpmath precision context.

Caller owns root selection, residual acceptance and precision refinement.
No global mpmath precision or library settings are modified.
"""
    c = context
    x, s = c.mpf(x), c.mpf(doppler)
    if azimuthal_m < 1 or not isinstance(azimuthal_m, int):
        raise ValueError("azimuthal_m must be a positive integer")
    if not (x > 0 and 0 < abs(s) < 2):
        raise ValueError("require x>0 and 0<|doppler|<2 on the real J branch")
    m = azimuthal_m
    lam = x*c.sqrt(4-s*s)/abs(s)
    jm = c.besselj(m, lam)
    km = c.besselk(m, x)
    jr = lam*c.besselj(m-1, lam)/jm-m
    kr = -x*c.besselk(m-1, x)/km-m
    return (4-s*s)*kr+s*s*jr-2*m*s
