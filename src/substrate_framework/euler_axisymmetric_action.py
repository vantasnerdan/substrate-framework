"""Exact axisymmetric Euler Clebsch equations and momentum surface terms.

The meridional measure is r dr dz and psi=K zeta is the whole-space Hodge
relation. Expressions here are local identities. Global chart coverage,
axial end jumps, finite relative action, and stability remain hypotheses.
"""

from dataclasses import dataclass

import sympy as sp


def axisymmetric_bracket(first, second, radius, axial_coordinate):
    """Meridional Poisson bracket for area form r dr wedge dz."""
    r, z = radius, axial_coordinate
    if not isinstance(r, sp.Symbol) or r.is_nonpositive is True:
        raise ValueError("radius must be a symbol on r>0")
    if not isinstance(z, sp.Symbol) or z == r:
        raise ValueError("axial coordinate must be a distinct symbol")
    f, g = map(sp.sympify, (first, second))
    return (sp.diff(f, r)*sp.diff(g, z)-sp.diff(f, z)*sp.diff(g, r))/r


@dataclass(frozen=True)
class AxisymmetricCanonicalEuler:
    poloidal_vorticity: sp.Expr
    swirl_rhs: sp.Expr
    clebsch_rhs: sp.Expr
    euler_vorticity_rhs: sp.Expr
    normalized_hamiltonian_density: sp.Expr
    normalized_translation_density: sp.Expr
    impulse_boundary_oneform: sp.Matrix


def axisymmetric_canonical_euler(swirl, clebsch, streamfunction,
                                radius, axial_coordinate):
    """Canonical equations for xi=r*u_theta, beta and zeta={xi,beta}.

    The action is 2*pi*rho * integral(xi*beta_t-H) r dr dz dt.
    Its symplectic form is 2*pi*rho * integral d beta wedge d xi r dr dz.
    The boundary one-form (dr,dz components) satisfies
    I_D=integral_D beta*xi_z r dr dz + integral_boundary_D oneform
    for I_D=integral_D r**2*zeta/2 r dr dz and positive orientation.
    Thus a nonzero end term is retained, not silently discarded.
    """
    r, z = radius, axial_coordinate
    xi, beta, psi = map(sp.sympify, (swirl, clebsch, streamfunction))
    def bracket(f, g):
        return axisymmetric_bracket(f, g, r, z)
    zeta = bracket(xi, beta)
    return AxisymmetricCanonicalEuler(
        zeta,
        -bracket(psi, xi),
        xi/r**2+bracket(beta, psi),
        -bracket(psi, zeta)+2*xi*sp.diff(xi, z)/r**4,
        (zeta*psi+xi**2/r**2)/2,
        beta*sp.diff(xi, z),
        -r**2*beta*sp.Matrix([sp.diff(xi, r), sp.diff(xi, z)])/2,
    )


def axisymmetric_linearized_euler(base_streamfunction, base_vorticity,
                                  base_swirl, streamfunction_variation,
                                  vorticity_variation, swirl_variation,
                                  radius, axial_coordinate, speed=0):
    """Full (eta_t, chi_t) Euler operator in a frame traveling at speed c.

    The inputs obey psi=K*zeta and delta_psi=K*eta on the physical domain.
    No core/exterior truncation or stationarity is certified by this API.
    Retaining delta_psi is essential: it contains the pressure/Hodge coupling
    and the base vorticity/label-gradient terms.
    """
    r, z = radius, axial_coordinate
    psi, zeta, xi, dpsi, eta, chi, c = map(sp.sympify, (
        base_streamfunction, base_vorticity, base_swirl,
        streamfunction_variation, vorticity_variation, swirl_variation, speed))
    if c.has(r, z):
        raise ValueError("frame speed must be spatially constant")
    def bracket(f, g):
        return axisymmetric_bracket(f, g, r, z)
    frame_psi = psi-c*r**2/2
    return sp.Matrix([
        -bracket(frame_psi, eta)-bracket(dpsi, zeta)+2*sp.diff(xi*chi, z)/r**4,
        -bracket(frame_psi, chi)-bracket(dpsi, xi),
    ])
