"""Exact Ertel-current identities for transported Euler material labels."""

from __future__ import annotations

import sympy as sp


def ertel_stretching_residual(vorticity, tag_gradient, velocity_gradient):
    """Return ``D_t(omega dot grad chi)`` after Euler/tag substitution.

    ``velocity_gradient[i, j]`` means ``partial_j u_i``.  Incompressible
    Euler gives ``D_t omega = velocity_gradient * omega`` and a transported
    scalar gives ``D_t grad(chi) = -velocity_gradient.T * grad(chi)``.
    """

    omega = sp.Matrix(vorticity)
    grad_chi = sp.Matrix(tag_gradient)
    grad_u = sp.Matrix(velocity_gradient)
    return sp.expand((grad_u * omega).dot(grad_chi) - omega.dot(grad_u.T * grad_chi))


def ertel_current_parity(tag_parity: str) -> dict[str, str]:
    """Return O(3) parity for ``q_E=omega dot grad chi`` and ``q_E*u``."""

    if tag_parity == "scalar":
        return {"density": "pseudoscalar", "spatial_current": "axial"}
    if tag_parity == "pseudoscalar":
        return {"density": "scalar", "spatial_current": "polar"}
    raise ValueError("tag_parity must be 'scalar' or 'pseudoscalar'")


def ertel_charge_from_flux(boundary_flux):
    """Represent ``integral q_E`` by the boundary flux of ``chi*omega``.

    The caller supplies the oriented boundary flux.  It is zero for a smooth
    globally single-valued label with compactly supported ``chi*omega`` on
    all of R^3, and on any boundaryless domain where Stokes applies.
    """

    return sp.sympify(boundary_flux)


def forced_ertel_source(force_curl, tag_gradient):
    """Return ``(curl f) dot grad(chi)`` in the forced Ertel equation."""

    return sp.Matrix(force_curl).dot(sp.Matrix(tag_gradient))


def forced_ertel_flux(density, velocity, force, tag_gradient):
    """Return the conserved flux ``q_E*u-f cross grad(chi)``."""

    q = sp.sympify(density)
    return q * sp.Matrix(velocity) - sp.Matrix(force).cross(sp.Matrix(tag_gradient))


def cao_forced_magnetization_flux(
    coupling, mass_density, tag, tag_derivative, grad_phi, grad_streamfunction
):
    """Forced Cao flux ``(g/rho) chi chi' grad(Phi) cross grad(P)``."""

    coefficient = (
        sp.sympify(coupling)
        * sp.sympify(tag)
        * sp.sympify(tag_derivative)
        / sp.sympify(mass_density)
    )
    return coefficient * sp.Matrix(grad_phi).cross(sp.Matrix(grad_streamfunction))


def transported_lock_residual(tag, material_lambda_derivative):
    """Return ``D_t(q_E-lambda*chi)`` for unforced material invariants."""

    return -sp.sympify(tag) * sp.sympify(material_lambda_derivative)


def forced_lock_residual(
    tag, material_lambda_derivative, force_curl, tag_gradient
):
    """Return the forced lock equation residual.

    Vanishing is ``chi D_t(lambda)=(curl f) dot grad(chi)``.
    """

    return (
        sp.sympify(tag) * sp.sympify(material_lambda_derivative)
        - sp.Matrix(force_curl).dot(sp.Matrix(tag_gradient))
    )


def closed_vorticity_line_multiplier(integrated_lambda):
    """Monodromy multiplier for ``d chi/d tau=lambda chi`` on a closed line."""

    return sp.exp(sp.sympify(integrated_lambda))


def cao_ertel_density(zeta, theta_derivative_tag):
    """Cao convention ``omega=r*zeta*e_theta=zeta*partial_theta``.

    Hence ``omega dot grad(chi)=zeta*partial_theta(chi)``.  Every smooth
    axisymmetric tag has zero Ertel density, including ``chi=F(I)``.
    """

    return sp.sympify(zeta) * sp.sympify(theta_derivative_tag)


def cao_azimuthal_phase_density(zeta):
    """Local Ertel density for the circle phase on R^3 minus the axis."""

    return sp.sympify(zeta)


def cao_azimuthal_phase_charge(circulation):
    """Return the equilibrium integral ``int zeta dV=2*pi*kappa``."""

    return 2 * sp.pi * sp.sympify(circulation)


def geometric_azimuth_material_derivative(azimuthal_velocity, radius):
    """Return ``D_t theta=u_theta/r`` for the fixed geometric azimuth."""

    return sp.sympify(azimuthal_velocity) / sp.sympify(radius)
