"""Conditional, unpromoted exact Euler displacement identities (P251/0037).

These APIs retain the complete background transport and pressure Hessian.
They supply neither an invariant finite-dimensional reduction nor a Cosserat
constitutive closure. Fields are smooth real symbolic fields in Euclidean
coordinates; density is constant and positive. The Jacobi action additionally
requires an on-shell stationary incompressible background, divergence-free
displacements, and periodic or boundary-compatible variations eliminating the
stated integration-by-parts flux. A material parcel with pressure traction
instead retains that boundary flux. No computation occurs on import.
"""

from __future__ import annotations

import sympy as sp


def material_derivative(field, velocity, coordinates, time):
    """Return componentwise ``partial_t field + velocity.grad(field)``."""
    values, flow = sp.Matrix(field), sp.Matrix(velocity)
    if values.cols != 1 or flow.shape != (len(coordinates), 1):
        raise ValueError("field and velocity must be compatible column vectors")
    return values.diff(time) + values.jacobian(coordinates) * flow


def euler_displacement_perturbation(displacement, velocity, coordinates, time):
    """Eulerian velocity perturbation ``Dt xi - (xi.grad) velocity``.

    This is the tangent map from material displacement to Eulerian velocity,
    not the assertion that an arbitrary displacement solves linearized Euler.
    """
    xi, flow = sp.Matrix(displacement), sp.Matrix(velocity)
    if xi.shape != flow.shape:
        raise ValueError("displacement and velocity must have the same shape")
    return material_derivative(xi, flow, coordinates, time) - flow.jacobian(
        coordinates
    ) * xi


def euler_jacobi_density(displacement, velocity, pressure, density, coordinates, time):
    """Return ``(rho*|Dt xi|² - xi.Hess(p).xi)/2``.

    This is the quadratic action on div(xi)=0 about an exact stationary Euler
    solution, after pressure-constraint integration by parts. It includes the
    positive time-kinetic term, gyroscopic mixed term and static transport term.
    Positivity of its time-kinetic metric does not imply spectral stability.
    Units: xi length, velocity length/time, p mass/(length*time²), rho
    mass/length³; the density has energy/volume units.
    """
    rho = sp.sympify(density)
    if rho.is_positive is False or (rho.is_number and rho.is_positive is not True):
        raise ValueError("density must be positive")
    xi = sp.Matrix(displacement)
    advected = material_derivative(xi, velocity, coordinates, time)
    return (rho * advected.dot(advected) - (xi.T * sp.hessian(
        pressure, coordinates
    ) * xi)[0]) / 2


def stationarize_planar(streamfunction, profile_pressure, rotation, density, coordinates):
    """Convert a planar rotating Euler profile to a changed stationary field.

    With J(x,y)=(-y,x), v=J grad(psi), and a rigidly rotating Euler solution
    R(Omega*t)v(R(-Omega*t)x), return w=v-Omega*Jx and its pressure.
    Assumes a simply connected planar domain (or consistent streamfunction).
    The physical vorticity changes by -2*Omega. The far field and boundary
    conditions change accordingly; the result is not a coordinate-only
    equivalence, a finite-energy guarantee, or a general 3-D transformation.
    The routine constructs expressions; the caller checks the rotating-profile
    Euler equation and compatibility of its domain/boundary conditions.
    """
    if len(coordinates) != 2:
        raise ValueError("planar stationarization requires two coordinates")
    x, y = coordinates
    psi, omega, rho = map(sp.sympify, (streamfunction, rotation, density))
    if rho.is_positive is False or (rho.is_number and rho.is_positive is not True):
        raise ValueError("density must be positive")
    flow = sp.Matrix([-sp.diff(psi, y) + omega * y, sp.diff(psi, x) - omega * x])
    pressure = profile_pressure - 2 * rho * omega * psi + rho * omega**2 * (x*x + y*y) / 2
    return flow, pressure
