"""Exact leading multipoles for localized incompressible Euler fields.

The helpers expose asymptotic coefficients.  They do not identify any Euler
invariant with electric charge or construct an electromagnetic sector.
"""

import sympy as sp


def _vector3(value, name):
    vector = sp.ImmutableMatrix(value)
    if vector.shape != (3, 1):
        raise ValueError(f"{name} must have three components")
    return vector


def vortex_impulse_far_velocity(impulse, displacement):
    """Return the leading compact-vorticity velocity dipole at displacement.

    The convention is ``I = (1/2) integral x cross omega dx`` and the returned
    coefficient includes the ``1/(4*pi)`` Biot--Savart normalization.
    """

    impulse_vector = _vector3(impulse, "impulse")
    position = _vector3(displacement, "displacement")
    radius_squared = sp.simplify(position.dot(position))
    if radius_squared == 0:
        raise ValueError("displacement must be nonzero")
    radius = sp.sqrt(radius_squared)
    return sp.ImmutableMatrix(
        (3 * position * impulse_vector.dot(position) / radius_squared - impulse_vector)
        / (4 * sp.pi * radius**3)
    )


def vortex_dipole_cross_energy(density, impulse_1, impulse_2, displacement):
    """Return the leading two-carrier kinetic cross energy.

    The result is ``rho I1_i I2_j partial_i partial_j(1/(4*pi*r))``.
    """

    if density <= 0:
        raise ValueError("density must be positive")
    first = _vector3(impulse_1, "impulse_1")
    second = _vector3(impulse_2, "impulse_2")
    position = _vector3(displacement, "displacement")
    radius_squared = sp.simplify(position.dot(position))
    if radius_squared == 0:
        raise ValueError("displacement must be nonzero")
    radius = sp.sqrt(radius_squared)
    numerator = 3 * first.dot(position) * second.dot(position) / radius_squared - first.dot(second)
    return sp.simplify(density * numerator / (4 * sp.pi * radius**3))


def pressure_quadrupole_far_field(density, kinetic_moment, displacement):
    """Return the leading physical-pressure field of localized Euler data.

    ``kinetic_moment[i,j]`` is ``integral u_i u_j dx``.  The pressure
    convention is ``p=rho*(-Delta)^-1 partial_i partial_j(u_i*u_j)``.
    """

    if density <= 0:
        raise ValueError("density must be positive")
    moment = sp.ImmutableMatrix(kinetic_moment)
    if moment.shape != (3, 3):
        raise ValueError("kinetic_moment must be 3 by 3")
    if moment != moment.T:
        raise ValueError("kinetic_moment must be symmetric")
    position = _vector3(displacement, "displacement")
    radius_squared = sp.simplify(position.dot(position))
    if radius_squared == 0:
        raise ValueError("displacement must be nonzero")
    radius = sp.sqrt(radius_squared)
    contraction = 3 * (position.T * moment * position)[0] / radius_squared - sp.trace(moment)
    return sp.simplify(density * contraction / (4 * sp.pi * radius**3))


def radial_monopole_flux(strength):
    """Return flux of ``strength*x/|x|^3`` through any centered sphere."""

    return sp.simplify(4 * sp.pi * strength)


def internal_stress_source(wavevector, stress):
    """Return the Fourier source ``i*k_j*stress_ij`` of an internal stress."""

    covector = _vector3(wavevector, "wavevector")
    tensor = sp.ImmutableMatrix(stress)
    if tensor.shape != (3, 3):
        raise ValueError("stress must be 3 by 3")
    return sp.ImmutableMatrix(sp.I * tensor * covector)


def transverse_static_green(shear_modulus, displacement):
    """Return the three-dimensional transverse ``-mu*Delta`` Green tensor."""

    if shear_modulus <= 0:
        raise ValueError("shear_modulus must be positive")
    position = _vector3(displacement, "displacement")
    radius_squared = sp.simplify(position.dot(position))
    if radius_squared == 0:
        raise ValueError("displacement must be nonzero")
    radius = sp.sqrt(radius_squared)
    direction = position / radius
    return sp.ImmutableMatrix(
        (sp.eye(3) + direction * direction.T) / (8 * sp.pi * shear_modulus * radius)
    )


def transverse_scalar_source(wavevector, form_factor=1):
    """Project the isotropic scalar source ``i*k*F(k^2)`` transversely."""

    covector = _vector3(wavevector, "wavevector")
    norm_squared = sp.simplify(covector.dot(covector))
    if norm_squared == 0:
        raise ValueError("wavevector must be nonzero")
    projector = sp.eye(3) - covector * covector.T / norm_squared
    return sp.ImmutableMatrix(sp.simplify(projector * (sp.I * form_factor * covector)))
