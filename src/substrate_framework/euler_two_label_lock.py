"""Exact algebra for the Cao two-label forced axial-lock supplier.

The electromagnetic material tag and the circle-valued axial phase are
distinct labels.  These helpers encode only pointwise and contour-projection
identities; they do not construct a moving defect, a Cao branch, or a flavor
interaction.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def exponential_column_profile(
    radius, center_value, profile_exponent, radial_scale, major_radius
):
    """Return an explicit radial Liouville column ``(P, zeta)``.

    The pair satisfies ``-Delta_2 P=R**2*zeta`` and, away from the center,
    ``d(zeta)/d(P)=profile_exponent*zeta``.  Its vorticity decays like
    ``radius**(-4)`` while its velocity has a ``radius**(-1)`` tail; this is
    a response-kernel profile, not a finite-energy infinite column.
    """

    s = sp.sympify(radius)
    p0 = sp.sympify(center_value)
    exponent = sp.sympify(profile_exponent)
    scale = sp.sympify(radial_scale)
    radius_major = sp.sympify(major_radius)
    if exponent.is_zero is True or radius_major.is_zero is True:
        raise ValueError("profile_exponent and major_radius must be nonzero")
    if scale.is_nonpositive is True:
        raise ValueError("radial_scale must be positive")
    denominator = 1 + scale * s**2
    streamfunction = p0 - 2 * sp.log(denominator) / exponent
    vorticity = 8 * scale / (
        exponent * radius_major**2 * denominator**2
    )
    return sp.simplify(streamfunction), sp.simplify(vorticity)


def two_label_forced_lock_residual(
    tag, material_lambda_derivative, force_curl, phase_gradient
):
    """Return ``chi*D_t(lambda)-(curl f).grad(Theta)``.

    Vanishing is the forced compatibility equation for
    ``omega.grad(Theta)=lambda*chi`` when both ``Theta`` and ``chi`` are
    transported.  The formula is meaningful on the zero set without dividing
    by ``chi``.
    """

    return (
        sp.sympify(tag) * sp.sympify(material_lambda_derivative)
        - sp.Matrix(force_curl).dot(sp.Matrix(phase_gradient))
    )


def constant_ratio_residual(axial_density, tag, constant_ratio):
    """Return the undivided pointwise residual ``q_A-lambda_0*chi``."""

    return (
        sp.sympify(axial_density)
        - sp.sympify(constant_ratio) * sp.sympify(tag)
    )


def cao_constant_ratio_residual(
    profile_value,
    coupling,
    mass_density,
    tag,
    tag_derivative,
    phi,
    toroidal_magnetic_ratio,
    constant_ratio,
):
    """Return the Cao pointwise constant-ratio equation residual.

    The modified Grad--Shafranov first integral is
    ``zeta=h(P)+(g/rho)*(chi'(P)*Phi-chi(P)*H)``.
    """

    correction = (
        sp.sympify(tag_derivative) * sp.sympify(phi)
        - sp.sympify(tag) * sp.sympify(toroidal_magnetic_ratio)
    )
    return (
        sp.sympify(profile_value)
        + sp.sympify(coupling) * correction / sp.sympify(mass_density)
        - sp.sympify(constant_ratio) * sp.sympify(tag)
    )


def weighted_contour_zero_mean(
    values: Sequence, weights: Sequence
) -> tuple[sp.Expr, ...]:
    """Project finitely represented contour data off its weighted mean.

    This is an exact algebraic interface for quadrature-free symbolic rows or
    already-integrated contour sectors.  A continuum proof must separately
    supply the action-period measure and mapping properties.
    """

    if len(values) == 0 or len(values) != len(weights):
        raise ValueError("values and weights must have the same nonzero length")
    vals = tuple(map(sp.sympify, values))
    wts = tuple(map(sp.sympify, weights))
    total = sp.Add(*wts)
    if total.is_zero is True:
        raise ValueError("contour weights must have nonzero total")
    average = sp.Add(*(w * v for w, v in zip(wts, vals, strict=True))) / total
    return tuple(sp.simplify(v - average) for v in vals)


def fixed_profile_maxwell_obstruction(
    tag, tag_derivative, phi_hat_values, magnetic_hat_values, weights
) -> tuple[sp.Expr, ...]:
    """Return the order-``g**2`` contour obstruction for a fixed tag.

    The raw response is ``chi_0'*Phi_hat-chi_0*H_hat``.  Profile, tag and
    constant-ratio corrections that depend only on the streamline label are
    annihilated by the zero-mean projection.
    """

    if len(phi_hat_values) != len(magnetic_hat_values):
        raise ValueError("Maxwell contour rows must have the same length")
    chi = sp.sympify(tag)
    chip = sp.sympify(tag_derivative)
    raw = tuple(
        chip * sp.sympify(phi) - chi * sp.sympify(magnetic)
        for phi, magnetic in zip(
            phi_hat_values, magnetic_hat_values, strict=True
        )
    )
    return weighted_contour_zero_mean(raw, weights)


def comoving_lorenz_scalar_source(
    tag, permittivity, speed, wave_speed, relative_axial_velocity
):
    """Return the source of the signed-charge-factored comoving potential.

    If ``phi=g*phi_hat``, ``A=g*A_hat`` and
    ``Phi=phi-speed*A_z=g*Phi_hat``, the traveling Lorenz equations give

    ``L_c Phi_hat = chi/epsilon *
    (1-speed**2/c_EM**2-speed*W_z/c_EM**2)``.

    In particular the prefactor contains no residual power of the signed
    coupling ``g``.  This helper fixes the relative normalization needed when
    the electric and toroidal-magnetic responses are combined.
    """

    chi = sp.sympify(tag)
    eps = sp.sympify(permittivity)
    c = sp.sympify(speed)
    c_em = sp.sympify(wave_speed)
    w_z = sp.sympify(relative_axial_velocity)
    if eps.is_zero is True or c_em.is_zero is True:
        raise ValueError("permittivity and wave_speed must be nonzero")
    return sp.simplify(
        chi / eps * (1 - c**2 / c_em**2 - c * w_z / c_em**2)
    )


def ampere_first_speed_magnetic_response(
    electric_radial_derivative, wave_speed, major_radius
):
    """Return the first speed derivative of ``H_hat`` at zero speed.

    This follows from the signed-charge-factored Ampere primitive

    ``(1/mu-epsilon*c**2)*r*H_hat + epsilon*c*d_r Phi_hat = K(P)/r``

    and ``epsilon*mu=wave_speed**(-2)`` while the fixed profile and stream
    function are held constant in this response calculation.
    """

    derivative = sp.sympify(electric_radial_derivative)
    c_em = sp.sympify(wave_speed)
    radius = sp.sympify(major_radius)
    if c_em.is_zero is True or radius.is_zero is True:
        raise ValueError("wave_speed and major_radius must be nonzero")
    return sp.simplify(-derivative / (c_em**2 * radius))


def ampere_eliminated_lock_response(
    tag,
    tag_derivative,
    phi_hat,
    phi_hat_radial_derivative,
    tag_primitive,
    permittivity,
    speed,
    wave_speed,
    radius,
):
    """Eliminate the toroidal magnetic field from the lock response.

    Dividing the exact Ampere primitive by signed ``g`` gives

    ``epsilon*a*c_EM**2*r*H_hat + epsilon*c*d_r Phi_hat = K(P)/r``.

    The returned expression is exactly ``chi'*Phi_hat-chi*H_hat``.  No
    contour projection is applied here; in a curved ring even the explicit
    ``K(P)/r**2`` row need not be a streamline function.
    """

    chi = sp.sympify(tag)
    chip = sp.sympify(tag_derivative)
    phi = sp.sympify(phi_hat)
    phi_r = sp.sympify(phi_hat_radial_derivative)
    primitive = sp.sympify(tag_primitive)
    eps = sp.sympify(permittivity)
    c = sp.sympify(speed)
    c_em = sp.sympify(wave_speed)
    r = sp.sympify(radius)
    if eps.is_zero is True or c_em.is_zero is True or r.is_zero is True:
        raise ValueError("permittivity, wave_speed, and radius must be nonzero")
    a = sp.simplify(1 - c**2 / c_em**2)
    if a.is_nonpositive is True:
        raise ValueError("the traveling Maxwell block must be strictly elliptic")
    magnetic = primitive / (eps * a * c_em**2 * r**2) - c * phi_r / (
        a * c_em**2 * r
    )
    return sp.simplify(chip * phi - chi * magnetic)


def annular_dipole_matching_shift(inner_second_moment, outer_zeroth_moment):
    """Return the ``(A,B)`` shift in an annular ``A*r+B/r`` dipole row.

    For ``L_1 v=f`` with the regular/decaying Green representation

    ``v(r)=-(r**(-1)*integral_0^r t**2*f(t)dt
             +r*integral_r^infinity f(t)dt)/2``,

    a source perturbation supported strictly inside an observation annulus
    changes only ``B``; one supported strictly outside changes only ``A``.
    The arguments are precisely those two moments.
    """

    inner = sp.sympify(inner_second_moment)
    outer = sp.sympify(outer_zeroth_moment)
    return sp.simplify(-outer / 2), sp.simplify(-inner / 2)
