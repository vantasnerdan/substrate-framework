"""Exact quadratic identities for a tagged Euler--Maxwell carrier.

The helpers describe the declared classical ``U(1)`` extension used by the
P253 campaign.  They do not prove stability, construct a charged carrier, or
derive the gauge constants from incompressible Euler.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ComovingMaxwellQuadratic:
    """Two algebraically equal forms of ``E_EM-c P_EM,z``."""

    relative_density: sp.Expr
    completed_density: sp.Expr
    speed_ratio_squared: sp.Expr


def comoving_maxwell_quadratic(permittivity, permeability, speed, electric, magnetic):
    """Return the exact comoving Maxwell quadratic density.

    Field momentum uses ``P_EM = epsilon * E cross B`` and
    ``c_EM**2 = 1/(epsilon*mu)``.  Positivity of the completed form requires
    ``speed**2 < c_EM**2``.
    """

    epsilon = sp.sympify(permittivity)
    mu = sp.sympify(permeability)
    c = sp.sympify(speed)
    if epsilon.is_positive is False or epsilon.is_zero is True:
        raise ValueError("permittivity must be positive")
    if mu.is_positive is False or mu.is_zero is True:
        raise ValueError("permeability must be positive")
    e = sp.ImmutableMatrix(electric)
    b = sp.ImmutableMatrix(magnetic)
    if e.shape != (3, 1) or b.shape != (3, 1):
        raise ValueError("electric and magnetic fields must have three components")

    e_z = sp.ImmutableMatrix([0, 0, 1])
    b_perp = sp.ImmutableMatrix([b[0], b[1], 0])
    relative = (
        epsilon * e.dot(e) / 2
        + b.dot(b) / (2 * mu)
        - c * epsilon * e.cross(b)[2]
    )
    beta_squared = sp.factor(epsilon * mu * c**2)
    completed = (
        epsilon * (e + c * e_z.cross(b)).dot(e + c * e_z.cross(b)) / 2
        + (1 - beta_squared) * b_perp.dot(b_perp) / (2 * mu)
        + b[2] ** 2 / (2 * mu)
    )
    return ComovingMaxwellQuadratic(
        relative_density=sp.expand(relative),
        completed_density=sp.expand(completed),
        speed_ratio_squared=beta_squared,
    )


def constrained_charge_schur_symbol(
    permittivity, permeability, speed, wavevector, charge_amplitude_squared
):
    """Return the minimized field quadratic for one Fourier charge mode.

    The minimization is over the transverse electric and divergence-free
    magnetic components at fixed Gauss row ``i k dot E = rho/epsilon``.
    The result is the anisotropic positive ``H^{-1}`` multiplier

    ``|rho|^2 (1-beta^2) / (2 epsilon (|k|^2-beta^2 k_z^2))``.

    It is a field-only constrained Schur complement.  Matter-current and
    coadjoint constraints must still be included in a carrier Hessian.
    """

    epsilon = sp.sympify(permittivity)
    mu = sp.sympify(permeability)
    c = sp.sympify(speed)
    if epsilon.is_positive is False or epsilon.is_zero is True:
        raise ValueError("permittivity must be positive")
    if mu.is_positive is False or mu.is_zero is True:
        raise ValueError("permeability must be positive")
    k = sp.ImmutableMatrix(wavevector)
    if k.shape != (3, 1):
        raise ValueError("wavevector must have three components")
    k_squared = sp.factor(k.dot(k))
    if k_squared == 0:
        raise ValueError("wavevector must be nonzero")
    beta_squared = sp.factor(epsilon * mu * c**2)
    denominator = sp.factor(k_squared - beta_squared * k[2] ** 2)
    return sp.factor(
        sp.sympify(charge_amplitude_squared)
        * (1 - beta_squared)
        / (2 * epsilon * denominator)
    )


@dataclass(frozen=True)
class MaterialTagLocking:
    """Modewise tag-to-vorticity coefficients on a regular toroidal band."""

    axisymmetric_from_toroidal_vorticity: sp.Expr
    nonaxisymmetric_from_radial_vorticity: sp.Expr


def material_tag_locking_coefficients(tag_derivative, vorticity, vorticity_derivative, harmonic):
    """Return exact coefficients for ``chi=F(I)``, ``Omega=zeta(I)d_theta``.

    For ``delta Omega=-[xi,Omega]`` and ``delta chi=-xi^I F'``, the
    axisymmetric relation is ``delta chi=(F'/zeta') delta Omega^theta``.
    For a nonzero toroidal Fourier harmonic ``n``, it is
    ``delta chi=-(F'/(i*n*zeta)) delta Omega^I``.
    """

    f_prime = sp.sympify(tag_derivative)
    zeta = sp.sympify(vorticity)
    zeta_prime = sp.sympify(vorticity_derivative)
    n = sp.sympify(harmonic)
    if zeta.is_zero is True:
        raise ValueError("vorticity must be nonzero on the tagged band")
    if zeta_prime.is_zero is True:
        raise ValueError("vorticity derivative must be nonzero on the tagged band")
    if n.is_zero is True:
        nonaxisymmetric = sp.nan
    else:
        nonaxisymmetric = sp.factor(-f_prime / (sp.I * n * zeta))
    return MaterialTagLocking(
        axisymmetric_from_toroidal_vorticity=sp.factor(f_prime / zeta_prime),
        nonaxisymmetric_from_radial_vorticity=nonaxisymmetric,
    )


@dataclass(frozen=True)
class ComovingRadiationShell:
    """Radial root and coarea data for one signed nonzero frequency."""

    radius: sp.Expr
    radial_derivative_magnitude: sp.Expr
    coarea_weight: sp.Expr


def comoving_radiation_shell(
    maxwell_speed, carrier_speed, frequency_magnitude, frequency_sign, direction_z
):
    """Return the exact star-shaped Maxwell shell in the translating frame.

    ``frequency_sign`` must be ``+1`` or ``-1`` and
    ``frequency_magnitude`` must be positive.  For a unit direction ``n`` the
    positive radial root is

    ``r=|omega|/(c_EM-sign(omega)*c*n_z)``.

    The returned coarea weight is ``r**2/|partial_r D|`` for
    ``D=|k|**2-(omega+c*k_z)**2/c_EM**2``.
    """

    c_em = sp.sympify(maxwell_speed)
    c = sp.sympify(carrier_speed)
    omega_abs = sp.sympify(frequency_magnitude)
    sign = sp.sympify(frequency_sign)
    n_z = sp.sympify(direction_z)
    if c_em.is_positive is False or c_em.is_zero is True:
        raise ValueError("maxwell_speed must be positive")
    if omega_abs.is_positive is False or omega_abs.is_zero is True:
        raise ValueError("frequency_magnitude must be positive")
    if sign not in (sp.Integer(-1), sp.Integer(1)):
        raise ValueError("frequency_sign must be +1 or -1")
    denominator = c_em - sign * c * n_z
    radius = sp.factor(omega_abs / denominator)
    radial_derivative = sp.factor(2 * omega_abs / c_em)
    return ComovingRadiationShell(
        radius=radius,
        radial_derivative_magnitude=radial_derivative,
        coarea_weight=sp.factor(radius**2 / radial_derivative),
    )


def comoving_maxwell_radiation_denominator(maxwell_speed, carrier_speed, frequency, wavevector):
    """Return the traveling-frame Maxwell wave denominator.

    A mode ``exp(i k dot (x-c t e_z)-i omega t)`` has laboratory frequency
    ``omega+c*k_z``.  Its free Maxwell denominator is therefore
    ``|k|^2-(omega+c*k_z)^2/c_EM^2``.  Zeros are the radiation shell on
    which a localized periodic source needs a transverse-current
    cancellation to remain nonradiating.
    """

    c_em = sp.sympify(maxwell_speed)
    c = sp.sympify(carrier_speed)
    omega = sp.sympify(frequency)
    if c_em.is_positive is False or c_em.is_zero is True:
        raise ValueError("maxwell_speed must be positive")
    k = sp.ImmutableMatrix(wavevector)
    if k.shape != (3, 1):
        raise ValueError("wavevector must have three components")
    return sp.factor(k.dot(k) - (omega + c * k[2]) ** 2 / c_em**2)
