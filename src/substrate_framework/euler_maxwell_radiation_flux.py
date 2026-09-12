"""Exact Maxwell shell-to-flux identities for a translating conserved current.

The helpers use the unitary spatial Fourier transform and real fields written
as ``Re(exp(-i*omega*t) profile(x-c*t*e_z))``.  They describe prescribed
current radiation in the explicit Euler--Maxwell extension.  They do not
construct a Cao response current, a resonance pole, or a particle.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


def _require_positive(name: str, value: sp.Expr) -> None:
    if value.is_positive is False or value.is_zero is True:
        raise ValueError(f"{name} must be positive")


def _require_subluminal(maxwell_speed: sp.Expr, carrier_speed: sp.Expr) -> None:
    gap = sp.simplify(maxwell_speed**2 - carrier_speed**2)
    if gap.is_positive is False or gap.is_zero is True:
        raise ValueError("carrier_speed must be strictly subluminal")


def _require_direction_component(direction_z: sp.Expr) -> None:
    excess = sp.simplify(sp.Abs(direction_z) - 1)
    if excess.is_positive is True:
        raise ValueError("direction_z must lie in the unit interval")


@dataclass(frozen=True)
class DopplerShellGeometry:
    """Star-shaped positive-radius component of a comoving Maxwell shell."""

    radius: sp.Expr
    temporal_frequency: sp.Expr
    radial_derivative_magnitude: sp.Expr
    gradient_squared: sp.Expr


def doppler_shell_geometry(
    maxwell_speed,
    carrier_speed,
    frequency,
    direction_z,
) -> DopplerShellGeometry:
    """Return exact geometry of ``c^2|k|^2-(omega+v*k_z)^2=0``.

    ``direction_z`` is the z component of a unit direction.  The returned
    shell has ``sign(Omega)=sign(omega)``.
    """

    c_em = sp.sympify(maxwell_speed)
    speed = sp.sympify(carrier_speed)
    omega = sp.sympify(frequency)
    n_z = sp.sympify(direction_z)
    _require_positive("maxwell_speed", c_em)
    if omega.is_zero is True:
        raise ValueError("frequency must be nonzero")
    _require_subluminal(c_em, speed)
    _require_direction_component(n_z)
    sign = sp.sign(omega)
    radius = sp.factor(sp.Abs(omega) / (c_em - sign * speed * n_z))
    omega_shell = sp.factor(sign * c_em * radius)
    radial = sp.factor(2 * c_em * sp.Abs(omega))
    gradient_squared = sp.factor(
        4
        * c_em**2
        * radius**2
        * (c_em**2 + speed**2 - 2 * sign * c_em * speed * n_z)
    )
    return DopplerShellGeometry(
        radius=radius,
        temporal_frequency=omega_shell,
        radial_derivative_magnitude=radial,
        gradient_squared=gradient_squared,
    )


def outgoing_power_prefactor(permittivity):
    """Return ``pi/(2*epsilon)`` in the outgoing shell power identity."""

    epsilon = sp.sympify(permittivity)
    _require_positive("permittivity", epsilon)
    return sp.pi / (2 * epsilon)


def shell_sphere_power_weight(
    permittivity,
    maxwell_speed,
    carrier_speed,
    frequency,
    direction_z,
):
    """Return the coefficient multiplying ``|P_T J_hat(r*n)|^2 dOmega``.

    This is

    ``pi*|omega|^2 / (4*epsilon*(c-sign(omega)*v*n_z)^3)``.
    """

    epsilon = sp.sympify(permittivity)
    c_em = sp.sympify(maxwell_speed)
    speed = sp.sympify(carrier_speed)
    omega = sp.sympify(frequency)
    n_z = sp.sympify(direction_z)
    _require_positive("permittivity", epsilon)
    _require_positive("maxwell_speed", c_em)
    if omega.is_zero is True:
        raise ValueError("frequency must be nonzero")
    _require_subluminal(c_em, speed)
    _require_direction_component(n_z)
    return sp.factor(
        sp.pi
        * sp.Abs(omega) ** 2
        / (4 * epsilon * (c_em - sp.sign(omega) * speed * n_z) ** 3)
    )


def switched_continuity_residual(
    envelope,
    envelope_derivative,
    mode_residual,
    charge_profile,
    completion_divergence,
):
    """Return the scalar residual of the conserved envelope completion.

    For ``rho_T=a*rho`` and ``J_T=a*J+a'*K`` in translating phasor
    variables, the residual is ``a*R_mode+a'*(rho+div K)``.
    """

    a = sp.sympify(envelope)
    da = sp.sympify(envelope_derivative)
    residual = sp.sympify(mode_residual)
    rho = sp.sympify(charge_profile)
    div_k = sp.sympify(completion_divergence)
    return sp.expand(a * residual + da * (rho + div_k))


def gaussian_curl_current_power(
    permittivity,
    maxwell_speed,
    frequency_magnitude,
    gaussian_width,
    current_scale_squared,
    polarization_norm_squared,
):
    """Exact outgoing power of a transverse Gaussian curl current at rest.

    The Fourier profile is
    ``J_hat=i*j0*(k cross a)*exp(-sigma**2*|k|**2/2)`` with ``c_g=0``.
    The final two arguments are ``j0**2`` and ``|a|**2``.
    """

    epsilon = sp.sympify(permittivity)
    c_em = sp.sympify(maxwell_speed)
    omega = sp.sympify(frequency_magnitude)
    width = sp.sympify(gaussian_width)
    j_sq = sp.sympify(current_scale_squared)
    a_sq = sp.sympify(polarization_norm_squared)
    for name, value in (
        ("permittivity", epsilon),
        ("maxwell_speed", c_em),
        ("frequency_magnitude", omega),
        ("gaussian_width", width),
        ("current_scale_squared", j_sq),
        ("polarization_norm_squared", a_sq),
    ):
        _require_positive(name, value)
    return sp.factor(
        2
        * sp.pi**2
        * j_sq
        * a_sq
        * omega**4
        * sp.exp(-(width * omega / c_em) ** 2)
        / (3 * epsilon * c_em**5)
    )
