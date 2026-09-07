"""Exact scaling identities for high-harmonic Cao-mode radiation.

The helpers expose the physical core clock, the sharp transverse Maxwell
shell radius, and the large-order Bessel ratios used by P253/0087.  They do
not construct a Cao mode, prove a dark-current condition, or turn a shell
trace into a resonance width.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


def _require_positive(name: str, value: sp.Expr) -> None:
    if value.is_positive is False or value.is_zero is True:
        raise ValueError(f"{name} must be positive")


def _require_subluminal(maxwell_speed: sp.Expr, carrier_speed: sp.Expr) -> sp.Expr:
    gap = sp.simplify(maxwell_speed**2 - carrier_speed**2)
    if gap.is_positive is False or gap.is_zero is True:
        raise ValueError("carrier_speed must be strictly subluminal")
    return gap


@dataclass(frozen=True)
class PhysicalModeLedger:
    """Core clock and physical frequency for one toroidal mode."""

    core_radius: sp.Expr
    core_clock: sp.Expr
    physical_frequency: sp.Expr


def physical_mode_ledger(circulation, radius, wave_number, harmonic, frequency_magnitude):
    """Return ``a=R*k/n``, ``Omega_a=kappa/(2*pi*a**2)``, and ``Omega_a*sigma``.

    The circulation convention is ``kappa=Gamma``.  ``frequency_magnitude``
    is the dimensionless positive column/ring frequency ``|sigma|``.
    """

    kappa = sp.sympify(circulation)
    ring_radius = sp.sympify(radius)
    k = sp.sympify(wave_number)
    n = sp.sympify(harmonic)
    sigma = sp.sympify(frequency_magnitude)
    for name, value in (
        ("circulation", kappa),
        ("radius", ring_radius),
        ("wave_number", k),
        ("harmonic", n),
        ("frequency_magnitude", sigma),
    ):
        _require_positive(name, value)
    core_radius = sp.factor(ring_radius * k / n)
    core_clock = sp.factor(kappa / (2 * sp.pi * core_radius**2))
    return PhysicalModeLedger(
        core_radius=core_radius,
        core_clock=core_clock,
        physical_frequency=sp.factor(core_clock * sigma),
    )


@dataclass(frozen=True)
class TransverseShellMaximum:
    """Sharp transverse radius of a translating Maxwell radiation shell."""

    maximizing_direction_z: sp.Expr
    transverse_wave_number: sp.Expr
    conservative_wave_number: sp.Expr


def transverse_shell_maximum(maxwell_speed, carrier_speed, frequency_magnitude, frequency_sign=1):
    """Return the sharp maximum of ``k_perp`` on the real Doppler shell.

    For ``|k|**2=(omega+c*k_z)**2/c_EM**2`` and
    ``t=k_z/|k|``, the maximum occurs at
    ``t=sign(omega)*c/c_EM`` and equals
    ``|omega|/sqrt(c_EM**2-c**2)``.  The final field is the weaker bound
    ``|omega|/(c_EM-|c|)``.
    """

    c_em = sp.sympify(maxwell_speed)
    c = sp.sympify(carrier_speed)
    omega = sp.sympify(frequency_magnitude)
    sign = sp.sympify(frequency_sign)
    _require_positive("maxwell_speed", c_em)
    _require_positive("frequency_magnitude", omega)
    if sign not in (sp.Integer(-1), sp.Integer(1)):
        raise ValueError("frequency_sign must be +1 or -1")
    _require_subluminal(c_em, c)
    return TransverseShellMaximum(
        maximizing_direction_z=sp.factor(sign * c / c_em),
        transverse_wave_number=sp.factor(omega / sp.sqrt(c_em**2 - c**2)),
        conservative_wave_number=sp.factor(omega / (c_em - sp.Abs(c))),
    )


def fixed_mode_bessel_ratio(
    circulation,
    frequency_magnitude,
    harmonic,
    wave_number,
    radius,
    maxwell_speed,
    carrier_speed,
    support_constant=1,
):
    """Return ``k_perp,max*(R+C*a)/(n-1)`` for one massive fixed-J mode."""

    ledger = physical_mode_ledger(
        circulation, radius, wave_number, harmonic, frequency_magnitude
    )
    n = sp.sympify(harmonic)
    if (n - 1).is_positive is False:
        raise ValueError("harmonic must exceed one")
    ring_radius = sp.sympify(radius)
    support_factor = sp.sympify(support_constant)
    if support_factor.is_nonnegative is False:
        raise ValueError("support_constant must be nonnegative")
    support = ring_radius + support_factor * ledger.core_radius
    shell = transverse_shell_maximum(
        maxwell_speed, carrier_speed, ledger.physical_frequency
    )
    return sp.factor(shell.transverse_wave_number * support / (n - 1))


def high_index_limiting_bessel_ratio(
    circulation,
    harmonic_step,
    column_coefficient,
    radius,
    wave_number,
    maxwell_speed,
    carrier_speed,
):
    """Return the finite-J leading predictor at the supplied carrier speed."""

    kappa = sp.sympify(circulation)
    ell = sp.sympify(harmonic_step)
    l_phi = sp.sympify(column_coefficient)
    ring_radius = sp.sympify(radius)
    k = sp.sympify(wave_number)
    c_em = sp.sympify(maxwell_speed)
    c = sp.sympify(carrier_speed)
    for name, value in (
        ("circulation", kappa),
        ("harmonic_step", ell),
        ("column_coefficient", l_phi),
        ("radius", ring_radius),
        ("wave_number", k),
        ("maxwell_speed", c_em),
    ):
        _require_positive(name, value)
    _require_subluminal(c_em, c)
    return sp.factor(
        kappa
        * ell
        * l_phi
        / (2 * sp.pi**2 * ring_radius * k * sp.sqrt(c_em**2 - c**2))
    )


def debye_speed_ceiling(maxwell_speed, zero_speed_ratio_numerator, margin):
    """Return the carrier-speed ceiling for ``q_rad <= 1-margin``.

    The public campaign helper requires the strict Debye margin
    ``0 < margin < 1``.  The mathematical boundary ``margin=0`` belongs to a
    lower-level turning-point analysis and is intentionally outside this API.
    """

    c_em = sp.sympify(maxwell_speed)
    c_0 = sp.sympify(zero_speed_ratio_numerator)
    eta = sp.sympify(margin)
    _require_positive("maxwell_speed", c_em)
    _require_positive("zero_speed_ratio_numerator", c_0)
    if eta.is_positive is False or (1 - eta).is_positive is False:
        raise ValueError("margin must lie strictly between zero and one")
    radicand = sp.simplify(c_em**2 - (c_0 / (1 - eta)) ** 2)
    if radicand.is_positive is False or radicand.is_zero is True:
        raise ValueError("Debye margin is impossible at the supplied Maxwell speed")
    return sp.factor(sp.sqrt(radicand))


def debye_exponential_rate(argument_order_ratio):
    """Return the optimized contour rate for ``0<q<1``.

    The generating-function contour bound gives
    ``|J_m(m*q)| <= exp(-m*(acosh(1/q)-sqrt(1-q**2)))``.
    """

    q = sp.sympify(argument_order_ratio)
    if q.is_positive is False or (1 - q).is_positive is False:
        raise ValueError("argument/order ratio must lie strictly between zero and one")
    return sp.factor(sp.acosh(1 / q) - sp.sqrt(1 - q**2))


def cylindrical_vector_bessel_orders(harmonic):
    """Return the three Fourier-Bessel orders from a cylindrical vector mode."""

    n = sp.sympify(harmonic)
    return (n - 1, n, n + 1)
