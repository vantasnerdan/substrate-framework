"""Exact conditional continuum density-of-states and mode-counting formulas.

The functions in this module use the explicitly supplied phase-space measure
``V*d**d k/(2*pi)**d`` and an isotropic gapped dispersion.  Spatial dimension,
branch degeneracy, volume, and cutoff are independent inputs.  A continuum
ball volume is not an exact finite lattice-point count, and matching it to a
target count does not derive a microscopic Brillouin zone, material spectrum,
participating-mode set, coupling, occupation, or rate.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive_integer(value: Any, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be a positive integer")
    return int(expression)


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    if expression.is_positive is False:
        raise ValueError(f"{name} must not be known nonpositive")
    return expression


def _nonnegative(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    if expression.is_nonnegative is False:
        raise ValueError(f"{name} must not be known negative")
    return expression


def unit_sphere_surface(dimension: Any) -> sp.Expr:
    r"""Return the surface area of the unit sphere in ``dimension`` dimensions.

    For positive integer ``d`` the result is
    ``2*pi**(d/2)/gamma(d/2)``.  Thus the first three values are ``2``,
    ``2*pi``, and ``4*pi``.
    """

    d = _positive_integer(dimension, "dimension")
    return sp.simplify(
        2 * sp.pi ** sp.Rational(d, 2) / sp.gamma(sp.Rational(d, 2))
    )


def isotropic_gapped_angular_frequency(
    wavenumber: Any,
    signal_speed: Any,
    gap_frequency: Any,
) -> sp.Expr:
    r"""Return ``sqrt(omega_0**2+c**2*k**2)`` for nonnegative radial ``k``."""

    k = _nonnegative(wavenumber, "wavenumber")
    c = _positive(signal_speed, "signal_speed")
    omega_0 = _nonnegative(gap_frequency, "gap_frequency")
    return sp.sqrt(omega_0**2 + c**2 * k**2)


def isotropic_continuum_dos_on_band(
    angular_frequency: Any,
    volume: Any,
    signal_speed: Any,
    gap_frequency: Any,
    dimension: Any,
    *,
    branches: Any = 1,
) -> sp.Expr:
    r"""Return the continuum frequency DOS on the open propagating band.

    For separately supplied positive integer spatial dimension ``d`` and
    branch degeneracy ``b``, positive ``V`` and ``c``, nonnegative
    ``omega_0``, and ``omega>omega_0``, the result is

    ``b*V*S_(d-1)*omega*(omega**2-omega_0**2)**((d-2)/2)``
    ``/ ((2*pi)**d*c**d)``.

    The caller owns the open-band condition for symbolic inputs.  The
    ``d=1`` expression has an integrable threshold divergence, so this helper
    deliberately does not invent a point value at the band edge.  The DOS is
    zero below the gap by definition, but this on-band helper returns no
    piecewise extension.
    """

    omega = _nonnegative(angular_frequency, "angular_frequency")
    V = _positive(volume, "volume")
    c = _positive(signal_speed, "signal_speed")
    omega_0 = _nonnegative(gap_frequency, "gap_frequency")
    d = _positive_integer(dimension, "dimension")
    b = _positive_integer(branches, "branches")
    if omega.is_number and omega_0.is_number and (omega > omega_0) is not sp.S.true:
        raise ValueError("angular_frequency must exceed gap_frequency on the open band")
    return sp.simplify(
        b
        * V
        * unit_sphere_surface(d)
        * omega
        * (omega**2 - omega_0**2) ** sp.Rational(d - 2, 2)
        / ((2 * sp.pi) ** d * c**d)
    )


def isotropic_continuum_mode_count(
    radial_cutoff: Any,
    volume: Any,
    dimension: Any,
    *,
    branches: Any = 1,
) -> sp.Expr:
    r"""Return the supplied-measure count inside a radial momentum ball.

    The exact phase-space expression is
    ``b*V*S_(d-1)*K**d/(d*(2*pi)**d)``.  It is exact for the declared
    continuum measure, not an exact count of discrete wavevectors in a finite
    box.
    """

    cutoff = _positive(radial_cutoff, "radial_cutoff")
    V = _positive(volume, "volume")
    d = _positive_integer(dimension, "dimension")
    b = _positive_integer(branches, "branches")
    return sp.simplify(
        b * V * unit_sphere_surface(d) * cutoff**d / (d * (2 * sp.pi) ** d)
    )


def isotropic_continuum_target_cutoff(
    target_count: Any,
    volume: Any,
    dimension: Any,
    *,
    branches: Any = 1,
) -> sp.Expr:
    r"""Return the radial cutoff that matches a supplied continuum target.

    Solving :func:`isotropic_continuum_mode_count` for positive ``K`` gives
    ``(N*d*(2*pi)**d/(b*V*S_(d-1)))**(1/d)``.  The target is an input; this
    inversion does not derive it or a microscopic cutoff geometry.
    """

    target = _positive(target_count, "target_count")
    V = _positive(volume, "volume")
    d = _positive_integer(dimension, "dimension")
    b = _positive_integer(branches, "branches")
    return sp.simplify(
        (
            target
            * d
            * (2 * sp.pi) ** d
            / (b * V * unit_sphere_surface(d))
        )
        ** sp.Rational(1, d)
    )

