"""Exact normalized whole-line overlap and conditional mass-map ledgers.

All integrals in this module use Cartesian Lebesgue measure ``dx`` on the real
line.  The hyperbolic-secant formulas require one declared common inverse width
for the mode and multiplier profile.  These mathematical expectations do not
derive a Yukawa interaction, a generation assignment, or an absolute mass.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and (
        result.is_real is not True or result.is_positive is not True
    ):
        raise ValueError(f"{name} must be positive and real")
    return result


def _real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and result.is_real is not True:
        raise ValueError(f"{name} must be real")
    return result


def _nonzero_real(value: Any, name: str) -> sp.Expr:
    result = _real(value, name)
    if result.is_number and sp.simplify(result) == 0:
        raise ValueError(f"{name} must be nonzero")
    return result


@dataclass(frozen=True)
class MatchedSechOverlap:
    """Exact normalized overlap of matched-width sech-power profiles."""

    mode_power: sp.Expr
    profile_power: sp.Expr
    amplitude: sp.Expr
    inverse_width: sp.Expr
    unnormalized_mode_norm: sp.Expr
    raw_overlap: sp.Expr
    normalized_overlap: sp.Expr


@dataclass(frozen=True)
class QuarticBoundModeOverlapLedger:
    """Normalized squared-density expectations for C-QBL-003's two modes."""

    even_mode_norm: sp.Expr
    odd_mode_norm: sp.Expr
    even_overlap: sp.Expr
    odd_overlap: sp.Expr
    weighted_cross_overlap: sp.Expr


@dataclass(frozen=True)
class ConditionalOverlapMassLedger:
    """A declared product map and its additive mass-dimension bookkeeping."""

    overlap: sp.Expr
    scale: sp.Expr
    mapped_mass: sp.Expr
    profile_mass_dimension: sp.Expr
    scale_mass_dimension: sp.Expr
    mapped_mass_dimension: sp.Expr


def normalized_expectation_bounds(
    profile_lower: Any,
    profile_upper: Any,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the sharp multiplier bounds for a normalized density.

    If ``integral |eta|**2 dx = 1`` and the supplied real multiplication
    profile obeys ``lower <= Phi(x) <= upper`` almost everywhere, then
    ``lower <= integral |eta|**2 Phi dx <= upper``.  This function records the
    theorem's supplied essential-range endpoints; callers remain responsible
    for proving the profile bound and normalization.
    """

    lower = _real(profile_lower, "profile_lower")
    upper = _real(profile_upper, "profile_upper")
    if lower.is_number and upper.is_number and bool(lower > upper):
        raise ValueError("profile_lower must not exceed profile_upper")
    return lower, upper


def sech_power_integral(power: Any, inverse_width: Any) -> sp.Expr:
    """Return ``integral_R sech(kappa*x)**s dx`` for ``s,kappa>0``.

    The exact beta-function value is
    ``sqrt(pi)*gamma(s/2)/(kappa*gamma((s+1)/2))``.
    """

    exponent = _positive_real(power, "power")
    kappa = _positive_real(inverse_width, "inverse_width")
    return sp.sqrt(sp.pi) * sp.gamma(exponent / 2) / (
        kappa * sp.gamma((exponent + 1) / 2)
    )


def matched_width_sech_overlap(
    mode_power: Any,
    profile_power: Any,
    amplitude: Any,
    inverse_width: Any,
) -> MatchedSechOverlap:
    """Return the exact normalized ``sech**p`` density overlap.

    For a real mode proportional to ``sech(kappa*x)**p`` and multiplier
    ``A*sech(kappa*x)**r``, the normalized expectation is

    ``A*I(2*p+r,kappa)/I(2*p,kappa)``.

    The common inverse width cancels from the normalized ratio.  This does not
    cover mismatched widths or a different integration measure.
    """

    p = _positive_real(mode_power, "mode_power")
    r = _positive_real(profile_power, "profile_power")
    coefficient = _real(amplitude, "amplitude")
    kappa = _positive_real(inverse_width, "inverse_width")
    norm = sech_power_integral(2 * p, kappa)
    raw = sp.simplify(coefficient * sech_power_integral(2 * p + r, kappa))
    return MatchedSechOverlap(
        mode_power=p,
        profile_power=r,
        amplitude=coefficient,
        inverse_width=kappa,
        unnormalized_mode_norm=norm,
        raw_overlap=raw,
        normalized_overlap=sp.simplify(raw / norm),
    )


def quartic_bound_mode_overlap_ledger(
    amplitude: Any,
    inverse_width: Any,
) -> QuarticBoundModeOverlapLedger:
    """Return exact overlaps for C-QBL-003's even and odd mode shapes.

    The unnormalized shapes are ``sech(z)**2`` and
    ``sech(z)*tanh(z)`` with ``z=kappa*x``.  The declared multiplier is
    ``A*sech(z)``.  Its even parity makes the normalized weighted cross overlap
    vanish, while the two squared-density expectations are respectively
    ``9*pi*A/32`` and ``3*pi*A/16``.
    """

    coefficient = _real(amplitude, "amplitude")
    kappa = _positive_real(inverse_width, "inverse_width")
    integral_2 = sech_power_integral(2, kappa)
    integral_3 = sech_power_integral(3, kappa)
    integral_4 = sech_power_integral(4, kappa)
    integral_5 = sech_power_integral(5, kappa)
    even_norm = integral_4
    odd_norm = sp.simplify(integral_2 - integral_4)
    even_overlap = sp.simplify(coefficient * integral_5 / even_norm)
    odd_overlap = sp.simplify(
        coefficient * (integral_3 - integral_5) / odd_norm
    )
    return QuarticBoundModeOverlapLedger(
        even_mode_norm=even_norm,
        odd_mode_norm=odd_norm,
        even_overlap=even_overlap,
        odd_overlap=odd_overlap,
        weighted_cross_overlap=sp.Integer(0),
    )


def conditional_overlap_mass_ledger(
    overlap: Any,
    scale: Any,
    *,
    profile_mass_dimension: Any = 0,
    scale_mass_dimension: Any = 1,
) -> ConditionalOverlapMassLedger:
    """Return a separately declared product ``m=y*v`` and its dimensions.

    L2 normalization makes the overlap ``y`` carry the multiplier profile's
    mass dimension.  The product has the sum of that dimension and the
    supplied scale dimension.  The function checks bookkeeping only; it does
    not derive an interaction or identify the product with a physical mass.
    """

    expectation = _real(overlap, "overlap")
    free_scale = _real(scale, "scale")
    profile_dimension = _real(profile_mass_dimension, "profile_mass_dimension")
    scale_dimension = _real(scale_mass_dimension, "scale_mass_dimension")
    return ConditionalOverlapMassLedger(
        overlap=expectation,
        scale=free_scale,
        mapped_mass=sp.simplify(expectation * free_scale),
        profile_mass_dimension=profile_dimension,
        scale_mass_dimension=scale_dimension,
        mapped_mass_dimension=sp.simplify(profile_dimension + scale_dimension),
    )


def overlap_mass_ratio(overlap: Any, reference_overlap: Any) -> sp.Expr:
    """Return the ratio for products sharing the same nonzero free scale."""

    value = _real(overlap, "overlap")
    reference = _nonzero_real(reference_overlap, "reference_overlap")
    return sp.simplify(value / reference)


def reciprocal_overlap_scale_rescaling(
    overlap: Any,
    scale: Any,
    factor: Any,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return ``(rho*y, v/rho, y*v)`` for a nonzero rescaling ``rho``.

    The invariant product exhibits one free rescaling direction whenever the
    overlap amplitude and external scale are not independently fixed.
    """

    expectation = _real(overlap, "overlap")
    free_scale = _real(scale, "scale")
    rho = _nonzero_real(factor, "factor")
    return (
        sp.simplify(rho * expectation),
        sp.simplify(free_scale / rho),
        sp.simplify(expectation * free_scale),
    )
