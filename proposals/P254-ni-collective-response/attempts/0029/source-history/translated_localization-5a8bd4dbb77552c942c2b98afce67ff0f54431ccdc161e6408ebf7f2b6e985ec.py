"""Exact translated-localization overlaps and identifiability ledgers.

The closed-form overlaps in this module use Cartesian Lebesgue measure ``dx``
on the real line.  A translated one-dimensional well is not thereby a radial
mode or a generation, and a supplied center spacing remains a free input.
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


@dataclass(frozen=True)
class TranslatedSechOverlap:
    """Exact normalized overlap for matched-width translated sech profiles."""

    mode_power: sp.Expr
    profile_power: sp.Expr
    amplitude: sp.Expr
    inverse_width: sp.Expr
    displacement: sp.Expr
    dimensionless_displacement: sp.Expr
    mode_density_norm: sp.Expr
    raw_dimensionless_convolution: sp.Expr
    normalized_overlap: sp.Expr


@dataclass(frozen=True)
class SechTailLedger:
    """Large-displacement class for an exact sech-density convolution."""

    mode_density_power: sp.Expr
    profile_power: sp.Expr
    dimensionless_decay_power: sp.Expr
    physical_decay_rate: sp.Expr
    polynomial_prefactor_power: int
    normalized_leading_coefficient: sp.Expr


@dataclass(frozen=True)
class PoschlTellerGroundLedger:
    """Exact ground state of one declared translated Pöschl--Teller well."""

    depth: sp.Expr
    width: sp.Expr
    center: sp.Expr
    index: sp.Expr
    eigenvalue: sp.Expr
    normalization: sp.Expr
    density_tail_rate: sp.Expr


@dataclass(frozen=True)
class PoschlSechTailLedger:
    """Exact tail class for a Pöschl density against a fixed sech core."""

    ground: PoschlTellerGroundLedger
    profile_amplitude: sp.Expr
    profile_inverse_width: sp.Expr
    overlap_tail_rate: sp.Expr
    polynomial_prefactor_power: int
    leading_coefficient: sp.Expr


@dataclass(frozen=True)
class TailSpacingLedger:
    """Conditional asymptotic ladder slope and its free product inputs."""

    mode_density_tail_rate: sp.Expr
    profile_tail_rate: sp.Expr
    spacing: sp.Expr
    overlap_tail_rate: sp.Expr
    asymptotic_log_ratio: sp.Expr
    resonant_equal_rates: bool | None


def sech_convolution(
    shifted_power: Any,
    core_power: Any,
    dimensionless_displacement: Any,
) -> sp.Expr:
    """Return ``integral sech(z-a)^alpha sech(z)^beta dz`` exactly.

    ``alpha`` and ``beta`` must be positive.  The displacement may be any real
    value; reflection symmetry makes the result depend on ``|a|``.  The
    hypergeometric representation is

    ``2^(alpha+beta-1) exp(-alpha|a|) B(l,l)
      2F1(alpha,l;2l;1-exp(-2|a|))``,

    where ``l=(alpha+beta)/2``.
    """

    alpha = _positive_real(shifted_power, "shifted_power")
    beta = _positive_real(core_power, "core_power")
    displacement = sp.Abs(_real(dimensionless_displacement, "displacement"))
    half_sum = (alpha + beta) / 2
    if displacement.is_zero is True:
        return sp.sqrt(sp.pi) * sp.gamma(half_sum) / sp.gamma(
            half_sum + sp.Rational(1, 2)
        )
    return (
        2 ** (alpha + beta - 1)
        * sp.exp(-alpha * displacement)
        * sp.beta(half_sum, half_sum)
        * sp.hyper(
            (alpha, half_sum),
            (2 * half_sum,),
            1 - sp.exp(-2 * displacement),
        )
    )


def translated_sech_overlap(
    mode_power: Any,
    profile_power: Any,
    amplitude: Any,
    inverse_width: Any,
    displacement: Any,
) -> TranslatedSechOverlap:
    """Return a matched-width translated normalized density expectation.

    The normalized mode is proportional to
    ``sech(kappa*(x-R))**p`` and the real multiplier is
    ``A*sech(kappa*x)**r``.  This is a whole-line Cartesian result; it is not a
    half-line radial integral.
    """

    p = _positive_real(mode_power, "mode_power")
    r = _positive_real(profile_power, "profile_power")
    coefficient = _real(amplitude, "amplitude")
    kappa = _positive_real(inverse_width, "inverse_width")
    offset = _real(displacement, "displacement")
    scaled_offset = sp.simplify(kappa * sp.Abs(offset))
    density_power = sp.simplify(2 * p)
    density_norm = sp.sqrt(sp.pi) * sp.gamma(p) / sp.gamma(p + sp.Rational(1, 2))
    convolution = sech_convolution(density_power, r, scaled_offset)
    expectation = sp.simplify(coefficient * convolution / density_norm)
    return TranslatedSechOverlap(
        mode_power=p,
        profile_power=r,
        amplitude=coefficient,
        inverse_width=kappa,
        displacement=offset,
        dimensionless_displacement=scaled_offset,
        mode_density_norm=density_norm,
        raw_dimensionless_convolution=convolution,
        normalized_overlap=expectation,
    )


def sech_overlap_tail_ledger(
    mode_power: Any,
    profile_power: Any,
    amplitude: Any,
    inverse_width: Any,
) -> SechTailLedger:
    """Return the exact large-separation class of a matched-width overlap.

    With mode-density power ``alpha=2*p`` and profile power ``beta=r``, the
    dimensionless decay exponent is ``min(alpha,beta)``.  Unequal exponents
    have a constant leading prefactor.  Equal exponents acquire one factor of
    ``kappa*|R|``; omitting that resonance is a false pure-geometric claim.
    """

    p = _positive_real(mode_power, "mode_power")
    beta = _positive_real(profile_power, "profile_power")
    coefficient = _real(amplitude, "amplitude")
    kappa = _positive_real(inverse_width, "inverse_width")
    alpha = sp.simplify(2 * p)
    norm = sp.sqrt(sp.pi) * sp.gamma(p) / sp.gamma(p + sp.Rational(1, 2))

    if alpha.is_number and beta.is_number:
        difference = sp.simplify(alpha - beta)
        if difference == 0:
            decay_power = alpha
            polynomial_power = 1
            leading = sp.simplify(coefficient * 2 ** (2 * alpha) / norm)
        else:
            decay_power = sp.Min(alpha, beta)
            polynomial_power = 0
            leading = sp.simplify(
                coefficient
                * 2 ** (alpha + beta - 1)
                * sp.beta(sp.Abs(difference) / 2, (alpha + beta) / 2)
                / norm
            )
    else:
        decay_power = sp.Min(alpha, beta)
        polynomial_power = -1
        leading = sp.nan

    return SechTailLedger(
        mode_density_power=alpha,
        profile_power=beta,
        dimensionless_decay_power=decay_power,
        physical_decay_rate=sp.simplify(kappa * decay_power),
        polynomial_prefactor_power=polynomial_power,
        normalized_leading_coefficient=leading,
    )


def poschl_teller_ground_ledger(
    depth: Any,
    width: Any,
    center: Any = 0,
) -> PoschlTellerGroundLedger:
    """Return the exact ground ledger for ``-d2-V0*sech((x-R)/w)^2``.

    The index ``s`` is the positive solution of ``s*(s+1)=V0*w**2``.
    Translation changes the center but neither ``s`` nor the eigenvalue.
    """

    well_depth = _positive_real(depth, "depth")
    well_width = _positive_real(width, "width")
    origin = _real(center, "center")
    index = sp.simplify(
        (sp.sqrt(1 + 4 * well_depth * well_width**2) - 1) / 2
    )
    normalization = sp.sqrt(
        sp.gamma(index + sp.Rational(1, 2))
        / (well_width * sp.sqrt(sp.pi) * sp.gamma(index))
    )
    return PoschlTellerGroundLedger(
        depth=well_depth,
        width=well_width,
        center=origin,
        index=index,
        eigenvalue=sp.simplify(-index**2 / well_width**2),
        normalization=sp.simplify(normalization),
        density_tail_rate=sp.simplify(2 * index / well_width),
    )


def poschl_teller_ground_state(
    coordinate: Any,
    depth: Any,
    width: Any,
    center: Any = 0,
) -> sp.Expr:
    """Return the normalized exact translated Pöschl--Teller ground state."""

    x = sp.sympify(coordinate)
    ledger = poschl_teller_ground_ledger(depth, width, center)
    return sp.simplify(
        ledger.normalization
        * sp.sech((x - ledger.center) / ledger.width) ** ledger.index
    )


def poschl_teller_operator(
    mode: Any,
    coordinate: sp.Symbol,
    depth: Any,
    width: Any,
    center: Any = 0,
) -> sp.Expr:
    """Apply one declared translated Pöschl--Teller Hamiltonian."""

    function = sp.sympify(mode)
    well_depth = _positive_real(depth, "depth")
    well_width = _positive_real(width, "width")
    origin = _real(center, "center")
    potential = -well_depth * sp.sech((coordinate - origin) / well_width) ** 2
    return sp.simplify(-sp.diff(function, coordinate, 2) + potential * function)


def poschl_sech_overlap_tail_ledger(
    depth: Any,
    width: Any,
    profile_amplitude: Any,
    profile_inverse_width: Any,
    center: Any = 0,
) -> PoschlSechTailLedger:
    """Return the exact large-center tail ledger for a Pöschl ground density.

    The normalized density is the square of the exact ground state of
    ``-d2-V0*sech((x-R)/w)^2`` and the core multiplier is
    ``A*sech(kappa*x)``.  For unequal density-tail rate ``mu=2*s/w`` and core
    rate ``kappa``, the slower rate controls the overlap and the displayed beta
    coefficient is exact.  Equal rates acquire one factor of ``R``.
    """

    ground = poschl_teller_ground_ledger(depth, width, center)
    amplitude = _real(profile_amplitude, "profile_amplitude")
    kappa = _positive_real(profile_inverse_width, "profile_inverse_width")
    mode_rate = ground.density_tail_rate
    difference = sp.simplify(mode_rate - kappa)
    if difference.is_number is not True:
        return PoschlSechTailLedger(
            ground=ground,
            profile_amplitude=amplitude,
            profile_inverse_width=kappa,
            overlap_tail_rate=sp.Min(mode_rate, kappa),
            polynomial_prefactor_power=-1,
            leading_coefficient=sp.nan,
        )

    if difference == 0:
        leading = sp.simplify(
            amplitude * ground.normalization**2 * 2 ** (2 * ground.index + 1)
        )
        return PoschlSechTailLedger(
            ground=ground,
            profile_amplitude=amplitude,
            profile_inverse_width=kappa,
            overlap_tail_rate=kappa,
            polynomial_prefactor_power=1,
            leading_coefficient=leading,
        )

    if bool(difference > 0):
        q = sp.simplify(kappa * ground.width)
        leading = sp.simplify(
            amplitude
            * ground.normalization**2
            * ground.width
            * 2 ** (2 * ground.index)
            * sp.beta(ground.index - q / 2, ground.index + q / 2)
        )
        selected_rate = kappa
    else:
        q = sp.simplify(mode_rate / kappa)
        leading = sp.simplify(
            amplitude
            * ground.normalization**2
            * 2 ** (2 * ground.index)
            / kappa
            * sp.beta((1 + q) / 2, (1 - q) / 2)
        )
        selected_rate = mode_rate
    return PoschlSechTailLedger(
        ground=ground,
        profile_amplitude=amplitude,
        profile_inverse_width=kappa,
        overlap_tail_rate=selected_rate,
        polynomial_prefactor_power=0,
        leading_coefficient=leading,
    )


def tail_spacing_ledger(
    mode_density_tail_rate: Any,
    profile_tail_rate: Any,
    spacing: Any,
) -> TailSpacingLedger:
    """Return the conditional asymptotic log ratio for a linear center ladder.

    For exact exponential tails with positive rates ``mu`` and ``nu``, the
    overlap rate is the slower ``min(mu,nu)``.  A ladder ``R_n=R_0+n*d`` has
    limiting log ratio ``-min(mu,nu)*d``.  This product is not predicted while
    ``d`` remains free.
    """

    mode_rate = _positive_real(mode_density_tail_rate, "mode_density_tail_rate")
    profile_rate = _positive_real(profile_tail_rate, "profile_tail_rate")
    ladder_spacing = _positive_real(spacing, "spacing")
    overlap_rate = sp.Min(mode_rate, profile_rate)
    resonance: bool | None = None
    if mode_rate.is_number and profile_rate.is_number:
        resonance = bool(sp.simplify(mode_rate - profile_rate) == 0)
    return TailSpacingLedger(
        mode_density_tail_rate=mode_rate,
        profile_tail_rate=profile_rate,
        spacing=ladder_spacing,
        overlap_tail_rate=overlap_rate,
        asymptotic_log_ratio=sp.simplify(-overlap_rate * ladder_spacing),
        resonant_equal_rates=resonance,
    )


def reciprocal_rate_spacing_rescaling(
    tail_rate: Any,
    spacing: Any,
    factor: Any,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return ``(rho*mu,d/rho,mu*d)`` for a nonzero real ``rho``."""

    rate = _positive_real(tail_rate, "tail_rate")
    ladder_spacing = _positive_real(spacing, "spacing")
    rho = _positive_real(factor, "factor")
    return (
        sp.simplify(rho * rate),
        sp.simplify(ladder_spacing / rho),
        sp.simplify(rate * ladder_spacing),
    )


def normalized_gaussian_overlap(
    mode_density_rate: Any,
    profile_rate: Any,
    amplitude: Any,
    displacement: Any,
) -> sp.Expr:
    """Return an exact Gaussian localization countermodel.

    For normalized density ``sqrt(a/pi)*exp(-a*(x-R)^2)`` and profile
    ``A*exp(-b*x^2)``, the overlap is
    ``A*sqrt(a/(a+b))*exp(-a*b*R^2/(a+b))``.  A linear center ladder therefore
    has quadratic, not geometric, log attenuation.
    """

    density_rate = _positive_real(mode_density_rate, "mode_density_rate")
    multiplier_rate = _positive_real(profile_rate, "profile_rate")
    coefficient = _real(amplitude, "amplitude")
    offset = _real(displacement, "displacement")
    return sp.simplify(
        coefficient
        * sp.sqrt(density_rate / (density_rate + multiplier_rate))
        * sp.exp(
            -density_rate
            * multiplier_rate
            * offset**2
            / (density_rate + multiplier_rate)
        )
    )
