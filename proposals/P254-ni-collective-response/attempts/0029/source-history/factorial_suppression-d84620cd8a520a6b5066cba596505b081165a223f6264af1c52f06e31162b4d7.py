"""Exact suppression ledgers for inverse-square factorial sequences.

The helpers compose the classical coefficient theorem in
:mod:`substrate_framework.cosine_vertices` with elementary exact factorial
bounds.  A coefficient square is not by itself a matrix element, probability,
transition rate, branching weight, or physical subdivision law.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Any

import sympy as sp

from .cosine_vertices import vacuum_one_high_coefficient


def _integer_at_least(value: Any, *, name: str, minimum: int) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.has(sp.Float)
    ):
        raise ValueError(f"{name} must be an exact integer")
    integer = int(expression)
    if integer < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return integer


def _exact_real(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float) or expression.is_real is not True:
        raise ValueError(f"{name} must be exact and provably real")
    return expression


def _exact_positive(value: Any, *, name: str) -> sp.Expr:
    expression = _exact_real(value, name=name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


def _exact_positive_rational(value: Any, *, name: str) -> sp.Rational:
    expression = sp.sympify(value)
    if (
        expression.has(sp.Float)
        or expression.is_Rational is not True
        or expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be an exact positive rational")
    return sp.Rational(expression)


def exact_rational_log10_floor(value: Any) -> int:
    """Return ``floor(log10(value))`` using exact integer comparisons.

    ``value`` must be a positive exact rational.  No logarithm or floating-
    point conversion is evaluated.
    """

    expression = sp.sympify(value)
    if expression.has(sp.Float) or expression.is_Rational is not True:
        raise ValueError("value must be an exact rational")
    rational = sp.Rational(expression)
    if rational <= 0:
        raise ValueError("value must be positive")
    numerator = int(sp.numer(rational))
    denominator = int(sp.denom(rational))

    def at_least_power_of_ten(exponent: int) -> bool:
        if exponent >= 0:
            return numerator >= denominator * 10**exponent
        return numerator * 10 ** (-exponent) >= denominator

    exponent = len(str(numerator)) - len(str(denominator))
    while not at_least_power_of_ten(exponent):
        exponent -= 1
    while at_least_power_of_ten(exponent + 1):
        exponent += 1
    return exponent


def cosine_one_high_coefficient_square(
    low_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
) -> sp.Expr:
    """Return the exact square of the declared real one-high coefficient.

    At zero background this is zero for even ``low_order``.  For odd ``n`` it
    is ``A**2*a_H**2*a_L**(2*n)/(n!)**2``.  The explicit normalization factors
    prevent the unit specialization from being mistaken for a universal
    transition weight.
    """

    order = _integer_at_least(low_order, name="low_order", minimum=0)
    amplitude_expression = _exact_real(amplitude, name="amplitude")
    high_expression = _exact_real(high_scale, name="high_scale")
    low_expression = _exact_real(low_scale, name="low_scale")
    coefficient = vacuum_one_high_coefficient(
        order,
        amplitude=amplitude_expression,
        high_scale=high_expression,
        low_scale=low_expression,
    )
    return sp.factor(coefficient**2)


@dataclass(frozen=True)
class FactorialSuppressionEvidence:
    """Exact data for ``q_n=1/(n!)**2`` at one positive integer order."""

    order: int
    inverse_square_factorial: sp.Rational
    next_inverse_square_factorial: sp.Rational
    recurrence_ratio: sp.Rational
    exact_log10_floor: int
    exponential_upper_bound: sp.Expr
    physical_rate_interpretation_is_separate_premise: bool


def factorial_suppression_evidence(order: Any) -> FactorialSuppressionEvidence:
    """Return exact sequence, recurrence, decimal, and exponential-bound data."""

    n = _integer_at_least(order, name="order", minimum=1)
    value = sp.Rational(1, int(sp.factorial(n)) ** 2)
    next_value = sp.Rational(1, int(sp.factorial(n + 1)) ** 2)
    return FactorialSuppressionEvidence(
        order=n,
        inverse_square_factorial=value,
        next_inverse_square_factorial=next_value,
        recurrence_ratio=sp.simplify(next_value / value),
        exact_log10_floor=exact_rational_log10_floor(value),
        exponential_upper_bound=sp.Pow(sp.E / n, 2 * n),
        physical_rate_interpretation_is_separate_premise=True,
    )


@dataclass(frozen=True)
class SuperpolynomialTailEvidence:
    """A geometric tail certificate for ``n**p/(n!)**2``."""

    power: int
    start_order: int
    ratio_ceiling: sp.Rational
    exact_ratio_at_start: sp.Rational


def factorial_superpolynomial_tail(power: Any) -> SuperpolynomialTailEvidence:
    """Return an exact tail where consecutive terms have ratio at most one-half.

    For fixed nonnegative integer ``p``, put ``a_n=n**p/(n!)**2``.  Since
    ``a_(n+1)/a_n <= 2**p/(n+1)**2`` for ``n>=1``, choosing
    ``n+1>=ceil(sqrt(2**(p+1)))`` gives a geometric tail with ratio at most
    one-half and hence ``a_n -> 0``.
    """

    p = _integer_at_least(power, name="power", minimum=0)
    target = 2 ** (p + 1)
    root = isqrt(target)
    ceiling_root = root if root * root == target else root + 1
    start = max(1, ceiling_root - 1)
    ratio = sp.Rational((start + 1) ** p, start**p * (start + 1) ** 2)
    return SuperpolynomialTailEvidence(
        power=p,
        start_order=start,
        ratio_ceiling=sp.Rational(1, 2),
        exact_ratio_at_start=ratio,
    )


@dataclass(frozen=True)
class FactorialDecadeBound:
    """An exact base-ten exponent ceiling at ``n=10**decade``."""

    decade: int
    order: int
    e_series_upper_bound: sp.Rational
    convenient_e_upper_bound: sp.Rational
    twentieth_power_left: int
    twentieth_power_right: int
    twentieth_power_margin: int
    positive_exponent: int
    log10_upper_bound: int
    exposed_order_is_separate_physical_premise: bool


def factorial_decade_bound(decade: Any) -> FactorialDecadeBound:
    r"""Return the exact conservative bound at ``n=10**decade``.

    The exponential series gives ``e<49/18<11/4`` and exact integer arithmetic
    gives ``(11/4)**20<10**9``.  Grouping ``2*n=20*(n/10)`` therefore proves

    ``1/(n!)**2 < 10**(-((20*decade-9)*n/10))``

    for every positive integer decade.  The returned exponent is an integer;
    the enormous decimal itself is never materialized.
    """

    d = _integer_at_least(decade, name="decade", minimum=1)
    n = 10**d
    left = 11**20
    right = 10**9 * 4**20
    positive_exponent = (20 * d - 9) * n // 10
    return FactorialDecadeBound(
        decade=d,
        order=n,
        e_series_upper_bound=sp.Rational(49, 18),
        convenient_e_upper_bound=sp.Rational(11, 4),
        twentieth_power_left=left,
        twentieth_power_right=right,
        twentieth_power_margin=right - left,
        positive_exponent=positive_exponent,
        log10_upper_bound=-positive_exponent,
        exposed_order_is_separate_physical_premise=True,
    )


def parity_thinned_factorial_mass(
    order: Any,
    *,
    activity: Any = 1,
) -> sp.Expr:
    r"""Return ``1_odd(n) * z**(2*n)/(n!)**2`` exactly.

    The ambient sample space is the positive integers with counting measure;
    the mass has support on the positive odd integers. ``activity`` is a
    declared positive dimensionless normalization parameter. The returned
    mathematical mass does not by itself define an occurrence law, matrix
    element, transition rate, or physical channel.
    """

    n = _integer_at_least(order, name="order", minimum=1)
    z = _exact_positive(activity, name="activity")
    if n % 2 == 0:
        return sp.Integer(0)
    return sp.factor(z ** (2 * n) / sp.factorial(n) ** 2)


def odd_factorial_total_mass(activity: Any = 1) -> sp.Expr:
    r"""Return the exact total positive-odd factorial mass.

    Parity filtering the defining power series gives

    ``sum_{n positive odd} z**(2*n)/(n!)**2``
    ``= (besseli(0, 2*z) - besselj(0, 2*z))/2``.
    """

    z = _exact_positive(activity, name="activity")
    return (sp.besseli(0, 2 * z) - sp.besselj(0, 2 * z)) / 2


def normalized_parity_factorial_mass(
    order: Any,
    *,
    activity: Any = 1,
) -> sp.Expr:
    """Return the exact normalized mathematical mass at one integer order."""

    mass = parity_thinned_factorial_mass(order, activity=activity)
    if mass == 0:
        return sp.Integer(0)
    return sp.factor(mass / odd_factorial_total_mass(activity))


@dataclass(frozen=True)
class OddFactorialMassEnclosure:
    """Exact rational enclosure for an odd factorial-mass normalization."""

    activity: sp.Rational
    maximum_odd_order: int
    partial_mass: sp.Rational
    first_omitted_order: int
    first_omitted_mass: sp.Rational
    tail_ratio_ceiling: sp.Rational
    tail_mass_upper_bound: sp.Rational
    total_mass_lower_bound: sp.Rational
    total_mass_upper_bound: sp.Rational
    normalized_tail_upper_bound: sp.Rational


def odd_factorial_mass_enclosure(
    maximum_odd_order: Any,
    *,
    activity: Any = 1,
) -> OddFactorialMassEnclosure:
    r"""Enclose the infinite mass by an exact geometric tail.

    For positive odd ``M``, the first omitted order is ``M+2``. Consecutive
    omitted odd terms have ratios bounded by
    ``z**4/((M+3)**2*(M+4)**2)``. The caller must choose ``M`` so this bound is
    below one. Rational ``activity`` keeps the complete certificate exact.
    """

    maximum = _integer_at_least(
        maximum_odd_order,
        name="maximum_odd_order",
        minimum=1,
    )
    if maximum % 2 == 0:
        raise ValueError("maximum_odd_order must be odd")
    z = _exact_positive_rational(activity, name="activity")
    partial = sp.Rational(
        sum(
            parity_thinned_factorial_mass(order, activity=z)
            for order in range(1, maximum + 1, 2)
        )
    )
    first_order = maximum + 2
    first = sp.Rational(
        parity_thinned_factorial_mass(first_order, activity=z)
    )
    ratio = sp.Rational(
        z**4 / ((maximum + 3) ** 2 * (maximum + 4) ** 2)
    )
    if ratio >= 1:
        raise ValueError("maximum_odd_order must make the tail ratio below one")
    tail_upper = sp.factor(first / (1 - ratio))
    total_upper = sp.factor(partial + tail_upper)
    return OddFactorialMassEnclosure(
        activity=z,
        maximum_odd_order=maximum,
        partial_mass=partial,
        first_omitted_order=first_order,
        first_omitted_mass=first,
        tail_ratio_ceiling=ratio,
        tail_mass_upper_bound=tail_upper,
        total_mass_lower_bound=partial,
        total_mass_upper_bound=total_upper,
        normalized_tail_upper_bound=sp.factor(tail_upper / partial),
    )
