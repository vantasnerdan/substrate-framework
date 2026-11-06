"""Exact allocation between two declared nonnegative rate-valued inputs.

This module normalizes two inputs that a caller has already established to
share one rate dimension. It does not derive physical states, interactions,
decay channels, final-state measures, kinetics, or material parameters merely
because the inputs are named rates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import sympy as sp


def _exact_real(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


def _exact_nonnegative(value: Any, *, name: str) -> sp.Expr:
    expression = _exact_real(value, name=name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be explicitly nonnegative")
    return expression


def _exact_positive(value: Any, *, name: str) -> sp.Expr:
    expression = _exact_real(value, name=name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


def _positive_integer(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_integer is not True or expression.is_positive is not True:
        raise ValueError(f"{name} must be a positive integer")
    return expression


@dataclass(frozen=True)
class TwoChannelAllocation:
    """Exact normalized shares of two declared common-dimension inputs."""

    first_rate: sp.Expr
    second_rate: sp.Expr
    total_rate: sp.Expr
    first_fraction: sp.Expr
    second_fraction: sp.Expr


@dataclass(frozen=True)
class WeightedChannelAllocation:
    """Exact specialization with one dimensionless weight and population."""

    baseline_weighted_rate: sp.Expr
    comparison_rate: sp.Expr
    weight: sp.Expr
    population: sp.Expr
    weighted_rate: sp.Expr
    baseline_ratio: sp.Expr
    weighted_fraction: sp.Expr
    comparison_fraction: sp.Expr


WeightedComparisonMonotonicity = Literal["decreasing", "stationary", "increasing"]


@dataclass(frozen=True)
class PopulationDependentWeightLedger:
    """Exact local derivative data for a population-dependent weight."""

    population: sp.Expr
    comparison_ratio: sp.Expr
    weight: sp.Expr
    weight_derivative: sp.Expr
    comparison_fraction: sp.Expr
    derivative_control: sp.Expr
    comparison_fraction_derivative: sp.Expr
    monotonicity: WeightedComparisonMonotonicity
    physical_weight_law_is_separate_premise: bool
    exhaustive_channel_interpretation_is_separate_premise: bool


def two_channel_allocation(
    first_rate: Any,
    second_rate: Any,
) -> TwoChannelAllocation:
    """Normalize two exact nonnegative common-dimension inputs.

    At least one input must be explicitly positive, which excludes the
    undefined double-zero denominator while retaining either zero endpoint.
    The caller is responsible for establishing that both inputs share one rate
    dimension and correspond to an exhaustive physical channel set if that
    interpretation is desired.
    """

    first = _exact_nonnegative(first_rate, name="first_rate")
    second = _exact_nonnegative(second_rate, name="second_rate")
    total = sp.simplify(first + second)
    if total.is_positive is not True:
        raise ValueError("first_rate + second_rate must be explicitly positive")
    return TwoChannelAllocation(
        first_rate=first,
        second_rate=second,
        total_rate=total,
        first_fraction=sp.simplify(first / total),
        second_fraction=sp.simplify(second / total),
    )


def channel_odds(numerator_rate: Any, denominator_rate: Any) -> sp.Expr:
    """Return an exact channel-rate ratio with a positive denominator."""

    numerator = _exact_nonnegative(numerator_rate, name="numerator_rate")
    denominator = _exact_positive(denominator_rate, name="denominator_rate")
    return sp.simplify(numerator / denominator)


def weighted_channel_allocation(
    baseline_weighted_rate: Any,
    comparison_rate: Any,
    weight: Any,
    population: Any,
) -> WeightedChannelAllocation:
    """Specialize to ``r_w * weight * population`` and ``r_c``.

    Both baseline inputs are exact and positive. ``weight`` is an exact
    positive dimensionless factor and ``population`` an exact positive integer.
    These typing premises do not derive a physical rate law.
    """

    baseline = _exact_positive(
        baseline_weighted_rate,
        name="baseline_weighted_rate",
    )
    comparison = _exact_positive(comparison_rate, name="comparison_rate")
    dimensionless_weight = _exact_positive(weight, name="weight")
    count = _positive_integer(population, name="population")
    weighted = sp.simplify(baseline * dimensionless_weight * count)
    allocation = two_channel_allocation(weighted, comparison)
    ratio = sp.simplify(comparison / baseline)
    return WeightedChannelAllocation(
        baseline_weighted_rate=baseline,
        comparison_rate=comparison,
        weight=dimensionless_weight,
        population=count,
        weighted_rate=weighted,
        baseline_ratio=ratio,
        weighted_fraction=allocation.first_fraction,
        comparison_fraction=allocation.second_fraction,
    )


def relative_weighted_odds_enhancement(
    weight: Any,
    population: Any,
    baseline_weight: Any,
) -> sp.Expr:
    """Return ``weight*population/baseline_weight`` exactly.

    This is the ratio of the weighted-to-comparison odds at the declared point
    to the odds with population one and ``baseline_weight``. Cancellation also
    assumes the same positive channel normalizations in numerator and
    denominator. The result does not determine ``weight`` or a physical
    enhancement magnitude.
    """

    current = _exact_positive(weight, name="weight")
    count = _positive_integer(population, name="population")
    baseline = _exact_positive(baseline_weight, name="baseline_weight")
    return sp.simplify(current * count / baseline)


def population_dependent_weight_ledger(
    population: Any,
    comparison_ratio: Any,
    weight: Any,
    weight_derivative: Any,
) -> PopulationDependentWeightLedger:
    r"""Return the exact local population derivative for ``w=w(N)``.

    For ``B_c(N)=rho/(N*w(N)+rho)``, the chain rule gives

    ``B_c'(N)=-rho*(w(N)+N*w'(N))/(N*w(N)+rho)**2``.

    Thus a positive weight alone does not determine monotonicity. The
    comparison fraction decreases, is stationary, or increases according as
    ``w+N*w'`` is positive, zero, or negative. The caller separately supplies
    the local weight and derivative and remains responsible for deriving a
    physical weight law, common rate dimensions, and exhaustive channels.
    """

    count = _exact_positive(population, name="population")
    rho = _exact_positive(comparison_ratio, name="comparison_ratio")
    current_weight = _exact_positive(weight, name="weight")
    slope = _exact_real(weight_derivative, name="weight_derivative")
    control = sp.simplify(current_weight + count * slope)
    if control.is_positive is True:
        monotonicity: WeightedComparisonMonotonicity = "decreasing"
    elif control.is_zero is True:
        monotonicity = "stationary"
    elif control.is_negative is True:
        monotonicity = "increasing"
    else:
        raise ValueError("w + N*w' must have an explicitly decidable sign")
    denominator = sp.simplify(count * current_weight + rho)
    return PopulationDependentWeightLedger(
        population=count,
        comparison_ratio=rho,
        weight=current_weight,
        weight_derivative=slope,
        comparison_fraction=sp.simplify(rho / denominator),
        derivative_control=control,
        comparison_fraction_derivative=sp.simplify(
            -rho * control / denominator**2
        ),
        monotonicity=monotonicity,
        physical_weight_law_is_separate_premise=True,
        exhaustive_channel_interpretation_is_separate_premise=True,
    )
