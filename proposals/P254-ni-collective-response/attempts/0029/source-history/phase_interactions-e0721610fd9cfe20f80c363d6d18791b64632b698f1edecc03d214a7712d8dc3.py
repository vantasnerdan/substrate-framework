"""Exact trial-profile phase interactions and scalar-circle packing.

The pair interaction in this module belongs to one separately declared
quartic energy evaluated on a superposition of two fixed profiles.  It is not
a common nonlinear two-soliton solution, a force law outside the displayed
collective coordinate, or a stability theorem.  The packing result concerns a
complete graph of scalar phases on one circle; changing the graph or internal
orientation space changes the capacity problem.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Any, Iterable

import sympy as sp


def _positive_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and (
        result.is_real is not True or result.is_positive is not True
    ):
        raise ValueError(f"{name} must be positive and real")
    return result


def _phase_cosine(value: Any) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number:
        if result.is_real is not True or not bool(-1 <= result <= 1):
            raise ValueError("phase_cosine must be real and lie in [-1, 1]")
    return result


def sech_pair_mixed_cubic_shape(scaled_separation: Any) -> sp.Expr:
    """Return the dimensionless ``sech^3*sech`` overlap shape.

    For ``s>0`` this is

    ``J31(s)=2*(sinh(s)*cosh(s)-s)/sinh(s)^3``.

    With equal profiles centered a distance ``d`` apart and
    ``s=kappa*d``, the dimensional overlap is ``A^4*J31(s)/kappa``.
    """

    s = _positive_real(scaled_separation, "scaled_separation")
    return sp.factor(
        2 * (sp.sinh(s) * sp.cosh(s) - s) / sp.sinh(s) ** 3
    )


def sech_pair_density_shape(scaled_separation: Any) -> sp.Expr:
    """Return the dimensionless ``sech^2*sech^2`` overlap shape.

    For ``s>0`` this is

    ``J22(s)=4*(s*cosh(s)-sinh(s))/sinh(s)^3``.
    """

    s = _positive_real(scaled_separation, "scaled_separation")
    return sp.factor(
        4 * (s * sp.cosh(s) - sp.sinh(s)) / sp.sinh(s) ** 3
    )


@dataclass(frozen=True)
class QuarticSechPairInteraction:
    """Exact cross-energy ledger for two fixed equal-width sech profiles."""

    separation: sp.Expr
    scaled_separation: sp.Expr
    phase_cosine: sp.Expr
    amplitude: sp.Expr
    inverse_width: sp.Expr
    mixed_cubic_overlap: sp.Expr
    density_overlap: sp.Expr
    linear_phase_energy: sp.Expr
    nonlinear_phase_energy: sp.Expr
    interaction_energy: sp.Expr


def quartic_sech_pair_interaction(
    separation: Any,
    phase_cosine: Any,
    amplitude: Any,
    inverse_width: Any,
) -> QuarticSechPairInteraction:
    """Return the exact two-profile cross energy of the declared functional.

    Let ``f_1=A*sech(kappa*(x+d/2))`` and
    ``f_2=A*sech(kappa*(x-d/2))`` with ``d,kappa,A>0``, and set
    ``Phi=f_1+exp(i*delta)*f_2`` and ``c=cos(delta)``.  For

    ``E[Phi]=integral (|Phi'|^2+kappa^2|Phi|^2-|Phi|^4/24) dx``,

    subtraction of both isolated energies and use of the profile equation
    gives

    ``E_int=-(c/6)*I31-((1+2*c^2)/12)*I22``.

    Thus only the leading large-separation term is linear in ``c``.  The
    phase-independent and ``c^2`` terms are retained rather than hidden in a
    fitted remainder.
    """

    d = _positive_real(separation, "separation")
    c = _phase_cosine(phase_cosine)
    coefficient = _positive_real(amplitude, "amplitude")
    kappa = _positive_real(inverse_width, "inverse_width")
    scaled = sp.simplify(kappa * d)
    common_scale = coefficient**4 / kappa
    cubic = sp.factor(common_scale * sech_pair_mixed_cubic_shape(scaled))
    density = sp.factor(common_scale * sech_pair_density_shape(scaled))
    linear = sp.factor(-c * cubic / 6)
    nonlinear = sp.factor(-(1 + 2 * c**2) * density / 12)
    return QuarticSechPairInteraction(
        separation=d,
        scaled_separation=scaled,
        phase_cosine=c,
        amplitude=coefficient,
        inverse_width=kappa,
        mixed_cubic_overlap=cubic,
        density_overlap=density,
        linear_phase_energy=linear,
        nonlinear_phase_energy=nonlinear,
        interaction_energy=sp.factor(linear + nonlinear),
    )


def pairwise_phase_cosines(phases: Iterable[Any]) -> tuple[sp.Expr, ...]:
    """Return all complete-graph pairwise scalar-phase cosines exactly."""

    values = tuple(sp.sympify(phase) for phase in phases)
    if len(values) < 2:
        raise ValueError("at least two phases are required")
    if any(value.is_number and value.is_real is not True for value in values):
        raise ValueError("phases must be real")
    return tuple(
        sp.simplify(sp.cos(values[i] - values[j]))
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )


@dataclass(frozen=True)
class ScalarCirclePacking:
    """Sharp complete-graph scalar-phase packing data at one count."""

    count: int
    nearest_gap_upper_bound: sp.Expr
    optimal_worst_pairwise_cosine: sp.Expr
    regular_phases: tuple[sp.Expr, ...]
    regular_pairwise_cosines: tuple[sp.Expr, ...]
    strictly_negative_possible: bool
    nonpositive_possible: bool


@dataclass(frozen=True)
class CompletePhaseCosineLedger:
    """Exact complete-graph cosine sum and resultant representation."""

    count: int
    phases: tuple[sp.Expr, ...]
    phase_resultant: sp.Expr
    resultant_squared: sp.Expr
    pairwise_cosine_sum: sp.Expr
    minimum: sp.Expr
    excess_above_minimum: sp.Expr


def complete_phase_cosine_ledger(
    phases: Iterable[Any],
) -> CompletePhaseCosineLedger:
    """Return the exact identity for a complete scalar-phase cosine sum.

    For ``N`` unit phasors,

    ``sum_{a<b} cos(theta_a-theta_b)=(|sum_a exp(i*theta_a)|^2-N)/2``.

    The minimum is therefore ``-N/2`` and is attained exactly when the
    phasor resultant vanishes.  Regular polygons attain it for every
    ``N>=2``, but the minimizing set need not be unique.  This static surrogate
    is not a physical energy unless a separate interaction claim identifies
    it as one.
    """

    values = tuple(sp.sympify(phase) for phase in phases)
    if len(values) < 2:
        raise ValueError("at least two phases are required")
    if any(value.is_number and value.is_real is not True for value in values):
        raise ValueError("phases must be real")
    resultant = sp.trigsimp(
        sp.expand_complex(sum(sp.exp(sp.I * value) for value in values))
    )
    resultant_squared = sp.trigsimp(
        sp.expand_complex(resultant * sp.conjugate(resultant))
    )
    direct_sum = sp.trigsimp(
        sum(
            sp.cos(values[i] - values[j])
            for i in range(len(values))
            for j in range(i + 1, len(values))
        )
    )
    minimum = -sp.Rational(len(values), 2)
    excess = sp.simplify(resultant_squared / 2)
    if sp.trigsimp(direct_sum - minimum - excess) != 0:
        raise ValueError("failed to establish the complete cosine identity")
    return CompletePhaseCosineLedger(
        count=len(values),
        phases=values,
        phase_resultant=resultant,
        resultant_squared=sp.simplify(resultant_squared),
        pairwise_cosine_sum=sp.simplify(direct_sum),
        minimum=minimum,
        excess_above_minimum=excess,
    )


def scalar_circle_packing(count: int) -> ScalarCirclePacking:
    """Return the sharp worst-pair cosine for ``count`` scalar phases.

    Sorting ``N`` phases gives a nearest circular gap at most ``2*pi/N``.
    Since cosine decreases on ``[0,pi]``, some pair has cosine at least
    ``cos(2*pi/N)``.  The regular ``N``-gon attains that bound, so the optimum
    possible value of the largest pairwise cosine is exactly
    ``cos(2*pi/N)``.  Consequently strict negativity is possible exactly for
    ``N<=3`` and nonpositivity exactly for ``N<=4``.  These are capacities,
    not occupancy or stability selection rules.
    """

    if isinstance(count, bool) or not isinstance(count, Integral):
        raise TypeError("count must be an integer")
    n = int(count)
    if n < 2:
        raise ValueError("count must be at least two")
    gap = sp.Rational(2, n) * sp.pi
    phases = tuple(sp.Rational(2 * index, n) * sp.pi for index in range(n))
    cosines = pairwise_phase_cosines(phases)
    return ScalarCirclePacking(
        count=n,
        nearest_gap_upper_bound=gap,
        optimal_worst_pairwise_cosine=sp.simplify(sp.cos(gap)),
        regular_phases=phases,
        regular_pairwise_cosines=cosines,
        strictly_negative_possible=n <= 3,
        nonpositive_possible=n <= 4,
    )
