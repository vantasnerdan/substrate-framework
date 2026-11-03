"""Exact one-mode coherent-state coefficients and number probabilities.

The state in this module is a mathematical vector in the standard one-mode
Fock representation governed by :mod:`substrate_framework.bosonic_fock`.
Declaring a displacement amplitude does not prepare a material state, identify
a classical phase excursion, change the vacua of a cosine potential, or supply
an interaction, channel, transition probability, or rate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import sympy as sp

from .bosonic_fock import factorial_one_modes


CoherentOccupationSupport = Literal["vacuum_only", "all_nonnegative"]


def _nonnegative_integer(value: Any, *, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be an exact nonnegative integer")
    return int(expression)


def _exact_complex(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float) or expression.is_complex is not True:
        raise ValueError(f"{name} must be exact and explicitly complex-valued")
    return expression


def coherent_state_intensity(displacement: Any) -> sp.Expr:
    r"""Return the exact nonnegative intensity ``S=conjugate(alpha)*alpha``.

    The displacement must be an exact explicitly complex-valued expression.
    A floating-point approximation is rejected so that normalization, ties,
    and zero support are never decided by tolerance.
    """

    alpha = _exact_complex(displacement, name="displacement")
    intensity = sp.expand_complex(sp.conjugate(alpha) * alpha)
    intensity = sp.simplify(intensity)
    if intensity.is_nonnegative is not True:
        raise ValueError("displacement intensity must be provably nonnegative")
    return intensity


def coherent_state_coefficient(order: Any, *, displacement: Any) -> sp.Expr:
    r"""Return ``<n|alpha>=exp(-S/2)*alpha**n/sqrt(n!)`` exactly.

    These coefficients define the normalized norm-convergent coherent vector

    ``|alpha> = exp(-S/2) sum_n alpha**n/sqrt(n!) |n>``.

    In the standard Weyl representation this vector is
    ``D(alpha)|0>`` for
    ``D(alpha)=exp(alpha*a_dagger-conjugate(alpha)*a)``.  The Gaussian half
    factor is load bearing: omitting it makes the squared norm ``exp(S)``.
    """

    n = _nonnegative_integer(order, name="order")
    alpha = _exact_complex(displacement, name="displacement")
    intensity = coherent_state_intensity(alpha)
    return sp.simplify(
        sp.exp(-intensity / 2) * alpha**n / sp.sqrt(sp.factorial(n))
    )


def coherent_state_number_probability(
    order: Any,
    *,
    displacement: Any,
) -> sp.Expr:
    r"""Return the exact number-measurement probability ``exp(-S)S**n/n!``.

    This is a genuine Born probability for a number-basis measurement on the
    separately declared coherent state.  It is not a probability that a
    material event, transition, reaction, or branching channel occurs.
    """

    n = _nonnegative_integer(order, name="order")
    intensity = coherent_state_intensity(displacement)
    return sp.factor(sp.exp(-intensity) * intensity**n / sp.factorial(n))


def coherent_state_overlap(alpha: Any, beta: Any) -> sp.Expr:
    r"""Return the exact overlap ``<alpha|beta>``.

    With normalized coherent vectors the result is
    ``exp(-(S_alpha+S_beta)/2+conjugate(alpha)*beta)``.  Its diagonal is one;
    nonorthogonality is a Hilbert-space property, not a physical coherence or
    transition-rate claim.
    """

    first = _exact_complex(alpha, name="alpha")
    second = _exact_complex(beta, name="beta")
    first_intensity = coherent_state_intensity(first)
    second_intensity = coherent_state_intensity(second)
    return sp.simplify(
        sp.exp(
            -(first_intensity + second_intensity) / 2
            + sp.conjugate(first) * second
        )
    )


def coherent_state_number_modes(*, displacement: Any) -> tuple[int, ...]:
    """Return every exact number-basis mode for rational intensity.

    Zero displacement has the vacuum as its sole mode.  Positive rational
    intensity reuses the complete C-CMB-003 mode convention, including the
    adjacent tie at every positive integer intensity.  Irrational or
    undecidable symbolic intensity is rejected by this constructive API.
    """

    intensity = coherent_state_intensity(displacement)
    if intensity == 0:
        return (0,)
    if intensity.is_Rational is not True or intensity.is_positive is not True:
        raise ValueError("displacement intensity must be an exact positive rational")
    return factorial_one_modes(
        intensity=sp.Rational(intensity),
        support="all_nonnegative",
    )


@dataclass(frozen=True)
class CoherentStateLedger:
    """Exact state, support, moment, and interpretation data."""

    displacement: sp.Expr
    intensity: sp.Expr
    vacuum_amplitude: sp.Expr
    vacuum_probability: sp.Expr
    mean_occupation: sp.Expr
    occupation_variance: sp.Expr
    annihilation_eigenvalue: sp.Expr
    occupation_support: CoherentOccupationSupport
    normalized: bool
    physical_preparation_is_separate_premise: bool
    classical_phase_map_is_separate_premise: bool
    physical_process_probability_is_separate_premise: bool


def coherent_state_ledger(*, displacement: Any) -> CoherentStateLedger:
    r"""Return exact coherent-state invariants for one declared ``alpha``.

    The coefficient recurrence
    ``sqrt(n+1)*<n+1|alpha> = alpha*<n|alpha>`` gives
    ``a|alpha>=alpha|alpha>``.  Number measurement has mean and variance
    ``S=|alpha|**2``.  Every nonnegative occupation has positive probability
    when ``alpha`` is nonzero, but this unbounded discrete support does not
    imply an unbounded classical excursion or additional potential vacua.
    """

    alpha = _exact_complex(displacement, name="displacement")
    intensity = coherent_state_intensity(alpha)
    support: CoherentOccupationSupport = (
        "vacuum_only" if intensity == 0 else "all_nonnegative"
    )
    return CoherentStateLedger(
        displacement=alpha,
        intensity=intensity,
        vacuum_amplitude=sp.exp(-intensity / 2),
        vacuum_probability=sp.exp(-intensity),
        mean_occupation=intensity,
        occupation_variance=intensity,
        annihilation_eigenvalue=alpha,
        occupation_support=support,
        normalized=True,
        physical_preparation_is_separate_premise=True,
        classical_phase_map_is_separate_premise=True,
        physical_process_probability_is_separate_premise=True,
    )
