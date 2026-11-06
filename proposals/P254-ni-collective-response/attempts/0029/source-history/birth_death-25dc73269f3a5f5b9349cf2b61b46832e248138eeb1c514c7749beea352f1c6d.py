"""Exact immigration--death dynamics with a factorial-one stationary law.

The process in this module is a separately declared continuous-time Markov
chain on the nonnegative integers.  Its stationary mass agrees algebraically
with the all-nonnegative factorial-one law, but that static law does not select
this generator, its rate scale, its initial law, or a material realization.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .bosonic_fock import normalized_factorial_one_mass


def _nonnegative_integer(value: Any, *, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.has(sp.Float)
        or expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be an exact nonnegative integer")
    return int(expression)


def _exact_real(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float) or expression.is_real is not True:
        raise ValueError(f"{name} must be exact and explicitly real")
    return expression


def _exact_positive(value: Any, *, name: str) -> sp.Expr:
    expression = _exact_real(value, name=name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


def _exact_nonnegative(value: Any, *, name: str) -> sp.Expr:
    expression = _exact_real(value, name=name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be explicitly nonnegative")
    return expression


def immigration_death_rates(
    state: Any,
    *,
    stationary_mean: Any,
    rate: Any,
) -> tuple[sp.Expr, sp.Expr]:
    r"""Return ``(lambda_n, mu_n)=(r*S, r*n)`` exactly.

    The boundary death rate is therefore ``mu_0=0``.  ``rate`` has inverse
    time units and ``stationary_mean`` is dimensionless; both are independent
    model inputs rather than consequences of a static probability mass.
    """

    n = _nonnegative_integer(state, name="state")
    S = _exact_positive(stationary_mean, name="stationary_mean")
    r = _exact_positive(rate, name="rate")
    return sp.simplify(r * S), sp.simplify(r * n)


def immigration_death_generator_action(
    state: Any,
    *,
    lower_value: Any,
    current_value: Any,
    upper_value: Any,
    stationary_mean: Any,
    rate: Any,
) -> sp.Expr:
    r"""Return ``L f(n)`` for three declared adjacent function values.

    At the reflecting state-space boundary ``n=0`` the lower value is ignored
    because ``mu_0=0``.  In particular, applying this generator to
    ``f(n)=n`` gives the local conditional drift ``r*(S-n)``.  A positive
    local drift is not monotone sample-path growth: deaths still occur at
    every positive state.
    """

    n = _nonnegative_integer(state, name="state")
    lower = sp.sympify(lower_value)
    current = sp.sympify(current_value)
    upper = sp.sympify(upper_value)
    if any(value.has(sp.Float) for value in (lower, current, upper)):
        raise ValueError("generator values must be exact")
    birth, death = immigration_death_rates(
        n,
        stationary_mean=stationary_mean,
        rate=rate,
    )
    return sp.simplify(birth * (upper - current) + death * (lower - current))


def immigration_death_local_drift(
    state: Any,
    *,
    stationary_mean: Any,
    rate: Any,
) -> sp.Expr:
    """Return the exact identity-function drift ``r*(S-n)``."""

    n = _nonnegative_integer(state, name="state")
    return immigration_death_generator_action(
        n,
        lower_value=n - 1,
        current_value=n,
        upper_value=n + 1,
        stationary_mean=stationary_mean,
        rate=rate,
    )


def immigration_death_stationary_mass(
    state: Any,
    *,
    stationary_mean: Any,
) -> sp.Expr:
    """Return ``pi_n=exp(-S)*S**n/n!`` on all nonnegative states."""

    return normalized_factorial_one_mass(
        state,
        intensity=_exact_positive(stationary_mean, name="stationary_mean"),
        support="all_nonnegative",
    )


def immigration_death_mean(
    time: Any,
    *,
    initial_mean: Any,
    stationary_mean: Any,
    rate: Any,
) -> sp.Expr:
    r"""Return ``S+(m0-S)*exp(-r*t)`` for a declared finite initial mean."""

    t = _exact_nonnegative(time, name="time")
    m0 = _exact_nonnegative(initial_mean, name="initial_mean")
    S = _exact_positive(stationary_mean, name="stationary_mean")
    r = _exact_positive(rate, name="rate")
    return sp.simplify(S + (m0 - S) * sp.exp(-r * t))


def immigration_death_probability_generating_function(
    variable: sp.Symbol,
    time: Any,
    *,
    initial_generating_function: Any,
    stationary_mean: Any,
    rate: Any,
) -> sp.Expr:
    r"""Return the exact time-dependent probability generating function.

    For initial PGF ``G0(z)``,

    ``G(z,t)=G0(1+(z-1)e**(-r*t))*exp(S*(z-1)*(1-e**(-r*t)))``.

    The caller remains responsible for supplying a normalized initial PGF.
    This function transforms it; it does not infer an initial physical state.
    """

    if not isinstance(variable, sp.Symbol):
        raise ValueError("variable must be a SymPy Symbol")
    t = _exact_nonnegative(time, name="time")
    S = _exact_positive(stationary_mean, name="stationary_mean")
    r = _exact_positive(rate, name="rate")
    initial = sp.sympify(initial_generating_function)
    if initial.has(sp.Float):
        raise ValueError("initial_generating_function must be exact")
    survival = sp.exp(-r * t)
    transported_variable = 1 + (variable - 1) * survival
    immigrants = sp.exp(S * (variable - 1) * (1 - survival))
    return sp.simplify(initial.subs(variable, transported_variable) * immigrants)


def immigration_death_transition_probability(
    final_state: Any,
    *,
    initial_state: Any,
    time: Any,
    stationary_mean: Any,
    rate: Any,
) -> sp.Expr:
    r"""Return the exact deterministic-initial transition probability.

    At time ``t``, the survivors of ``n0`` initial individuals have law
    ``Binomial(n0, exp(-r*t))`` and surviving immigrants have the independent
    law ``Poisson(S*(1-exp(-r*t)))``.  This finite convolution is their sum.
    """

    destination = _nonnegative_integer(final_state, name="final_state")
    origin = _nonnegative_integer(initial_state, name="initial_state")
    t = _exact_nonnegative(time, name="time")
    S = _exact_positive(stationary_mean, name="stationary_mean")
    r = _exact_positive(rate, name="rate")
    survival = sp.exp(-r * t)
    immigrant_mean = S * (1 - survival)
    probability = sp.Integer(0)
    for survivors in range(min(origin, destination) + 1):
        binomial_mass = (
            sp.binomial(origin, survivors)
            * survival**survivors
            * (1 - survival) ** (origin - survivors)
        )
        immigrant_order = destination - survivors
        poisson_mass = (
            sp.exp(-immigrant_mean)
            * immigrant_mean**immigrant_order
            / sp.factorial(immigrant_order)
        )
        probability += binomial_mass * poisson_mass
    return sp.simplify(probability)


def reversible_factorial_one_rates(
    state: Any,
    *,
    stationary_mean: Any,
    rate: Any,
) -> tuple[sp.Expr, sp.Expr]:
    r"""Return a distinct reversible chain with the same stationary mass.

    The alternative rates are ``lambda_n=r*S/(n+1)`` and ``mu_n=r`` for
    ``n>=1``, with ``mu_0=0``.  They have the same detailed-balance ratio
    ``lambda_n/mu_(n+1)=S/(n+1)`` as the immigration--death chain but different
    holding times and transients.  This is a constructive nonuniqueness
    witness against inferring dynamics from the stationary mass alone.
    """

    n = _nonnegative_integer(state, name="state")
    S = _exact_positive(stationary_mean, name="stationary_mean")
    r = _exact_positive(rate, name="rate")
    death = sp.Integer(0) if n == 0 else r
    return sp.simplify(r * S / (n + 1)), sp.simplify(death)


@dataclass(frozen=True)
class ImmigrationDeathLedger:
    """Exact rates, stationary data, drift, and interpretation ceilings."""

    state: int
    stationary_mean: sp.Expr
    rate: sp.Expr
    birth_rate: sp.Expr
    death_rate: sp.Expr
    local_drift: sp.Expr
    stationary_mass: sp.Expr
    adjacent_stationary_ratio: sp.Expr
    boundary_death_rate_is_zero: bool
    static_mass_selects_unique_dynamics: bool
    positive_drift_implies_monotone_sample_paths: bool
    material_process_is_separate_premise: bool


def immigration_death_ledger(
    state: Any,
    *,
    stationary_mean: Any,
    rate: Any,
) -> ImmigrationDeathLedger:
    """Return one exact state ledger for the declared Markov chain."""

    n = _nonnegative_integer(state, name="state")
    S = _exact_positive(stationary_mean, name="stationary_mean")
    r = _exact_positive(rate, name="rate")
    birth, death = immigration_death_rates(n, stationary_mean=S, rate=r)
    return ImmigrationDeathLedger(
        state=n,
        stationary_mean=S,
        rate=r,
        birth_rate=birth,
        death_rate=death,
        local_drift=immigration_death_local_drift(
            n,
            stationary_mean=S,
            rate=r,
        ),
        stationary_mass=immigration_death_stationary_mass(
            n,
            stationary_mean=S,
        ),
        adjacent_stationary_ratio=sp.simplify(S / (n + 1)),
        boundary_death_rate_is_zero=immigration_death_rates(
            0,
            stationary_mean=S,
            rate=r,
        )[1]
        == 0,
        static_mass_selects_unique_dynamics=False,
        positive_drift_implies_monotone_sample_paths=False,
        material_process_is_separate_premise=True,
    )
