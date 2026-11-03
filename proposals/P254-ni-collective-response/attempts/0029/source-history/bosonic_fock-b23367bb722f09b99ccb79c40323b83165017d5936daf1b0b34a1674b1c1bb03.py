"""Exact one-mode bosonic Fock algebra and conditional cosine composition.

The infinite ladder identities are stated on the algebraic finite-support
span of an orthonormal occupation basis.  Finite matrices are explicitly
truncations and carry a top-state commutator defect.  The conditional cosine
helpers require the declared low-coordinate convention
``Q=q_0*(a+a_dagger)``.  They return the low-sector matrix element of the
H-linear coefficient, not a complete high-to-low transition amplitude; no
density of states, transition probability, or physical rate follows.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import sympy as sp

from .cosine_vertices import vacuum_one_high_coefficient


FactorialOneSupport = Literal[
    "all_nonnegative",
    "positive",
    "positive_odd",
]
_FACTORIAL_ONE_SUPPORTS = {
    "all_nonnegative",
    "positive",
    "positive_odd",
}


def _nonnegative_integer(value: Any, *, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be an exact nonnegative integer")
    return int(expression)


def _positive_integer(value: Any, *, name: str) -> int:
    integer = _nonnegative_integer(value, name=name)
    if integer == 0:
        raise ValueError(f"{name} must be a positive integer")
    return integer


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


def _exact_positive_rational(value: Any, *, name: str) -> sp.Rational:
    expression = sp.sympify(value)
    if (
        expression.has(sp.Float)
        or expression.is_Rational is not True
        or expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be an exact positive rational")
    return sp.Rational(expression)


def _support(value: str) -> FactorialOneSupport:
    if value not in _FACTORIAL_ONE_SUPPORTS:
        raise ValueError(
            "support must be 'all_nonnegative', 'positive', or 'positive_odd'"
        )
    return value  # type: ignore[return-value]


@dataclass(frozen=True)
class BosonicFockRung:
    """Exact one-level data on the normalized infinite Fock basis.

    The basis convention is ``<m|n>=delta_mn`` with
    ``a|n>=sqrt(n)|n-1>`` and ``a_dagger|n>=sqrt(n+1)|n+1>``.  These
    identities define mathematical operators on the algebraic finite-support
    domain; they do not establish a physical oscillator realization.
    """

    level: int
    annihilation_coefficient: sp.Expr
    creation_coefficient: sp.Expr
    number_coefficient: sp.Integer
    reverse_number_coefficient: sp.Integer
    commutator_coefficient: sp.Integer
    repeated_creation_vacuum_coefficient: sp.Expr
    repeated_creation_vacuum_norm_squared: sp.Integer
    physical_rate_interpretation_is_separate_premise: bool


def bosonic_fock_rung(level: Any) -> BosonicFockRung:
    r"""Return the exact normalized ladder data at occupation ``n``.

    In particular,
    ``(a_dagger)**n|0> = sqrt(n!)|n>`` and its squared norm is ``n!``.
    The normalized occupation state therefore divides the repeated-creation
    vector by ``sqrt(n!)``; the factorial is not a count of distinct
    single-mode final states.
    """

    n = _nonnegative_integer(level, name="level")
    return BosonicFockRung(
        level=n,
        annihilation_coefficient=sp.sqrt(n),
        creation_coefficient=sp.sqrt(n + 1),
        number_coefficient=sp.Integer(n),
        reverse_number_coefficient=sp.Integer(n + 1),
        commutator_coefficient=sp.Integer(1),
        repeated_creation_vacuum_coefficient=sp.sqrt(sp.factorial(n)),
        repeated_creation_vacuum_norm_squared=sp.factorial(n),
        physical_rate_interpretation_is_separate_premise=True,
    )


@dataclass(frozen=True)
class TruncatedBosonicFockLadder:
    """Exact matrices and top-state defect for a D-level truncation."""

    dimension: int
    annihilation: sp.ImmutableMatrix
    creation: sp.ImmutableMatrix
    commutator: sp.ImmutableMatrix
    identity: sp.ImmutableMatrix
    top_projector: sp.ImmutableMatrix
    expected_commutator: sp.ImmutableMatrix
    identity_minus_commutator: sp.ImmutableMatrix
    commutator_trace: sp.Integer
    full_identity_commutator_is_impossible: bool


def truncated_bosonic_fock_ladder(
    dimension: Any,
) -> TruncatedBosonicFockLadder:
    r"""Return the exact D-level ladder and its unavoidable edge defect.

    The matrices obey
    ``[a_D,a_D_dagger] = I_D-D|D-1><D-1|``.  Thus the block below the top
    state is the identity, while the full finite commutator is not.  This also
    realizes the trace obstruction: every finite commutator has trace zero,
    whereas ``trace(I_D)=D``.
    """

    count = _positive_integer(dimension, name="dimension")
    annihilation = sp.zeros(count)
    for source_level in range(1, count):
        annihilation[source_level - 1, source_level] = sp.sqrt(source_level)
    creation = annihilation.T
    identity = sp.eye(count)
    top_projector = sp.zeros(count)
    top_projector[count - 1, count - 1] = 1
    commutator = sp.simplify(annihilation * creation - creation * annihilation)
    expected = identity - count * top_projector
    return TruncatedBosonicFockLadder(
        dimension=count,
        annihilation=sp.ImmutableMatrix(annihilation),
        creation=sp.ImmutableMatrix(creation),
        commutator=sp.ImmutableMatrix(commutator),
        identity=sp.ImmutableMatrix(identity),
        top_projector=sp.ImmutableMatrix(top_projector),
        expected_commutator=sp.ImmutableMatrix(expected),
        identity_minus_commutator=sp.ImmutableMatrix(identity - commutator),
        commutator_trace=sp.Integer(sp.trace(commutator)),
        full_identity_commutator_is_impossible=True,
    )


def bosonic_cosine_matrix_element(
    low_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
) -> sp.Expr:
    r"""Conditionally compose a vacuum cosine coefficient with ``sqrt(n!)``.

    The caller declares ``Q=low_scale*(a+a_dagger)`` on the normalized
    one-mode basis.  In ``<n|Q**n|0>`` only the all-creation word can reach
    level ``n``, giving ``low_scale**n*sqrt(n!)``.  The H-linear factor remains
    a formal classical coefficient: a high-sector operator and state element
    are separately required for a complete transition amplitude.  At zero
    background this returns zero for even ``n`` and
    ``amplitude*(-1)**((n-1)/2)*high_scale*low_scale**n/sqrt(n!)`` for odd
    ``n``.  It is a conditional low-sector algebraic element, not a rate.
    """

    n = _nonnegative_integer(low_order, name="low_order")
    amplitude_expression = _exact_real(amplitude, name="amplitude")
    high_expression = _exact_real(high_scale, name="high_scale")
    low_expression = _exact_real(low_scale, name="low_scale")
    coefficient = vacuum_one_high_coefficient(
        n,
        amplitude=amplitude_expression,
        high_scale=high_expression,
        low_scale=low_expression,
    )
    return sp.simplify(coefficient * sp.sqrt(sp.factorial(n)))


def bosonic_cosine_matrix_element_square(
    low_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
) -> sp.Expr:
    """Return the exact square of the declared real algebraic element."""

    element = bosonic_cosine_matrix_element(
        low_order,
        amplitude=amplitude,
        high_scale=high_scale,
        low_scale=low_scale,
    )
    return sp.factor(element**2)


def factorial_one_mass(
    order: Any,
    *,
    intensity: Any,
    support: FactorialOneSupport,
) -> sp.Expr:
    r"""Return ``S**n/n!`` on one explicitly declared integer support.

    The ambient domain is the nonnegative integers.  The support is either
    all nonnegative integers, positive integers, or positive odd integers;
    off-support values are exactly zero.  This is a mathematical mass, not a
    physical occurrence law or rate.
    """

    n = _nonnegative_integer(order, name="order")
    S = _exact_positive(intensity, name="intensity")
    declared_support = _support(support)
    included = (
        declared_support == "all_nonnegative"
        or (declared_support == "positive" and n >= 1)
        or (declared_support == "positive_odd" and n >= 1 and n % 2 == 1)
    )
    if not included:
        return sp.Integer(0)
    return sp.factor(S**n / sp.factorial(n))


def factorial_one_total_mass(
    *,
    intensity: Any,
    support: FactorialOneSupport,
) -> sp.Expr:
    """Return the exact total mass for the declared factorial-one support."""

    S = _exact_positive(intensity, name="intensity")
    declared_support = _support(support)
    if declared_support == "all_nonnegative":
        return sp.exp(S)
    if declared_support == "positive":
        return sp.exp(S) - 1
    return sp.sinh(S)


def normalized_factorial_one_mass(
    order: Any,
    *,
    intensity: Any,
    support: FactorialOneSupport,
) -> sp.Expr:
    """Return one exact normalized point mass on the declared support."""

    mass = factorial_one_mass(order, intensity=intensity, support=support)
    if mass == 0:
        return sp.Integer(0)
    total = factorial_one_total_mass(intensity=intensity, support=support)
    return sp.factor(mass / total)


def factorial_one_modes(
    *,
    intensity: Any,
    support: FactorialOneSupport,
) -> tuple[int, ...]:
    """Return every exact mode for a positive rational intensity.

    For the all-nonnegative and positive supports, integer intensity produces
    the familiar adjacent tie except that excluding zero leaves intensity one
    with the single positive mode one.  On positive odd support the ratio of
    adjacent supported masses is ``S**2/((n+1)*(n+2))``; the returned tuple
    preserves any exact tie rather than choosing the first grid occurrence.
    """

    S = _exact_positive_rational(intensity, name="intensity")
    declared_support = _support(support)
    if declared_support in {"all_nonnegative", "positive"}:
        if S.q == 1:
            integer = int(S)
            if declared_support == "positive" and integer == 1:
                return (1,)
            return (integer - 1, integer)
        floor = int(sp.floor(S))
        if declared_support == "positive":
            floor = max(1, floor)
        return (floor,)

    squared = S**2
    mode = 1
    while squared > (mode + 1) * (mode + 2):
        mode += 2
    if squared == (mode + 1) * (mode + 2):
        return (mode, mode + 2)
    return (mode,)


def factorial_one_log_concavity_ratio(order: Any) -> sp.Rational:
    r"""Return ``w(n)**2/(w(n-1)*w(n+1))=(n+1)/n`` for ``n>=1``.

    The ratio is strictly greater than one and independent of the positive
    intensity.  Strict discrete log-concavity implies that the modes are
    consecutive, but it does not exclude the adjacent tie already returned by
    :func:`factorial_one_modes` when the intensity is an integer.
    """

    n = _positive_integer(order, name="order")
    return sp.Rational(n + 1, n)


def factorial_one_probability_generating_function(
    *,
    intensity: Any,
    variable: Any,
) -> sp.Expr:
    r"""Return the exact PGF ``exp(S*(t-1))`` of the all-order mass.

    This is the generating function of the normalized mathematical mass on
    the nonnegative integers.  It does not assert a physical Poisson process,
    event count, rate, time interval, or medium realization.
    """

    S = _exact_positive(intensity, name="intensity")
    t = _exact_real(variable, name="variable")
    return sp.exp(S * (t - 1))


def factorial_one_falling_factorial_moment(
    moment_order: Any,
    *,
    intensity: Any,
) -> sp.Expr:
    r"""Return the exact all-order falling-factorial moment ``E[(N)_r]``.

    For every nonnegative integer ``r`` this mathematical moment is ``S**r``.
    In particular the mean and variance are both ``S`` after converting the
    first two falling-factorial moments to raw moments.
    """

    r = _nonnegative_integer(moment_order, name="moment_order")
    S = _exact_positive(intensity, name="intensity")
    return sp.factor(S**r)


def _geometric_tail_inputs(
    *,
    intensity: Any,
    alpha: Any,
    starting_order: Any,
) -> tuple[sp.Expr, sp.Expr, int]:
    S = _exact_positive(intensity, name="intensity")
    exponent = _exact_positive(alpha, name="alpha")
    start = _nonnegative_integer(starting_order, name="starting_order")
    threshold_residual = sp.simplify(start + 1 - S * sp.exp(exponent))
    if threshold_residual.is_nonnegative is not True:
        raise ValueError(
            "starting_order + 1 must be at least intensity*exp(alpha) exactly"
        )
    return S, exponent, start


def factorial_one_geometric_point_bound(
    steps: Any,
    *,
    intensity: Any,
    alpha: Any,
    starting_order: Any,
) -> sp.Expr:
    r"""Bound ``p(N+k)`` by ``p(N)*exp(-alpha*k)`` after the threshold.

    The exact premises are ``S>0``, ``alpha>0``, ``N,k>=0`` integral, and
    ``N+1>=S*exp(alpha)``.  Every later adjacent ratio is then at most
    ``exp(-alpha)``.  The result concerns the normalized all-nonnegative
    mathematical mass only.
    """

    k = _nonnegative_integer(steps, name="steps")
    S, exponent, start = _geometric_tail_inputs(
        intensity=intensity,
        alpha=alpha,
        starting_order=starting_order,
    )
    initial = normalized_factorial_one_mass(
        start,
        intensity=S,
        support="all_nonnegative",
    )
    return sp.factor(initial * sp.exp(-exponent * k))


def factorial_one_geometric_tail_bound(
    *,
    intensity: Any,
    alpha: Any,
    starting_order: Any,
) -> sp.Expr:
    r"""Bound the mass strictly above ``N`` by ``p(N)/(exp(alpha)-1)``.

    This sums the geometric point bounds for steps one through infinity under
    the same exact threshold premises as
    :func:`factorial_one_geometric_point_bound`.
    """

    S, exponent, start = _geometric_tail_inputs(
        intensity=intensity,
        alpha=alpha,
        starting_order=starting_order,
    )
    initial = normalized_factorial_one_mass(
        start,
        intensity=S,
        support="all_nonnegative",
    )
    return sp.factor(initial / (sp.exp(exponent) - 1))


@dataclass(frozen=True)
class FactorialOnePolynomialTailCertificate:
    """An explicit eventual contraction for ``n**r*p_S(n)``.

    The certificate uses ``(1+1/n)**r<=2**r`` for ``n>=1``.  It is a
    conservative exact threshold, not a fitted asymptotic estimate.
    """

    intensity: sp.Rational
    power: int
    contraction: sp.Rational
    starting_order: int
    threshold: sp.Rational
    ratio_bound: sp.Rational
    scaled_mass_tends_to_zero: bool
    physical_power_law_interpretation_is_separate_premise: bool


def factorial_one_polynomial_tail_certificate(
    power: Any,
    *,
    intensity: Any,
    contraction: Any = sp.Rational(1, 2),
) -> FactorialOnePolynomialTailCertificate:
    r"""Certify geometric decay of ``n**r*p_S(n)`` beyond an exact index.

    For positive rational ``S``, nonnegative integer ``r``, and rational
    ``q`` in ``(0,1)``, choosing ``N+1>=S*2**r/q`` makes the consecutive
    scaled-mass ratio at most ``q`` for every ``n>=N``.  Therefore
    ``n**r*p_S(n)`` tends to zero for every fixed ``r``.
    """

    r = _nonnegative_integer(power, name="power")
    S = _exact_positive_rational(intensity, name="intensity")
    q = _exact_positive_rational(contraction, name="contraction")
    if q >= 1:
        raise ValueError("contraction must be strictly less than one")
    threshold = sp.Rational(S * 2**r / q)
    start = max(1, int(sp.ceiling(threshold)) - 1)
    return FactorialOnePolynomialTailCertificate(
        intensity=S,
        power=r,
        contraction=q,
        starting_order=start,
        threshold=threshold,
        ratio_bound=q,
        scaled_mass_tends_to_zero=True,
        physical_power_law_interpretation_is_separate_premise=True,
    )
