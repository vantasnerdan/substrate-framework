"""Exact finite-dimensional SU(3) representation data from Dynkin labels.

The convention is the standard isospin embedding with Hermitian generators
``T_a=lambda_a/2`` and hypercharge ``Y=2*T_8/sqrt(3)``.  A representation with
Dynkin labels ``(p, q)`` uses the U(3) Gelfand--Tsetlin top row
``(p + q, q, 0)``.  No collective-coordinate or particle interpretation is
attached to the resulting mathematical weights.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from typing import Any

import sympy as sp


def _nonnegative_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, sp.Integer)):
        raise TypeError(f"{name} must be an integer")
    normalized = int(value)
    if normalized < 0:
        raise ValueError(f"{name} must be nonnegative")
    return normalized


def _exact_rational(value: Any, name: str) -> sp.Rational:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact rational number")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    if isinstance(value, (int, sp.Integer, sp.Rational)):
        result = sp.Rational(value)
        if result.is_Rational:
            return result
    raise TypeError(f"{name} must be an exact rational number")


@dataclass(frozen=True)
class SU3GelfandTsetlinState:
    """One Gelfand--Tsetlin basis state and its exact SU(2)xU(1) labels."""

    p: int
    q: int
    m12: int
    m22: int
    m11: int
    isospin: sp.Rational
    isospin_projection: sp.Rational
    hypercharge: sp.Rational


@dataclass(frozen=True)
class SU3WeightMultiplicity:
    """Multiplicity of one ``(I3, Y)`` weight in an SU(3) irrep."""

    isospin_projection: sp.Rational
    hypercharge: sp.Rational
    multiplicity: int


@dataclass(frozen=True)
class SU3IsospinMultiplet:
    """One SU(2) isospin multiplet at fixed hypercharge."""

    isospin: sp.Rational
    hypercharge: sp.Rational


@dataclass(frozen=True)
class SU3Irrep:
    """Exact mathematical SU(3) irrep specified by nonnegative ``(p, q)``."""

    p: int
    q: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "p", _nonnegative_integer(self.p, "p"))
        object.__setattr__(self, "q", _nonnegative_integer(self.q, "q"))

    @property
    def dimension(self) -> int:
        return su3_irrep_dimension(self.p, self.q)

    @property
    def quadratic_casimir(self) -> sp.Rational:
        return su3_irrep_quadratic_casimir(self.p, self.q)

    @property
    def triality(self) -> int:
        return su3_irrep_triality(self.p, self.q)

    @property
    def states(self) -> tuple[SU3GelfandTsetlinState, ...]:
        return su3_gelfand_tsetlin_states(self.p, self.q)

    @property
    def weights(self) -> tuple[SU3WeightMultiplicity, ...]:
        return su3_weight_multiplicities(self.p, self.q)

    @property
    def isospin_multiplets(self) -> tuple[SU3IsospinMultiplet, ...]:
        return su3_isospin_multiplets(self.p, self.q)

    def multiplets_at_hypercharge(
        self, hypercharge: Any
    ) -> tuple[SU3IsospinMultiplet, ...]:
        target = _exact_rational(hypercharge, "hypercharge")
        return tuple(
            multiplet
            for multiplet in self.isospin_multiplets
            if multiplet.hypercharge == target
        )

    def contains_hypercharge(self, hypercharge: Any) -> bool:
        return bool(self.multiplets_at_hypercharge(hypercharge))


@dataclass(frozen=True)
class SU3HyperchargeMatch:
    """One irrep selected by a supplied mathematical hypercharge filter."""

    p: int
    q: int
    dimension: int
    quadratic_casimir: sp.Rational
    triality: int
    hypercharge: sp.Rational
    isospins: tuple[sp.Rational, ...]


def su3_irrep_dimension(p: Any, q: Any) -> int:
    """Return the Weyl dimension of the SU(3) irrep ``(p, q)``."""

    p_value = _nonnegative_integer(p, "p")
    q_value = _nonnegative_integer(q, "q")
    numerator = (p_value + 1) * (q_value + 1) * (p_value + q_value + 2)
    if numerator % 2:
        raise ArithmeticError("SU(3) Weyl numerator must be even")
    return numerator // 2


def su3_irrep_quadratic_casimir(p: Any, q: Any) -> sp.Rational:
    """Return ``C2(p,q)`` in the ``T_a=lambda_a/2`` normalization."""

    p_value = _nonnegative_integer(p, "p")
    q_value = _nonnegative_integer(q, "q")
    return sp.Rational(
        p_value**2
        + p_value * q_value
        + q_value**2
        + 3 * p_value
        + 3 * q_value,
        3,
    )


def su3_irrep_triality(p: Any, q: Any) -> int:
    """Return the center triality ``p + 2*q (mod 3)``."""

    p_value = _nonnegative_integer(p, "p")
    q_value = _nonnegative_integer(q, "q")
    return (p_value + 2 * q_value) % 3


@cache
def _gelfand_tsetlin_states(
    p: int, q: int
) -> tuple[SU3GelfandTsetlinState, ...]:
    states: list[SU3GelfandTsetlinState] = []
    hypercharge_offset = sp.Rational(2 * (p + 2 * q), 3)
    for m12 in range(q, p + q + 1):
        for m22 in range(0, q + 1):
            isospin = sp.Rational(m12 - m22, 2)
            hypercharge = sp.Rational(m12 + m22) - hypercharge_offset
            for m11 in range(m22, m12 + 1):
                states.append(
                    SU3GelfandTsetlinState(
                        p=p,
                        q=q,
                        m12=m12,
                        m22=m22,
                        m11=m11,
                        isospin=isospin,
                        isospin_projection=sp.Rational(
                            2 * m11 - m12 - m22, 2
                        ),
                        hypercharge=hypercharge,
                    )
                )
    expected = su3_irrep_dimension(p, q)
    if len(states) != expected:
        raise ArithmeticError(
            f"Gelfand--Tsetlin count {len(states)} does not equal dimension {expected}"
        )
    return tuple(states)


def su3_gelfand_tsetlin_states(
    p: Any, q: Any
) -> tuple[SU3GelfandTsetlinState, ...]:
    """Enumerate every interlacing Gelfand--Tsetlin basis state."""

    return _gelfand_tsetlin_states(
        _nonnegative_integer(p, "p"), _nonnegative_integer(q, "q")
    )


@cache
def _weight_multiplicities(
    p: int, q: int
) -> tuple[SU3WeightMultiplicity, ...]:
    counts = Counter(
        (state.isospin_projection, state.hypercharge)
        for state in _gelfand_tsetlin_states(p, q)
    )
    ordered = sorted(counts.items(), key=lambda item: (-item[0][1], -item[0][0]))
    return tuple(
        SU3WeightMultiplicity(
            isospin_projection=weight[0],
            hypercharge=weight[1],
            multiplicity=multiplicity,
        )
        for weight, multiplicity in ordered
    )


def su3_weight_multiplicities(
    p: Any, q: Any
) -> tuple[SU3WeightMultiplicity, ...]:
    """Return all distinct ``(I3,Y)`` weights with exact multiplicity."""

    return _weight_multiplicities(
        _nonnegative_integer(p, "p"), _nonnegative_integer(q, "q")
    )


@cache
def _isospin_multiplets(p: int, q: int) -> tuple[SU3IsospinMultiplet, ...]:
    hypercharge_offset = sp.Rational(2 * (p + 2 * q), 3)
    multiplets = tuple(
        SU3IsospinMultiplet(
            isospin=sp.Rational(m12 - m22, 2),
            hypercharge=sp.Rational(m12 + m22) - hypercharge_offset,
        )
        for m12 in range(q, p + q + 1)
        for m22 in range(0, q + 1)
    )
    return tuple(
        sorted(
            multiplets,
            key=lambda multiplet: (-multiplet.hypercharge, -multiplet.isospin),
        )
    )


def su3_isospin_multiplets(
    p: Any, q: Any
) -> tuple[SU3IsospinMultiplet, ...]:
    """Return the multiplicity-free SU(3) to SU(2)xU(1) branching rows."""

    return _isospin_multiplets(
        _nonnegative_integer(p, "p"), _nonnegative_integer(q, "q")
    )


def irreps_containing_hypercharge(
    hypercharge: Any, *, max_p: Any, max_q: Any
) -> tuple[SU3HyperchargeMatch, ...]:
    """Enumerate every matching irrep in the declared rectangular label domain.

    This is a mathematical weight filter only.  It neither chooses among equal
    dimensions nor assigns collective, spin, statistics, or particle meaning.
    """

    target = _exact_rational(hypercharge, "hypercharge")
    p_limit = _nonnegative_integer(max_p, "max_p")
    q_limit = _nonnegative_integer(max_q, "max_q")
    matches: list[SU3HyperchargeMatch] = []
    for p in range(p_limit + 1):
        for q in range(q_limit + 1):
            irrep = SU3Irrep(p, q)
            selected = irrep.multiplets_at_hypercharge(target)
            if selected:
                matches.append(
                    SU3HyperchargeMatch(
                        p=p,
                        q=q,
                        dimension=irrep.dimension,
                        quadratic_casimir=irrep.quadratic_casimir,
                        triality=irrep.triality,
                        hypercharge=target,
                        isospins=tuple(
                            multiplet.isospin for multiplet in selected
                        ),
                    )
                )
    return tuple(sorted(matches, key=lambda match: (match.dimension, match.p, match.q)))


def minimal_dimension_hypercharge_matches(
    hypercharge: Any, *, max_p: Any, max_q: Any
) -> tuple[SU3HyperchargeMatch, ...]:
    """Return all minimum-dimension ties inside an explicit search domain."""

    matches = irreps_containing_hypercharge(
        hypercharge, max_p=max_p, max_q=max_q
    )
    if not matches:
        return ()
    minimum = matches[0].dimension
    return tuple(match for match in matches if match.dimension == minimum)
