"""Exact conditional geometry for rectangular flat tori.

This module is unpromoted infrastructure harvested by proposal P237 while
auditing ``Viktar-Pi/FlatIrrationalTorus`` at commit
``ba487f6361f1f8740a4fe518bf86584847ab518b`` (MIT).  It cleanly separates
the reusable flat-torus mathematics from that source's particle, gravity, and
cosmology interpretations.

The declared space is the rectangular quotient
``R^d / (L_1 Z x ... x L_d Z)``, with every ``L_i > 0``.  A complex scalar
with boundary phase
``phi(x + L_i e_i) = exp(2*pi*I*theta_i) phi(x)`` has reciprocal modes
``k_i = 2*pi*(n_i + theta_i)/L_i`` for ``n_i in Z``.  The scalar operator
``-Delta`` therefore has eigenvalue ``sum_i k_i**2``.  No spin structure,
Dirac multiplicity, spectral action, or particle-mass map follows from that
scalar identity.

The matched-circle helpers use only Euclidean sphere intersections in the
universal cover.  They do not specify an observer likelihood, transfer
function, last-scattering dynamics, or CMB isotropy.  In particular, this
module implements ``T^d`` only: a notation such as ``T^3/Z_2`` requires a
separately supplied group action, fixed-point treatment, field lift, and
operator domain.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from itertools import product
from typing import Any

from mpmath import mp, mpf
import sympy as sp


@dataclass(frozen=True)
class FlatTorusMode:
    """One integer-labelled reciprocal mode and its scalar eigenvalue."""

    index: tuple[int, ...]
    wavevector: tuple[sp.Expr, ...]
    eigenvalue: sp.Expr


@dataclass(frozen=True)
class LaplacianEigenspace:
    """Exactly coincident scalar-Laplacian modes in a finite index cube."""

    eigenvalue: sp.Expr
    modes: tuple[tuple[int, ...], ...]

    @property
    def multiplicity(self) -> int:
        """Return the number of enumerated integer modes in the eigenspace."""

        return len(self.modes)


@dataclass(frozen=True)
class LatticeTranslation:
    """One deck translation of a rectangular torus."""

    index: tuple[int, ...]
    vector: tuple[sp.Expr, ...]
    squared_length: sp.Expr

    @property
    def length(self) -> sp.Expr:
        """Return the nonnegative Euclidean translation length."""

        return sp.sqrt(self.squared_length)


@dataclass(frozen=True)
class MatchedCircleGeometry:
    """A nondegenerate sphere-intersection circle for one deck translation."""

    translation: LatticeTranslation
    angular_radius: sp.Expr


@dataclass(frozen=True)
class EpsteinRefinementPoint:
    """One cubic-cutoff partial sum and its increment from the prior cutoff."""

    max_index: int
    partial_sum: mpf
    increment: mpf | None


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and (
        expression.is_real is not True or expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be positive and real")
    if expression.is_real is False or expression.is_positive is False:
        raise ValueError(f"{name} must not be known nonpositive or nonreal")
    return expression


def _positive_numeric(value: Any, name: str) -> sp.Expr:
    expression = _positive(value, name)
    if expression.is_number is not True:
        raise ValueError(f"{name} must be numeric")
    return expression


def _side_lengths(side_lengths: Sequence[Any]) -> tuple[sp.Expr, ...]:
    if isinstance(side_lengths, (str, bytes)) or len(side_lengths) == 0:
        raise ValueError("side_lengths must be a nonempty sequence")
    return tuple(
        _positive(value, f"side_lengths[{axis}]")
        for axis, value in enumerate(side_lengths)
    )


def _numeric_side_lengths(side_lengths: Sequence[Any]) -> tuple[sp.Expr, ...]:
    sides = _side_lengths(side_lengths)
    if any(side.is_number is not True for side in sides):
        raise ValueError("side_lengths must be numeric for enumeration")
    return sides


def _integer_vector(
    values: Sequence[Any], dimension: int, name: str
) -> tuple[int, ...]:
    if isinstance(values, (str, bytes)) or len(values) != dimension:
        raise ValueError(f"{name} must have length {dimension}")
    result: list[int] = []
    for axis, value in enumerate(values):
        expression = sp.sympify(value)
        if expression.is_number is not True or expression.is_integer is not True:
            raise ValueError(f"{name}[{axis}] must be an integer")
        result.append(int(expression))
    return tuple(result)


def _twist_vector(
    twist: Sequence[Any] | None, dimension: int
) -> tuple[sp.Expr, ...]:
    if twist is None:
        return (sp.S.Zero,) * dimension
    if isinstance(twist, (str, bytes)) or len(twist) != dimension:
        raise ValueError(f"twist must have length {dimension}")
    result: list[sp.Expr] = []
    for axis, value in enumerate(twist):
        expression = sp.sympify(value)
        if expression.is_number and expression.is_real is not True:
            raise ValueError(f"twist[{axis}] must be real")
        if expression.is_real is False:
            raise ValueError(f"twist[{axis}] must not be known nonreal")
        result.append(expression)
    return tuple(result)


def _nonnegative_integer(value: Any, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(expression)


def _strictly_less(left: sp.Expr, right: sp.Expr, description: str) -> bool:
    relation = sp.StrictLessThan(left, right)
    if relation is sp.S.true:
        return True
    if relation is sp.S.false:
        return False
    raise ValueError(f"{description} must be numerically decidable")


def rectangular_torus_volume(side_lengths: Sequence[Any]) -> sp.Expr:
    """Return the exact fundamental-cell volume ``prod_i L_i``."""

    return sp.prod(_side_lengths(side_lengths))


def reciprocal_wavevector(
    index: Sequence[Any],
    side_lengths: Sequence[Any],
    *,
    twist: Sequence[Any] | None = None,
) -> tuple[sp.Expr, ...]:
    """Return ``2*pi*(n_i+theta_i)/L_i`` component by component.

    ``twist=None`` gives periodic modes.  A half twist must be supplied as
    ``Rational(1, 2)`` (or an equivalent real value); it is not obtained by
    iterating over all integer and half-integer labels in one spectrum.
    """

    sides = _side_lengths(side_lengths)
    modes = _integer_vector(index, len(sides), "index")
    phases = _twist_vector(twist, len(sides))
    return tuple(
        sp.simplify(2 * sp.pi * (mode + phase) / side)
        for mode, phase, side in zip(modes, phases, sides, strict=True)
    )


def scalar_laplacian_eigenvalue(
    index: Sequence[Any],
    side_lengths: Sequence[Any],
    *,
    twist: Sequence[Any] | None = None,
) -> sp.Expr:
    """Return the exact ``-Delta`` eigenvalue for one supplied scalar mode."""

    wavevector = reciprocal_wavevector(index, side_lengths, twist=twist)
    return sp.simplify(sum(component**2 for component in wavevector))


def enumerate_laplacian_modes(
    side_lengths: Sequence[Any],
    max_index: Any,
    *,
    twist: Sequence[Any] | None = None,
) -> tuple[FlatTorusMode, ...]:
    """Enumerate the complete integer cube ``-N <= n_i <= N``.

    The result includes the periodic zero mode when it is present.  This is a
    finite enumeration API, not a claim about a continuum density of states.
    """

    sides = _side_lengths(side_lengths)
    cutoff = _nonnegative_integer(max_index, "max_index")
    phases = _twist_vector(twist, len(sides))
    modes: list[FlatTorusMode] = []
    for index in product(range(-cutoff, cutoff + 1), repeat=len(sides)):
        wavevector = reciprocal_wavevector(index, sides, twist=phases)
        modes.append(
            FlatTorusMode(
                index=tuple(index),
                wavevector=wavevector,
                eigenvalue=sp.simplify(
                    sum(component**2 for component in wavevector)
                ),
            )
        )
    return tuple(modes)


def laplacian_eigenspaces(
    side_lengths: Sequence[Any],
    max_index: Any,
    *,
    twist: Sequence[Any] | None = None,
) -> tuple[LaplacianEigenspace, ...]:
    """Group exactly equal eigenvalues in a finite integer-mode cube."""

    grouped: dict[sp.Expr, list[tuple[int, ...]]] = {}
    for mode in enumerate_laplacian_modes(
        side_lengths, max_index, twist=twist
    ):
        key = sp.simplify(mode.eigenvalue)
        grouped.setdefault(key, []).append(mode.index)
    def eigenvalue_sort_key(item: tuple[sp.Expr, list[tuple[int, ...]]]) -> Any:
        value = item[0]
        if value.is_number is True and value.is_real is True:
            return (0, float(sp.N(value, 50)))
        return (1, sp.default_sort_key(value))

    return tuple(
        LaplacianEigenspace(eigenvalue=value, modes=tuple(sorted(indices)))
        for value, indices in sorted(grouped.items(), key=eigenvalue_sort_key)
    )


def translation_vector(
    index: Sequence[Any], side_lengths: Sequence[Any]
) -> tuple[sp.Expr, ...]:
    """Return the deck-translation vector ``(m_i L_i)``."""

    sides = _side_lengths(side_lengths)
    lattice_index = _integer_vector(index, len(sides), "index")
    return tuple(
        sp.simplify(mode * side)
        for mode, side in zip(lattice_index, sides, strict=True)
    )


def translation_squared_length(
    index: Sequence[Any], side_lengths: Sequence[Any]
) -> sp.Expr:
    """Return ``sum_i (m_i L_i)^2`` for one deck translation."""

    vector = translation_vector(index, side_lengths)
    return sp.simplify(sum(component**2 for component in vector))


def shortest_translation_length(side_lengths: Sequence[Any]) -> sp.Expr:
    """Return the shortest nonzero deck-translation length ``min_i L_i``."""

    return sp.Min(*_side_lengths(side_lengths))


def enumerate_translations_below(
    side_lengths: Sequence[Any],
    maximum_length: Any,
    *,
    unique_up_to_sign: bool = False,
) -> tuple[LatticeTranslation, ...]:
    """Enumerate every nonzero deck translation with length strictly below a bound.

    The coordinate search bounds are complete because
    ``|m_i| L_i <= |m L| < maximum_length``.  If ``unique_up_to_sign`` is
    true, the first nonzero integer component is required to be positive, so
    exactly one member of each ``m, -m`` pair is returned.
    """

    sides = _numeric_side_lengths(side_lengths)
    bound = _positive_numeric(maximum_length, "maximum_length")
    component_bounds = tuple(int(sp.floor(bound / side)) for side in sides)
    translations: list[LatticeTranslation] = []
    ranges = tuple(range(-limit, limit + 1) for limit in component_bounds)
    for index in product(*ranges):
        if not any(index):
            continue
        if unique_up_to_sign:
            first_nonzero = next(component for component in index if component)
            if first_nonzero < 0:
                continue
        squared_length = translation_squared_length(index, sides)
        if _strictly_less(
            squared_length,
            bound**2,
            "translation cutoff comparison",
        ):
            translations.append(
                LatticeTranslation(
                    index=tuple(index),
                    vector=translation_vector(index, sides),
                    squared_length=squared_length,
                )
            )
    translations.sort(
        key=lambda item: (float(sp.N(item.squared_length, 30)), item.index)
    )
    return tuple(translations)


def has_nondegenerate_matched_circles(
    side_lengths: Sequence[Any], observation_radius: Any
) -> bool:
    """Return whether any torus image sphere intersects in a true circle.

    For a rectangular torus the shortest image separation is ``min_i L_i``.
    Radius-``chi`` image spheres have a nondegenerate intersection exactly
    when that separation is strictly less than ``2*chi``.  Equality is a
    tangent point, not a circle.
    """

    sides = _numeric_side_lengths(side_lengths)
    radius = _positive_numeric(observation_radius, "observation_radius")
    return _strictly_less(
        shortest_translation_length(sides),
        2 * radius,
        "matched-circle existence comparison",
    )


def matched_circle_angular_radius(
    translation_length: Any, observation_radius: Any
) -> sp.Expr:
    """Return ``acos(d/(2*chi))`` for two intersecting radius-``chi`` spheres."""

    distance = _positive_numeric(translation_length, "translation_length")
    radius = _positive_numeric(observation_radius, "observation_radius")
    if not _strictly_less(
        distance,
        2 * radius,
        "matched-circle intersection comparison",
    ):
        raise ValueError(
            "translation_length must be strictly less than twice the "
            "observation_radius"
        )
    return sp.acos(sp.simplify(distance / (2 * radius)))


def matched_circle_geometries(
    side_lengths: Sequence[Any],
    observation_radius: Any,
    *,
    unique_pairs: bool = True,
) -> tuple[MatchedCircleGeometry, ...]:
    """Return all nondegenerate image-sphere circles for the radius.

    With the default ``unique_pairs=True``, translations related by sign are
    represented once because they define the two members of one circle pair.
    """

    radius = _positive_numeric(observation_radius, "observation_radius")
    translations = enumerate_translations_below(
        side_lengths,
        2 * radius,
        unique_up_to_sign=unique_pairs,
    )
    return tuple(
        MatchedCircleGeometry(
            translation=translation,
            angular_radius=matched_circle_angular_radius(
                translation.length, radius
            ),
        )
        for translation in translations
    )


def rectangular_epstein_partial_sum(
    side_lengths: Sequence[Any],
    exponent: Any,
    max_index: Any,
    *,
    precision: int = 50,
) -> mpf:
    r"""Return a cubic-cutoff sum of a convergent rectangular Epstein zeta.

    The implemented quantity is
    ``sum_{n in [-N,N]^d, n != 0} (sum_i (n_i L_i)^2)^(-s)``.
    The domain ``s > d/2`` makes the infinite positive-term series absolutely
    convergent.  This helper returns only the stated finite partial sum; it
    performs no analytic continuation, subtraction, Casimir normalization,
    vacuum-energy map, or error extrapolation.
    """

    sides = _numeric_side_lengths(side_lengths)
    s = _positive_numeric(exponent, "exponent")
    dimension = len(sides)
    if not _strictly_less(
        sp.Rational(dimension, 2), s, "Epstein convergence comparison"
    ):
        raise ValueError("exponent must be greater than dimension/2")
    cutoff = _nonnegative_integer(max_index, "max_index")
    if cutoff == 0:
        raise ValueError("max_index must be positive")
    if not isinstance(precision, int) or precision < 15:
        raise ValueError("precision must be an integer at least 15")

    with mp.workdps(precision):
        numeric_sides = tuple(mpf(str(sp.N(side, precision))) for side in sides)
        numeric_exponent = mpf(str(sp.N(s, precision)))
        terms: list[mpf] = []
        for index in product(range(-cutoff, cutoff + 1), repeat=dimension):
            if not any(index):
                continue
            quadratic_form = mp.fsum(
                (integer * side) ** 2
                for integer, side in zip(index, numeric_sides, strict=True)
            )
            terms.append(mp.power(quadratic_form, -numeric_exponent))
        return +mp.fsum(terms)


def rectangular_epstein_refinement(
    side_lengths: Sequence[Any],
    exponent: Any,
    cutoffs: Iterable[Any],
    *,
    precision: int = 50,
) -> tuple[EpsteinRefinementPoint, ...]:
    """Evaluate increasing cubic partial sums without inventing an extrapolation."""

    parsed = tuple(_nonnegative_integer(value, "cutoff") for value in cutoffs)
    if not parsed or any(value == 0 for value in parsed):
        raise ValueError("cutoffs must contain positive integers")
    if any(right <= left for left, right in zip(parsed, parsed[1:])):
        raise ValueError("cutoffs must be strictly increasing")
    points: list[EpsteinRefinementPoint] = []
    prior: mpf | None = None
    for cutoff in parsed:
        value = rectangular_epstein_partial_sum(
            side_lengths,
            exponent,
            cutoff,
            precision=precision,
        )
        increment = None if prior is None else value - prior
        points.append(
            EpsteinRefinementPoint(
                max_index=cutoff,
                partial_sum=value,
                increment=increment,
            )
        )
        prior = value
    return tuple(points)
