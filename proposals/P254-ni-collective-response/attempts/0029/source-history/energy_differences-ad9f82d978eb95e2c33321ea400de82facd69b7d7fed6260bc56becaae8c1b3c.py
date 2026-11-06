"""Pure signed linear-energy differences with explicit normalization.

The helpers in this module know no particle, state, reaction, or empirical
scale.  They only implement a declared linear combination and its monotone
rectangular interval image.  Physical interpretations require separately
governed inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


def _finite_real(name: str, value: float) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_multiplicity(multiplicity: int) -> int:
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int):
        raise TypeError("multiplicity must be a positive integer")
    if multiplicity <= 0:
        raise ValueError("multiplicity must be a positive integer")
    return multiplicity


def _positive_real(name: str, value: float) -> float:
    result = _finite_real(name, value)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


@dataclass(frozen=True)
class LinearDifferenceInterval:
    """Closed interval obtained from independent coefficient intervals."""

    lower: float
    upper: float

    def __post_init__(self) -> None:
        lower = _finite_real("lower", self.lower)
        upper = _finite_real("upper", self.upper)
        if lower > upper:
            raise ValueError("lower must not exceed upper")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def width(self) -> float:
        """Return the interval width."""

        return self.upper - self.lower

    def contains(self, value: float) -> bool:
        """Return whether a finite real value lies in the closed interval."""

        candidate = _finite_real("value", value)
        return self.lower <= candidate <= self.upper


def normalized_linear_difference(
    initial_coefficient: float,
    final_coefficient: float,
    *,
    multiplicity: int,
) -> float:
    """Return ``multiplicity*initial_coefficient-final_coefficient``."""

    count = _positive_multiplicity(multiplicity)
    initial = _finite_real("initial_coefficient", initial_coefficient)
    final = _finite_real("final_coefficient", final_coefficient)
    return count * initial - final


def linear_difference_coefficient(
    initial_coefficient: float,
    final_coefficient: float,
    *,
    multiplicity: int,
    normalization: float,
) -> float:
    """Scale a normalized signed difference by a declared positive factor."""

    factor = _positive_real("normalization", normalization)
    return factor * normalized_linear_difference(
        initial_coefficient,
        final_coefficient,
        multiplicity=multiplicity,
    )


def linear_energy_difference(
    initial_coefficient: float,
    final_coefficient: float,
    *,
    multiplicity: int,
    normalization: float,
    energy_scale: float,
) -> float:
    """Return the declared energy-scale multiple of a linear coefficient."""

    scale = _positive_real("energy_scale", energy_scale)
    return scale * linear_difference_coefficient(
        initial_coefficient,
        final_coefficient,
        multiplicity=multiplicity,
        normalization=normalization,
    )


def linear_difference_interval(
    initial_interval: tuple[float, float],
    final_interval: tuple[float, float],
    *,
    multiplicity: int,
    normalization: float,
) -> LinearDifferenceInterval:
    """Map independent closed intervals through a positive linear difference.

    For ``initial in [li, ui]`` and ``final in [lf, uf]``, positive
    normalization and multiplicity make the image
    ``[normalization*(multiplicity*li-uf),
    normalization*(multiplicity*ui-lf)]``.  The result is only as rigorous as
    the supplied input intervals.
    """

    count = _positive_multiplicity(multiplicity)
    factor = _positive_real("normalization", normalization)
    initial = LinearDifferenceInterval(*initial_interval)
    final = LinearDifferenceInterval(*final_interval)
    return LinearDifferenceInterval(
        factor * (count * initial.lower - final.upper),
        factor * (count * initial.upper - final.lower),
    )
