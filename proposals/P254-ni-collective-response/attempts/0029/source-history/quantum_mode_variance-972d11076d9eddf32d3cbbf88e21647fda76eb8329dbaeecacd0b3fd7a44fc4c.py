"""Conditional scalar ground-state mode variances and continuum moments.

Every function in this module assumes a separately declared canonical scalar
quadratic action and product ground state. The accepted classical
one-dimensional medium does not supply that quantization or a
three-dimensional lift. The radial cutoff, branch count, stiffness, and state
are inputs; the returned identities do not derive a microscopic cutoff,
material variance, participating-mode set, growth law, channel, or rate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp

from .mode_counting import isotropic_gapped_angular_frequency


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    if expression.is_positive is False:
        raise ValueError(f"{name} must not be known nonpositive")
    return expression


def _nonnegative(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    if expression.is_nonnegative is False:
        raise ValueError(f"{name} must not be known negative")
    return expression


def _positive_integer(value: Any, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be a positive integer")
    return int(expression)


def scalar_mode_ground_state_variance(
    action_scale: Any,
    quantization_volume: Any,
    stiffness: Any,
    signal_speed: Any,
    angular_frequency: Any,
) -> sp.Expr:
    r"""Return ``hbar*c**2/(2*V*kappa*omega)`` exactly.

    This is the coordinate variance of one oscillator only under the declared
    normalization in which its effective mass is ``V*kappa/c**2`` and it is in
    its ground state. The function does not derive that normalization or state
    from a classical dispersion relation.
    """

    hbar = _positive(action_scale, "action_scale")
    volume = _positive(quantization_volume, "quantization_volume")
    kappa = _positive(stiffness, "stiffness")
    speed = _positive(signal_speed, "signal_speed")
    frequency = _positive(angular_frequency, "angular_frequency")
    return sp.simplify(hbar * speed**2 / (2 * volume * kappa * frequency))


def gapped_vacuum_kernel(argument: Any) -> sp.Expr:
    r"""Return ``J(X)=(X*sqrt(1+X**2)-asinh(X))/2`` for ``X>=0``.

    Its derivative is ``X**2/sqrt(1+X**2)`` and therefore it equals the
    corresponding integral from zero to ``X``.
    """

    value = _nonnegative(argument, "argument")
    return sp.simplify(
        (value * sp.sqrt(1 + value**2) - sp.asinh(value)) / 2
    )


def scalar_continuum_vacuum_variance_3d(
    action_scale: Any,
    stiffness: Any,
    signal_speed: Any,
    gap_frequency: Any,
    radial_cutoff: Any,
    *,
    branches: Any = 1,
) -> sp.Expr:
    r"""Return the cutoff scalar ground-state variance in three dimensions.

    With one branch, dispersion ``sqrt(omega_0**2+c**2*k**2)``, per-mode
    variance ``hbar*c**2/(2*V*kappa*omega_k)``, and continuum measure
    ``V*4*pi*k**2*dk/(2*pi)**3``, the exact positive-gap result is

    ``hbar/(8*pi**2*kappa) *``
    ``(K*sqrt(omega_0**2+c**2*K**2)``
    `` - omega_0**2*asinh(c*K/omega_0)/c)``.

    Independent identical branches multiply this result. An exactly zero gap
    uses its continuous limit ``hbar*c*K**2/(8*pi**2*kappa)``. A symbolic gap
    must therefore be provably positive or exactly zero.
    """

    hbar = _positive(action_scale, "action_scale")
    kappa = _positive(stiffness, "stiffness")
    speed = _positive(signal_speed, "signal_speed")
    gap = _nonnegative(gap_frequency, "gap_frequency")
    cutoff = _positive(radial_cutoff, "radial_cutoff")
    branch_count = _positive_integer(branches, "branches")

    if gap.is_zero is True:
        return sp.simplify(
            branch_count * hbar * speed * cutoff**2 / (8 * sp.pi**2 * kappa)
        )
    if gap.is_positive is not True:
        raise ValueError("gap_frequency must be exactly zero or provably positive")

    upper_frequency = isotropic_gapped_angular_frequency(cutoff, speed, gap)
    return sp.simplify(
        branch_count
        * hbar
        / (8 * sp.pi**2 * kappa)
        * (
            cutoff * upper_frequency
            - gap**2 * sp.asinh(speed * cutoff / gap) / speed
        )
    )


@dataclass(frozen=True)
class ModeVarianceLedger:
    """Exact total and arithmetic mean for one fixed finite mode set."""

    variances: tuple[sp.Expr, ...]
    count: sp.Integer
    total: sp.Expr
    arithmetic_mean: sp.Expr

    @property
    def factorization_residual(self) -> sp.Expr:
        """Return ``count*arithmetic_mean-total`` exactly."""

        return sp.simplify(self.count * self.arithmetic_mean - self.total)


def mode_variance_ledger(variances: Iterable[Any]) -> ModeVarianceLedger:
    r"""Return the total and arithmetic mean of a fixed nonempty mode set.

    The identity ``M*mean=S`` holds for the returned fixed sequence. It does
    not assert that adding, removing, or retuning modes leaves ``S`` fixed.
    """

    values = tuple(_nonnegative(value, "variance") for value in variances)
    if not values:
        raise ValueError("variances must be a nonempty finite iterable")
    count = sp.Integer(len(values))
    total = sp.simplify(sum(values, sp.Integer(0)))
    return ModeVarianceLedger(
        variances=values,
        count=count,
        total=total,
        arithmetic_mean=sp.simplify(total / count),
    )
