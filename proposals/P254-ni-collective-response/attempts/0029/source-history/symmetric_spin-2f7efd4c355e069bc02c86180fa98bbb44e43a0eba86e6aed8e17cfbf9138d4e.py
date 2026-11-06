"""Exact collective-spin algebra on a declared symmetric two-state sector.

The helpers in this module concern normalized vectors in the permutation-
symmetric subspace of a finite tensor product.  They do not derive physical
two-level constituents, prepare a Dicke state, define an interaction
Hamiltonian, or turn a squared ladder coefficient into a transition rate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Sequence

import sympy as sp


LadderDirection = Literal["raise", "lower"]


def _nonnegative_integer(value: Any, *, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(expression)


def _positive_integer(value: Any, *, name: str) -> int:
    integer = _nonnegative_integer(value, name=name)
    if integer == 0:
        raise ValueError(f"{name} must be a positive integer")
    return integer


def _magnitude_squared(value: Any) -> sp.Expr:
    expression = sp.sympify(value)
    return sp.simplify(sp.conjugate(expression) * expression)


def _real_expression(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


@dataclass(frozen=True)
class SymmetricSpinRung:
    """Exact data for one normalized symmetric excitation rung."""

    particle_count: int
    excitation_count: int
    total_spin: sp.Rational
    magnetic_number: sp.Rational
    raising_coefficient: sp.Expr
    lowering_coefficient: sp.Expr
    raising_coefficient_squared: sp.Expr
    lowering_coefficient_squared: sp.Expr


@dataclass(frozen=True)
class GroundCouplingLedger:
    """Ground-to-one-excitation data for declared complex site couplings.

    ``symmetric_amplitude`` is the projection on the normalized equal-weight
    one-excitation vector.  ``total_norm_squared`` is the norm of the complete
    one-excitation image, and ``dark_norm_squared`` is the norm left in its
    orthogonal complement.  These are vector-space quantities, not rates.
    """

    particle_count: int
    couplings: tuple[sp.Expr, ...]
    symmetric_amplitude: sp.Expr
    symmetric_norm_squared: sp.Expr
    total_norm_squared: sp.Expr
    dark_norm_squared: sp.Expr


def symmetric_spin_rung(
    particle_count: Any,
    excitation_count: Any,
    *,
    operator_scale: Any = 1,
) -> SymmetricSpinRung:
    """Return exact raising and lowering coefficients on ``|D_N^k>``.

    ``|D_N^k>`` is the normalized equal superposition of all computational
    basis vectors with exactly ``k`` excited factors.  With
    ``J_+=operator_scale*sum_i sigma_i^+``, the raising coefficient is
    ``operator_scale*sqrt((N-k)*(k+1))``.  The lowering coefficient is
    ``operator_scale*sqrt(k*(N-k+1))``.  A top or bottom edge is represented
    by an exact zero rather than by an out-of-domain target vector.
    """

    count = _positive_integer(particle_count, name="particle_count")
    excitation = _nonnegative_integer(
        excitation_count,
        name="excitation_count",
    )
    if excitation > count:
        raise ValueError("excitation_count must not exceed particle_count")
    scale = _real_expression(operator_scale, name="operator_scale")
    raising_factor = sp.Integer((count - excitation) * (excitation + 1))
    lowering_factor = sp.Integer(excitation * (count - excitation + 1))
    raising = sp.simplify(scale * sp.sqrt(raising_factor))
    lowering = sp.simplify(scale * sp.sqrt(lowering_factor))
    return SymmetricSpinRung(
        particle_count=count,
        excitation_count=excitation,
        total_spin=sp.Rational(count, 2),
        magnetic_number=sp.Rational(2 * excitation - count, 2),
        raising_coefficient=raising,
        lowering_coefficient=lowering,
        raising_coefficient_squared=sp.simplify(
            _magnitude_squared(scale) * raising_factor
        ),
        lowering_coefficient_squared=sp.simplify(
            _magnitude_squared(scale) * lowering_factor
        ),
    )


def symmetric_spin_ladder_coefficient(
    particle_count: Any,
    excitation_count: Any,
    *,
    direction: LadderDirection = "raise",
    operator_scale: Any = 1,
) -> sp.Expr:
    """Return one exact normalized symmetric-ladder coefficient."""

    rung = symmetric_spin_rung(
        particle_count,
        excitation_count,
        operator_scale=operator_scale,
    )
    if direction == "raise":
        return rung.raising_coefficient
    if direction == "lower":
        return rung.lowering_coefficient
    raise ValueError("direction must be 'raise' or 'lower'")


def ground_coupling_ledger(
    couplings: Sequence[Any],
    *,
    operator_scale: Any = 1,
) -> GroundCouplingLedger:
    """Project a weighted ground-state raise into bright and dark sectors.

    For ``A_+=operator_scale*sum_i g_i sigma_i^+``, the normalized symmetric
    amplitude is ``operator_scale*sum_i(g_i)/sqrt(N)``.  Its equality with
    ``operator_scale*g*sqrt(N)`` requires all site couplings to equal the same
    ``g`` with the same phase.  The complete image norm instead depends on
    ``sum_i |g_i|^2``.
    """

    try:
        exact_couplings = tuple(sp.sympify(value) for value in couplings)
    except TypeError as error:
        raise ValueError("couplings must be a nonempty finite sequence") from error
    if not exact_couplings:
        raise ValueError("couplings must be a nonempty finite sequence")
    count = len(exact_couplings)
    scale = sp.sympify(operator_scale)
    symmetric_amplitude = sp.simplify(
        scale * sum(exact_couplings, sp.Integer(0)) / sp.sqrt(count)
    )
    symmetric_norm_squared = _magnitude_squared(symmetric_amplitude)
    total_norm_squared = sp.simplify(
        _magnitude_squared(scale)
        * sum((_magnitude_squared(value) for value in exact_couplings), sp.Integer(0))
    )
    dark_norm_squared = sp.simplify(total_norm_squared - symmetric_norm_squared)
    return GroundCouplingLedger(
        particle_count=count,
        couplings=exact_couplings,
        symmetric_amplitude=symmetric_amplitude,
        symmetric_norm_squared=symmetric_norm_squared,
        total_norm_squared=total_norm_squared,
        dark_norm_squared=dark_norm_squared,
    )
