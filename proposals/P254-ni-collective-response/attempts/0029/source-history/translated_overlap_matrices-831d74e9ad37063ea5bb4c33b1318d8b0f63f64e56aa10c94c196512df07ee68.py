"""Exact limit ledgers for identical translated overlap matrices.

Let one normalized mode shape and one bounded multiplier profile be translated
to finitely many mutually separating centers.  If the multiplier vanishes at
spatial infinity, every mixed translated overlap vanishes while each matched
self-overlap remains the same.  A phase-weighted overlap matrix therefore
approaches a diagonal matrix whose entries have one common magnitude.  This
module records that limit and the finite-residual singular-value bound; it does
not identify the modes, profiles, or matrix singular values with physical
generations or masses.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp


def _positive_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and (
        result.is_real is not True or result.is_positive is not True
    ):
        raise ValueError(f"{name} must be positive and real")
    return result


def _nonnegative_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and (
        result.is_real is not True or result.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be nonnegative and real")
    return result


def _real_phases(phases: Iterable[Any]) -> tuple[sp.Expr, ...]:
    values = tuple(sp.sympify(phase) for phase in phases)
    if not values:
        raise ValueError("at least one phase is required")
    if any(value.is_number and value.is_real is not True for value in values):
        raise ValueError("phases must be real")
    return values


def phase_weighted_self_overlap_limit(
    self_overlap: Any,
    phases: Iterable[Any],
) -> sp.ImmutableMatrix:
    """Return ``alpha*diag(exp(i*theta_a))`` for ``alpha>0``.

    This is the separated-center limit of the matrix

    ``Y_ab=sum_c exp(i*theta_c) integral conjugate(psi_a)*psi_b*Phi_c dx``

    when all ``psi_a`` and ``Phi_c`` are identical translates, the distinct
    mixed overlaps vanish, and the matched self-overlap is ``alpha``.  The
    function records the exact limiting matrix, not the convergence premise.
    """

    alpha = _positive_real(self_overlap, "self_overlap")
    angles = _real_phases(phases)
    return sp.ImmutableMatrix.diag(
        *(alpha * sp.exp(sp.I * angle) for angle in angles)
    )


@dataclass(frozen=True)
class SingularValueClusterBound:
    """Weyl-type cluster bound around a common positive singular value."""

    count: int
    self_overlap: sp.Expr
    residual_operator_norm: sp.Expr
    singular_value_lower_bound: sp.Expr
    singular_value_upper_bound: sp.Expr
    condition_number_upper_bound: sp.Expr


def singular_value_cluster_bound(
    self_overlap: Any,
    residual_operator_norm: Any,
    *,
    count: int,
) -> SingularValueClusterBound:
    """Bound every singular value of ``alpha*D+E`` for unitary diagonal ``D``.

    If ``alpha>0`` and ``||E||_2=epsilon<alpha``, singular-value perturbation
    gives ``alpha-epsilon <= sigma_j <= alpha+epsilon`` for every ``j`` and
    hence ``condition_2 <= (alpha+epsilon)/(alpha-epsilon)``.  Requiring the
    strict residual gate keeps the condition bound finite and prevents a
    near-zero lower singular value from being hidden.
    """

    if isinstance(count, bool) or not isinstance(count, int):
        raise TypeError("count must be an integer")
    if count < 1:
        raise ValueError("count must be positive")
    alpha = _positive_real(self_overlap, "self_overlap")
    epsilon = _nonnegative_real(
        residual_operator_norm, "residual_operator_norm"
    )
    if alpha.is_number and epsilon.is_number and not bool(epsilon < alpha):
        raise ValueError("residual_operator_norm must be smaller than self_overlap")
    lower = sp.simplify(alpha - epsilon)
    upper = sp.simplify(alpha + epsilon)
    return SingularValueClusterBound(
        count=count,
        self_overlap=alpha,
        residual_operator_norm=epsilon,
        singular_value_lower_bound=lower,
        singular_value_upper_bound=upper,
        condition_number_upper_bound=sp.simplify(upper / lower),
    )
