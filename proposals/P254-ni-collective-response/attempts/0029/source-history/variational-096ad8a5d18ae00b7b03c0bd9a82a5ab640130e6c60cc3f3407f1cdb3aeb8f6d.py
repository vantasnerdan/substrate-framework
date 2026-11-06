"""Pure symbolic variational utilities with explicit functional scope."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp


@dataclass(frozen=True)
class FiniteFunctionalInfimumLedger:
    """Exact pointwise excess ledger relative to supplied component infima.

    The caller supplies the true finite component infima on one common
    nonempty admissible set.  The ledger does not infer those infima or claim
    that they are attained.  For any common configuration, every component
    excess is nonnegative and their sum is the excess of the functional sum
    above the sum of the separate infima.
    """

    component_values: tuple[sp.Expr, ...]
    component_infima: tuple[sp.Expr, ...]
    component_excesses: tuple[sp.Expr, ...]
    summed_value: sp.Expr
    separate_infimum_sum: sp.Expr
    total_excess: sp.Expr

    @property
    def identity_residual(self) -> sp.Expr:
        """Return the exact sum-of-excesses identity residual."""

        return sp.simplify(
            self.total_excess - sum(self.component_excesses, sp.S.Zero)
        )


@dataclass(frozen=True)
class FiniteFunctionalInteractionLedger:
    """Exact mixed difference of four supplied variational infima.

    For a base functional ``A`` and additions ``P`` and ``Q`` on one common
    nonempty domain, callers supply the true finite infima of ``A``, ``A+P``,
    ``A+Q``, and ``A+P+Q``.  The interaction is the joint increment minus the
    two separate increments.  Its sign is not inferred: optimization can make
    it positive, negative, or zero even when all three functionals are
    nonnegative continuous coercive quadratics.

    The ledger does not infer any infimum, minimizer, field decomposition, or
    physical sector interpretation.
    """

    base_infimum: sp.Expr
    first_augmented_infimum: sp.Expr
    second_augmented_infimum: sp.Expr
    joint_augmented_infimum: sp.Expr
    first_increment: sp.Expr
    second_increment: sp.Expr
    joint_increment: sp.Expr
    interaction: sp.Expr

    @property
    def identity_residual(self) -> sp.Expr:
        """Return the mixed-increment identity residual."""

        return sp.simplify(
            self.interaction
            - (self.joint_increment - self.first_increment - self.second_increment)
        )


def _nonempty_expressions(values: Iterable[Any], name: str) -> tuple[sp.Expr, ...]:
    result = tuple(sp.sympify(value) for value in values)
    if not result:
        raise ValueError(f"{name} must be nonempty")
    return result


def finite_functional_infimum_ledger(
    component_values: Iterable[Any],
    component_infima: Iterable[Any],
) -> FiniteFunctionalInfimumLedger:
    """Return an exact common-configuration functional excess ledger.

    ``component_values`` must all be evaluated on the same admissible
    configuration. ``component_infima`` are separately supplied lower bounds
    that are assumed to be the corresponding exact infima.  Provably negative
    component excesses are rejected because they contradict that premise.
    """

    values = _nonempty_expressions(component_values, "component_values")
    infima = _nonempty_expressions(component_infima, "component_infima")
    if len(values) != len(infima):
        raise ValueError("component_values and component_infima must have equal length")
    excesses = tuple(sp.simplify(value - infimum) for value, infimum in zip(values, infima))
    if any(excess.is_negative is True for excess in excesses):
        raise ValueError("component value cannot lie below its supplied infimum")
    summed_value = sp.simplify(sum(values, sp.S.Zero))
    separate_sum = sp.simplify(sum(infima, sp.S.Zero))
    return FiniteFunctionalInfimumLedger(
        component_values=values,
        component_infima=infima,
        component_excesses=excesses,
        summed_value=summed_value,
        separate_infimum_sum=separate_sum,
        total_excess=sp.simplify(summed_value - separate_sum),
    )


def finite_functional_interaction_ledger(
    base_infimum: Any,
    first_augmented_infimum: Any,
    second_augmented_infimum: Any,
    joint_augmented_infimum: Any,
) -> FiniteFunctionalInteractionLedger:
    """Return the exact mixed interaction of four supplied finite infima.

    Writing the inputs as ``m_A``, ``m_AP``, ``m_AQ``, and ``m_APQ``, the
    result is ``m_APQ + m_A - m_AP - m_AQ``.  Additive constants in any one of
    ``A``, ``P``, or ``Q`` cancel when callers transform all four infima
    consistently.
    """

    infima = tuple(
        sp.sympify(value)
        for value in (
            base_infimum,
            first_augmented_infimum,
            second_augmented_infimum,
            joint_augmented_infimum,
        )
    )
    if any(
        value.is_finite is False
        or value.is_real is False
        or value in {sp.oo, -sp.oo, sp.zoo, sp.nan}
        for value in infima
    ):
        raise ValueError("all supplied infima must be finite real expressions")
    base, first_augmented, second_augmented, joint_augmented = infima
    first_increment = sp.simplify(first_augmented - base)
    second_increment = sp.simplify(second_augmented - base)
    joint_increment = sp.simplify(joint_augmented - base)
    interaction = sp.simplify(
        joint_augmented + base - first_augmented - second_augmented
    )
    return FiniteFunctionalInteractionLedger(
        base_infimum=base,
        first_augmented_infimum=first_augmented,
        second_augmented_infimum=second_augmented,
        joint_augmented_infimum=joint_augmented,
        first_increment=first_increment,
        second_increment=second_increment,
        joint_increment=joint_increment,
        interaction=interaction,
    )


def euler_lagrange_expression(
    lagrangian: Any,
    coordinate: sp.Expr,
    parameter: sp.Symbol,
) -> sp.Expr:
    """Return ``d/dt(dL/dqdot) - dL/dq`` for one coordinate."""

    expression = sp.sympify(lagrangian)
    velocity = sp.diff(coordinate, parameter)
    return sp.simplify(
        sp.diff(sp.diff(expression, velocity), parameter)
        - sp.diff(expression, coordinate)
    )


def solve_euler_lagrange_acceleration(
    lagrangian: Any,
    coordinate: sp.Expr,
    parameter: sp.Symbol,
) -> sp.Expr:
    """Solve a one-coordinate Euler-Lagrange equation for its acceleration."""

    acceleration = sp.diff(coordinate, parameter, 2)
    placeholder = sp.Dummy("acceleration")
    residual = euler_lagrange_expression(lagrangian, coordinate, parameter)
    solutions = sp.solve(residual.subs(acceleration, placeholder), placeholder)
    if len(solutions) != 1:
        raise ValueError(
            "Euler-Lagrange equation must have one explicit acceleration solution"
        )
    return sp.simplify(solutions[0])
