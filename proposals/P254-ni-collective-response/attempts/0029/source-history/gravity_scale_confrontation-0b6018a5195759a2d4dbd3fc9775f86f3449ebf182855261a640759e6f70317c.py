"""Exact ledgers for confronting induced gravity with a length hierarchy.

The APIs in this module compose accepted conditional relations while keeping
their empirical rows and free coordinates visible.  They do not derive a
Newton constant, cutoff, field-count interval, beta function, coupling,
conversion factor, observed length, or physical identification of either
length.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .induced_gravity import (
    cutoff_length_from_pure_induced_newton,
    induced_inverse_newton_ledger,
)
from .scale_constraints import LogConstraintDiagnostics, diagnose_log_constraints


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _provenance(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value


@dataclass(frozen=True)
class PureGravityCutoffIntervalLedger:
    """Image of a supplied positive coefficient interval on the pure branch."""

    newton_constant: sp.Expr
    signal_speed: sp.Expr
    action_scale: sp.Expr
    coefficient_lower: sp.Expr
    coefficient_upper: sp.Expr
    unit_coefficient_cutoff: sp.Expr
    cutoff_lower: sp.Expr
    cutoff_upper: sp.Expr
    coefficient_interval_width: sp.Expr
    cutoff_interval_width: sp.Expr


@dataclass(frozen=True)
class JointGravityTransmutationLogLedger:
    """Two provenance-labelled log rows on ``(log a, log s, 1/g^2)``."""

    reduced_newton_ratio: sp.Expr
    reduced_long_length_ratio: sp.Expr
    conversion_ratio: sp.Expr
    beta_coefficient: sp.Expr
    exponent_coefficient: sp.Expr
    system: LogConstraintDiagnostics


@dataclass(frozen=True)
class FixedCoefficientJointSolution:
    """Unique solve obtained only after supplying a coefficient ratio."""

    joint: JointGravityTransmutationLogLedger
    coefficient_ratio: sp.Expr
    system: LogConstraintDiagnostics
    log_solution: sp.ImmutableMatrix
    cutoff_ratio: sp.Expr
    inverse_coupling_squared: sp.Expr
    inferred_coupling_squared: sp.Expr | None
    positive_coupling_admissible: bool | None
    residuals: tuple[sp.Expr, ...]


def pure_gravity_cutoff_interval(
    newton_constant: Any,
    coefficient_lower: Any,
    coefficient_upper: Any,
    signal_speed: Any,
    action_scale: Any,
) -> PureGravityCutoffIntervalLedger:
    """Map a supplied closed coefficient interval through ``a=sqrt(s hbar G/c^3)``.

    The result is conditional on a zero additive inverse-coupling baseline.
    The coefficient bounds must be exact, positive, and provably ordered.
    """

    newton = _exact_positive(newton_constant, "newton_constant")
    lower = _exact_positive(coefficient_lower, "coefficient_lower")
    upper = _exact_positive(coefficient_upper, "coefficient_upper")
    speed = _exact_positive(signal_speed, "signal_speed")
    action = _exact_positive(action_scale, "action_scale")
    width = sp.simplify(upper - lower)
    if width.is_nonnegative is not True:
        raise ValueError("coefficient bounds must be provably ordered")

    unit = cutoff_length_from_pure_induced_newton(newton, 1, speed, action)
    cutoff_lower = cutoff_length_from_pure_induced_newton(
        newton, lower, speed, action
    )
    cutoff_upper = cutoff_length_from_pure_induced_newton(
        newton, upper, speed, action
    )
    return PureGravityCutoffIntervalLedger(
        newton_constant=newton,
        signal_speed=speed,
        action_scale=action,
        coefficient_lower=lower,
        coefficient_upper=upper,
        unit_coefficient_cutoff=unit,
        cutoff_lower=cutoff_lower,
        cutoff_upper=cutoff_upper,
        coefficient_interval_width=width,
        cutoff_interval_width=sp.simplify(cutoff_upper - cutoff_lower),
    )


def induced_coefficient_for_target_cutoff(
    target_cutoff: Any,
    newton_constant: Any,
    signal_speed: Any,
    action_scale: Any,
) -> sp.Expr:
    """Return the pure-branch coefficient required by a supplied target cutoff."""

    target = _exact_positive(target_cutoff, "target_cutoff")
    newton = _exact_positive(newton_constant, "newton_constant")
    speed = _exact_positive(signal_speed, "signal_speed")
    action = _exact_positive(action_scale, "action_scale")
    return sp.simplify(target**2 * speed**3 / (action * newton))


def inverse_newton_baseline_for_target(
    target_newton_constant: Any,
    cutoff_length: Any,
    induced_coefficient: Any,
    signal_speed: Any,
    action_scale: Any,
) -> sp.Expr:
    """Return the additive baseline needed to realize a supplied total ``G``."""

    target_newton = _exact_positive(
        target_newton_constant, "target_newton_constant"
    )
    cutoff = _exact_positive(cutoff_length, "cutoff_length")
    coefficient = sp.sympify(induced_coefficient)
    if coefficient.has(sp.Float):
        raise ValueError("induced_coefficient must be exact rather than floating")
    if coefficient.is_real is not True or coefficient.is_zero is not False:
        raise ValueError("induced_coefficient must be provably nonzero and real")
    speed = _exact_positive(signal_speed, "signal_speed")
    action = _exact_positive(action_scale, "action_scale")
    shift = induced_inverse_newton_ledger(
        cutoff,
        coefficient,
        speed,
        action,
    ).induced_inverse_newton
    return sp.simplify(1 / target_newton - shift)


def joint_gravity_transmutation_log_ledger(
    reduced_newton_ratio: Any,
    reduced_long_length_ratio: Any,
    beta_coefficient: Any,
    *,
    conversion_ratio: Any = 1,
    gravity_provenance: str,
    length_provenance: str,
) -> JointGravityTransmutationLogLedger:
    r"""Compose two exact rows on ``u=log(a/a0), v=log(s/s0), y=1/g^2``.

    Compatible pure-gravity references give ``2u-v=log(G/G0)``.  A separately
    supplied relation ``L/a=C*exp(K*y)``, with ``K=8*pi^2/b0``, gives
    ``u+K*y=log((L/a0)/C)``.  The two rows do not identify all three
    coordinates; row provenance establishes neither empirical independence nor
    a physical dictionary between the sectors.
    """

    newton_ratio = _exact_positive(
        reduced_newton_ratio, "reduced_newton_ratio"
    )
    length_ratio = _exact_positive(
        reduced_long_length_ratio, "reduced_long_length_ratio"
    )
    b0 = _exact_positive(beta_coefficient, "beta_coefficient")
    conversion = _exact_positive(conversion_ratio, "conversion_ratio")
    gravity_label = _provenance(gravity_provenance, "gravity_provenance")
    length_label = _provenance(length_provenance, "length_provenance")
    exponent = sp.simplify(8 * sp.pi**2 / b0)
    system = diagnose_log_constraints(
        [[2, -1, 0], [1, 0, exponent]],
        [
            sp.log(newton_ratio),
            sp.simplify(sp.log(length_ratio) - sp.log(conversion)),
        ],
        provenance=[gravity_label, length_label],
    )
    return JointGravityTransmutationLogLedger(
        reduced_newton_ratio=newton_ratio,
        reduced_long_length_ratio=length_ratio,
        conversion_ratio=conversion,
        beta_coefficient=b0,
        exponent_coefficient=exponent,
        system=system,
    )


def solve_joint_with_fixed_coefficient_ratio(
    reduced_newton_ratio: Any,
    reduced_long_length_ratio: Any,
    coefficient_ratio: Any,
    beta_coefficient: Any,
    *,
    conversion_ratio: Any = 1,
    gravity_provenance: str,
    length_provenance: str,
    coefficient_provenance: str,
) -> FixedCoefficientJointSolution:
    """Solve the joint rows after a third input fixes ``s/s0``.

    A positive inferred coupling exists only when the solved coordinate
    ``y=1/g^2`` is positive.  The returned admissibility flag is ``None`` when
    exact symbolic assumptions do not decide that sign; no empirical value is
    selected or fitted by this routine.
    """

    joint = joint_gravity_transmutation_log_ledger(
        reduced_newton_ratio,
        reduced_long_length_ratio,
        beta_coefficient,
        conversion_ratio=conversion_ratio,
        gravity_provenance=gravity_provenance,
        length_provenance=length_provenance,
    )
    coefficient = _exact_positive(coefficient_ratio, "coefficient_ratio")
    coefficient_label = _provenance(
        coefficient_provenance, "coefficient_provenance"
    )
    design = sp.Matrix.vstack(sp.Matrix(joint.system.design), sp.Matrix([[0, 1, 0]]))
    rhs = sp.Matrix.vstack(
        sp.Matrix(joint.system.rhs), sp.Matrix([sp.log(coefficient)])
    )
    system = diagnose_log_constraints(
        design,
        rhs,
        provenance=(*joint.system.provenance, coefficient_label),
    )
    solution = sp.ImmutableMatrix(design.inv() * rhs)
    residual_vector = design * sp.Matrix(solution) - rhs
    inverse_coupling = sp.simplify(solution[2])
    inferred = None
    if inverse_coupling.is_zero is not True:
        inferred = sp.simplify(1 / inverse_coupling)
    return FixedCoefficientJointSolution(
        joint=joint,
        coefficient_ratio=coefficient,
        system=system,
        log_solution=solution,
        cutoff_ratio=sp.simplify(sp.exp(solution[0])),
        inverse_coupling_squared=inverse_coupling,
        inferred_coupling_squared=inferred,
        positive_coupling_admissible=inverse_coupling.is_positive,
        residuals=tuple(sp.simplify(value) for value in residual_vector),
    )
