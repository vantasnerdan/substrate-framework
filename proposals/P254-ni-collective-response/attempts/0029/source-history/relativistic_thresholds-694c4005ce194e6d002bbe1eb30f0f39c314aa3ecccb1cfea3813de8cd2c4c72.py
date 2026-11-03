"""Exact two-body threshold four-momentum ledgers.

Subtracting one on-shell particle from a center-of-mass threshold vector
usually does not leave a second on-shell particle.  This module exposes the
exact residual and its mass-shell defect.  It does not derive a scattering
channel, boundary dynamics, detector response, interaction, or particle
identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return sp.simplify(expression)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


@dataclass(frozen=True)
class TwoBodyThresholdLedger:
    """Exact residual after observing one particle at two-body threshold.

    Components use signature ``(+,-)`` and ordering ``(energy, momentum)``.
    The total is the center-of-mass threshold vector ``(m1+m2,0)``; the
    observed particle has rapidity ``theta`` and is on shell with mass ``m1``.
    """

    observed_mass: sp.Expr
    residual_target_mass: sp.Expr
    observed_rapidity: sp.Expr
    threshold_four_momentum: sp.ImmutableMatrix
    observed_four_momentum: sp.ImmutableMatrix
    residual_four_momentum: sp.ImmutableMatrix
    observed_mass_shell_residual: sp.Expr
    four_momentum_closure: sp.ImmutableMatrix
    residual_invariant_mass_squared: sp.Expr
    residual_mass_shell_defect: sp.Expr


def two_body_threshold_ledger(
    observed_mass: Any,
    residual_target_mass: Any,
    observed_rapidity: Any,
) -> TwoBodyThresholdLedger:
    """Return the exact threshold residual for one observed on-shell particle.

    For positive exact masses ``m1,m2`` and exact real rapidity ``theta``, the
    returned defect is

    ``(P-p1)^2-m2^2 = 2*m1*(m1+m2)*(1-cosh(theta))``.

    Since ``cosh(theta) >= 1`` for real ``theta``, the residual is on shell at
    the target mass only at ``theta=0``.  Nonzero recoil therefore requires an
    above-threshold total or a separately modeled channel.
    """

    mass1 = _positive_exact(observed_mass, "observed_mass")
    mass2 = _positive_exact(residual_target_mass, "residual_target_mass")
    rapidity = _exact_real(observed_rapidity, "observed_rapidity")

    threshold = sp.ImmutableMatrix((mass1 + mass2, sp.S.Zero))
    observed = sp.ImmutableMatrix(
        (mass1 * sp.cosh(rapidity), mass1 * sp.sinh(rapidity))
    )
    residual = sp.ImmutableMatrix(
        tuple(sp.simplify(value) for value in threshold - observed)
    )
    observed_shell = sp.trigsimp(
        observed[0] ** 2 - observed[1] ** 2 - mass1**2
    )
    closure = sp.ImmutableMatrix(
        tuple(sp.simplify(value) for value in threshold - observed - residual)
    )
    residual_invariant = sp.trigsimp(residual[0] ** 2 - residual[1] ** 2)
    defect = sp.factor(sp.trigsimp(sp.expand(residual_invariant - mass2**2)))

    return TwoBodyThresholdLedger(
        observed_mass=mass1,
        residual_target_mass=mass2,
        observed_rapidity=rapidity,
        threshold_four_momentum=threshold,
        observed_four_momentum=observed,
        residual_four_momentum=residual,
        observed_mass_shell_residual=observed_shell,
        four_momentum_closure=closure,
        residual_invariant_mass_squared=residual_invariant,
        residual_mass_shell_defect=defect,
    )
