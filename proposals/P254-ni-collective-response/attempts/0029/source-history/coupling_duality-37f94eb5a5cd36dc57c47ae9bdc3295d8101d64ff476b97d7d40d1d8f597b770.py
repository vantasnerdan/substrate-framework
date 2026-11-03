"""Exact ledgers for conditional reciprocal coupling maps.

For separately supplied positive coordinates ``x`` and coefficients ``A``,
the map ``D_A(x)=A/x`` is an involution with positive fixed coordinate
``sqrt(A)``.  These algebraic facts do not derive ``A``, select the fixed
subfamily, or establish that any action, observable, or physical theory is
dual under the map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


@dataclass(frozen=True)
class ReciprocalCouplingLedger:
    """Exact orbit and fixed-point data for one supplied reciprocal map."""

    coupling_coordinate: sp.Expr
    duality_coefficient: sp.Expr
    dual_coordinate: sp.Expr
    double_dual_coordinate: sp.Expr
    orbit_product: sp.Expr
    positive_fixed_point: sp.Expr
    fixed_point_image: sp.Expr
    fixed_point_residual: sp.Expr


@dataclass(frozen=True)
class ReciprocalCoordinateChangeLedger:
    """Conjugation of ``D_A`` under the coordinate change ``x'=rho*x``."""

    original: ReciprocalCouplingLedger
    coordinate_rescaling: sp.Expr
    rescaled_coordinate: sp.Expr
    rescaled_duality_coefficient: sp.Expr
    rescaled_dual_coordinate: sp.Expr
    conjugated_dual_coordinate: sp.Expr
    rescaled_positive_fixed_point: sp.Expr
    conjugated_positive_fixed_point: sp.Expr
    coefficient_rescaling_ratio: sp.Expr
    fixed_point_rescaling_ratio: sp.Expr


def reciprocal_coupling_ledger(
    coupling_coordinate: Any,
    duality_coefficient: Any,
) -> ReciprocalCouplingLedger:
    """Return exact data for ``D_A(x)=A/x`` on the positive half-line.

    ``A`` is a load-bearing premise.  The returned fixed point describes the
    additional restriction ``x=D_A(x)``; generic positive ``x`` instead forms
    an off-fixed two-element orbit (or a singleton only at the fixed point).
    """

    coordinate = _positive(coupling_coordinate, "coupling_coordinate")
    coefficient = _positive(duality_coefficient, "duality_coefficient")
    dual = sp.simplify(coefficient / coordinate)
    double_dual = sp.simplify(coefficient / dual)
    fixed_point = sp.sqrt(coefficient)
    fixed_image = sp.simplify(coefficient / fixed_point)
    return ReciprocalCouplingLedger(
        coupling_coordinate=coordinate,
        duality_coefficient=coefficient,
        dual_coordinate=dual,
        double_dual_coordinate=double_dual,
        orbit_product=sp.simplify(coordinate * dual),
        positive_fixed_point=fixed_point,
        fixed_point_image=fixed_image,
        fixed_point_residual=sp.simplify(fixed_image - fixed_point),
    )


def reciprocal_coefficient_for_fixed_target(target_coordinate: Any) -> sp.Expr:
    """Return the coefficient that makes a supplied positive target fixed.

    This inverse construction is ``A=target**2``.  It shows that a desired
    fixed coordinate can be encoded by choosing ``A``; it does not predict the
    target or the coefficient.
    """

    target = _positive(target_coordinate, "target_coordinate")
    return sp.simplify(target**2)


def reciprocal_coordinate_change_ledger(
    coupling_coordinate: Any,
    duality_coefficient: Any,
    coordinate_rescaling: Any,
) -> ReciprocalCoordinateChangeLedger:
    """Conjugate a reciprocal map by ``x'=rho*x`` for positive ``rho``.

    The equivalent primed map is ``D'_(rho**2*A)(x')``.  Consequently its
    coefficient and numeric fixed coordinate become ``rho**2*A`` and
    ``rho*sqrt(A)`` while the involution and orbit-product structure persist.
    """

    original = reciprocal_coupling_ledger(
        coupling_coordinate,
        duality_coefficient,
    )
    rho = _positive(coordinate_rescaling, "coordinate_rescaling")
    rescaled_coordinate = sp.simplify(rho * original.coupling_coordinate)
    rescaled_coefficient = sp.simplify(
        rho**2 * original.duality_coefficient
    )
    rescaled = reciprocal_coupling_ledger(
        rescaled_coordinate,
        rescaled_coefficient,
    )
    return ReciprocalCoordinateChangeLedger(
        original=original,
        coordinate_rescaling=rho,
        rescaled_coordinate=rescaled_coordinate,
        rescaled_duality_coefficient=rescaled_coefficient,
        rescaled_dual_coordinate=rescaled.dual_coordinate,
        conjugated_dual_coordinate=sp.simplify(
            rho * original.dual_coordinate
        ),
        rescaled_positive_fixed_point=rescaled.positive_fixed_point,
        conjugated_positive_fixed_point=sp.simplify(
            rho * original.positive_fixed_point
        ),
        coefficient_rescaling_ratio=sp.simplify(
            rescaled_coefficient / original.duality_coefficient
        ),
        fixed_point_rescaling_ratio=sp.simplify(
            rescaled.positive_fixed_point / original.positive_fixed_point
        ),
    )
