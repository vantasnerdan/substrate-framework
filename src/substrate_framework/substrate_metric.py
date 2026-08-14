"""Conditional exact algebra for a declared flowing 1+1 metric.

Authority status: conditional, unpromoted infrastructure linked to issue #62
and source PR #63.  Callers declare a positive refractive index ``n(t, x)``,
a real flow field ``V(t, x)``, and a positive constant signal speed ``c0``.
The resulting covariant metric in laboratory coordinates ``(t, x)`` is

``g = [[-c0**2/n + n*V**2, n*V], [n*V, n]]``.

Because ``t`` and ``x`` carry time and length units, respectively, the three
independent components have the coordinate-dependent units needed to make
every term of ``ds**2`` a squared length.  At ``V=0`` this metric is exactly
``c0**2`` times the accepted C-OG-001 metric; it is not literal component
equality.  The constant rescaling preserves the Levi-Civita connection and
unparametrized geodesics while rescaling the Ricci scalar.

This module derives consequences of the displayed metric ansatz.  It does not
derive that ansatz from an accepted substrate action, promote C-OG-006, prove
reflection-free propagation, identify physical gravity, or supply an Einstein
source.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .pseudo_riemannian import (
    coordinate_symbols,
    exact_metric_matrix,
    metric_christoffel_symbols,
    metric_inverse,
    metric_ricci_scalar,
)


__all__ = [
    "FlowingSubstrateMetric",
    "flowing_christoffel_lab_frame",
    "flowing_metric_ricci_scalar",
    "flowing_substrate_metric",
]


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be declared real")
    return sp.simplify(expression)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be declared positive")
    return expression


def _constant_in_chart(
    expression: sp.Expr,
    coordinates: tuple[sp.Symbol, sp.Symbol],
    name: str,
) -> None:
    if any(
        sp.simplify(sp.diff(expression, coordinate)) != 0
        for coordinate in coordinates
    ):
        raise ValueError(f"{name} must be constant in the supplied coordinates")


@dataclass(frozen=True)
class FlowingSubstrateMetric:
    """Exact data for the declared flowing metric.

    Coordinates are ordered ``(t, x)`` and the signature is ``(-, +)``.
    ``massless_speeds`` are the two ``dx/dt`` roots in increasing branch-sign
    order, ``(-V-c0/n, -V+c0/n)``.
    """

    refractive_index: sp.Expr
    flow_velocity: sp.Expr
    signal_speed: sp.Expr
    covariant: sp.Matrix
    inverse: sp.Matrix
    determinant: sp.Expr
    massless_speeds: tuple[sp.Expr, sp.Expr]


def flowing_substrate_metric(
    index: Any,
    flow_velocity: Any,
    signal_speed: Any,
) -> FlowingSubstrateMetric:
    """Return exact algebraic data for the declared flowing metric.

    ``index`` must be declared positive, ``flow_velocity`` real, and all
    inputs exact.  ``signal_speed`` is a positive constant model parameter;
    chart-dependent uses enforce that constancy when coordinates are known.
    """

    n_value = _positive_exact(index, "index")
    velocity = _exact_real(flow_velocity, "flow_velocity")
    c0 = _positive_exact(signal_speed, "signal_speed")
    covariant = exact_metric_matrix(
        sp.Matrix(
            [
                [-c0**2 / n_value + n_value * velocity**2, n_value * velocity],
                [n_value * velocity, n_value],
            ]
        ),
        dimension=2,
    )
    determinant = sp.simplify(covariant.det())
    if sp.simplify(determinant + c0**2) != 0:
        raise AssertionError("declared flowing metric lost its determinant invariant")

    inverse = metric_inverse(covariant)
    massless_speeds = (
        sp.simplify(-velocity - c0 / n_value),
        sp.simplify(-velocity + c0 / n_value),
    )
    return FlowingSubstrateMetric(
        refractive_index=n_value,
        flow_velocity=velocity,
        signal_speed=c0,
        covariant=covariant,
        inverse=inverse,
        determinant=determinant,
        massless_speeds=massless_speeds,
    )


def flowing_christoffel_lab_frame(
    index: Any,
    flow_velocity: Any,
    signal_speed: Any,
    time: sp.Symbol,
    space: sp.Symbol,
) -> dict[str, sp.Expr]:
    """Return the six independent ``Gamma^a_bc`` components in ``(t, x)``.

    The fields may depend on both supplied coordinates.  The signal speed is
    the constant ``c0`` of the declared ansatz; a coordinate-dependent value
    is rejected rather than silently omitting its derivatives.
    """

    coordinates = coordinate_symbols((time, space))
    metric = flowing_substrate_metric(index, flow_velocity, signal_speed)
    _constant_in_chart(metric.signal_speed, coordinates, "signal_speed")
    gamma = metric_christoffel_symbols(metric.covariant, coordinates)
    return {
        "Gamma^t_tt": gamma[0][0][0],
        "Gamma^t_tX": gamma[0][0][1],
        "Gamma^t_XX": gamma[0][1][1],
        "Gamma^X_tt": gamma[1][0][0],
        "Gamma^X_tX": gamma[1][0][1],
        "Gamma^X_XX": gamma[1][1][1],
    }


def flowing_metric_ricci_scalar(
    index: Any,
    flow_velocity: Any,
    signal_speed: Any,
    time: sp.Symbol,
    space: sp.Symbol,
) -> sp.Expr:
    """Return the exact Ricci scalar of the declared metric in ``(t, x)``."""

    coordinates = coordinate_symbols((time, space))
    metric = flowing_substrate_metric(index, flow_velocity, signal_speed)
    _constant_in_chart(metric.signal_speed, coordinates, "signal_speed")
    return metric_ricci_scalar(metric.covariant, coordinates)
