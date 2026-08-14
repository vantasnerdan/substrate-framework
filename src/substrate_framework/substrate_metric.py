"""Exact flowing-substrate effective metric and impedance matching.

This module extends the accepted static 1+1 optical metric family
(C-OG-001) to the case of a flowing substrate. The constitutive fields are
a positive refractive-index field ``n(x, t)`` and a substrate flow-velocity
field ``V(x, t)`` defined on the background Euclidean space and Newtonian
time of the continuum. The impedance-matching identity
(``rho*Theta = Z_0**2 = const``) is a separate algebraic statement that the
module records as a building block for any subsequent no-reflection
propagation theorem.

The metric follows the paper's eq. (69):

    g_tt = -c0**2/n + n * V**2
    g_tX = n * V
    g_XX = n

with determinant ``-c0**2``, inverse

    g^tt = -n/c0**2
    g^tX = n*V/c0**2
    g^XX = V**2 - n*V**2/c0**2

The ``n=1, V=0`` limit reduces to the homogeneous non-flowing Minkowski
metric. The ``V=0`` limit reduces to the static 1+1 optical metric family
in the paper's ``c0**2``-explicit convention. The massless null constraint
yields two characteristic speeds ``v = -V +/- c0/n`` in laboratory time.

The Christoffel-symbol path delegates to
:func:`substrate_framework.pseudo_riemannian.metric_christoffel_from_derivatives`
so the same Levi-Civita connection machinery used by the unrelated accepted
covariant sine-Gordon and Gordon modules is the single source of truth for
the connection. The pseudo-Riemannian module also supplies the metric
inverse and Ricci scalar hooks; this module does not re-derive them.

This module does not derive an Einstein source, physical gravity, a
material action, or a no-reflection propagation theorem. It records the
metric algebra and the impedance identity as reusable exact objects.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .pseudo_riemannian import (
    exact_metric_matrix,
    metric_christoffel_from_derivatives,
    metric_inverse,
    metric_ricci_scalar,
)


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) > 0.0:
            raise ValueError(f"{name} must be real and positive")
    return expression


def _real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be declared real")
    return expression


def impedance_product(index: Any, density_reference: Any, stiffness_reference: Any) -> sp.Expr:
    """Return ``rho * Theta`` under the constitutive map
    ``rho = rho_0 * n``, ``Theta = Theta_0 / n``.
    """

    n_value = _positive(index, "index")
    rho_0 = _positive(density_reference, "density_reference")
    theta_0 = _positive(stiffness_reference, "stiffness_reference")
    return sp.simplify((rho_0 * n_value) * (theta_0 / n_value))


def impedance_product_residual(index: Any, density_reference: Any, stiffness_reference: Any) -> sp.Expr:
    """Return the deviation of ``rho*Theta`` from ``rho_0*Theta_0``.

    Under the declared constitutive map this residual vanishes identically
    in any positive ``n``. The vanishing is the load-bearing identity of the
    impedance-matching condition.
    """

    product = impedance_product(index, density_reference, stiffness_reference)
    rho_0 = _positive(density_reference, "density_reference")
    theta_0 = _positive(stiffness_reference, "stiffness_reference")
    return sp.simplify(product - rho_0 * theta_0)


@dataclass(frozen=True)
class FlowingSubstrateMetric:
    """Exact 1+1 metric components and structural invariants.

    Coordinates are ordered ``(t, x)``. Signature is ``(-, +)`` with the
    paper's ``c0**2``-explicit convention. The metric is symmetric; the
    off-diagonal component is stored as a single scalar.
    """

    refractive_index: sp.Expr
    flow_velocity: sp.Expr
    signal_speed: sp.Expr
    covariant: sp.Matrix
    inverse: sp.Matrix
    determinant: sp.Expr
    massless_speeds: tuple


def _covariant_metric_components(index: Any, flow_velocity: Any, signal_speed: Any) -> tuple:
    """Return the three independent covariant components
    ``(g_tt, g_tX, g_XX)`` of the flowing substrate metric.
    """

    n_value = _positive(index, "index")
    v_value = _real(flow_velocity, "flow_velocity")
    c0 = _positive(signal_speed, "signal_speed")
    g_tt = -c0**2 / n_value + n_value * v_value**2
    g_tX = n_value * v_value
    g_XX = n_value
    return g_tt, g_tX, g_XX, c0


def flowing_substrate_metric(index: Any, flow_velocity: Any, signal_speed: Any) -> FlowingSubstrateMetric:
    """Return the exact flowing effective metric with inverse and structural data.

    The metric components follow the paper's eq. (69):
    ``g_tt = -c0**2/n + n*V**2``, ``g_tX = n*V``, ``g_XX = n``.
    Determinant is exactly ``-c0**2``. The two massless characteristic
    speeds in laboratory time are ``v = -V +/- c0/n`` (signed to match
    the paper's eq. (105) under the paper's sign convention).
    """

    g_tt, g_tX, g_XX, c0 = _covariant_metric_components(index, flow_velocity, signal_speed)
    n_value = _positive(index, "index")
    v_value = _real(flow_velocity, "flow_velocity")

    covariant = sp.Matrix([[g_tt, g_tX], [g_tX, g_XX]])
    metric_input = exact_metric_matrix(covariant, dimension=2)
    inverse = metric_inverse(metric_input)
    det = sp.simplify(metric_input.det())
    if det != -c0**2:
        raise ValueError(f"flowing metric determinant is {det}, expected -c0**2")

    v_symbol = sp.symbols("v_local", real=True)
    null_constraint = g_tt + 2 * g_tX * v_symbol + g_XX * v_symbol**2
    speeds = tuple(
        sp.simplify(s) for s in sp.solve(sp.Eq(null_constraint, 0), v_symbol)
    )

    return FlowingSubstrateMetric(
        refractive_index=n_value,
        flow_velocity=v_value,
        signal_speed=c0,
        covariant=metric_input,
        inverse=inverse,
        determinant=det,
        massless_speeds=speeds,
    )


def flowing_metric_determinant(index: Any, flow_velocity: Any, signal_speed: Any) -> sp.Expr:
    """Return the determinant of the flowing substrate metric."""

    return flowing_substrate_metric(index, flow_velocity, signal_speed).determinant


def flowing_metric_inverse(index: Any, flow_velocity: Any, signal_speed: Any) -> sp.Matrix:
    """Return the inverse of the flowing substrate metric."""

    return flowing_substrate_metric(index, flow_velocity, signal_speed).inverse


def flowing_metric_massless_speeds(index: Any, flow_velocity: Any, signal_speed: Any) -> tuple:
    """Return the two massless characteristic speeds in laboratory time.

    These are the values of ``v = dx/dt`` that satisfy the null constraint
    ``g_munu xdot^mu xdot^nu = 0`` under the laboratory-time parametrization
    ``t_dot = 1``. The result is ``v = -V + c0/n`` and ``v = -V - c0/n``
    in the paper's sign convention.
    """

    return flowing_substrate_metric(index, flow_velocity, signal_speed).massless_speeds


def flowing_metric_minkowski_limit(signal_speed: Any) -> FlowingSubstrateMetric:
    """Return the metric at ``n=1, V=0``: the homogeneous Minkowski metric."""

    c0 = _positive(signal_speed, "signal_speed")
    return flowing_substrate_metric(sp.Integer(1), sp.Integer(0), c0)


def flowing_metric_static_limit(index: Any, signal_speed: Any) -> FlowingSubstrateMetric:
    """Return the metric at ``V=0``: the substrate-framework's static
    1+1 optical metric in the paper's ``c0**2``-explicit convention.
    """

    c0 = _positive(signal_speed, "signal_speed")
    return flowing_substrate_metric(index, sp.Integer(0), c0)


def flowing_metric_ricci_scalar(
    index: Any,
    flow_velocity: Any,
    signal_speed: Any,
    time: sp.Symbol,
    space: sp.Symbol,
) -> sp.Expr:
    """Return the exact Ricci scalar of the flowing substrate metric in the
    laboratory ``(t, x)`` frame.

    Delegates to
    :func:`substrate_framework.pseudo_riemannian.metric_ricci_scalar` so the
    accepted connection machinery is the single source of truth.
    """

    metric_obj = flowing_substrate_metric(index, flow_velocity, signal_speed)
    return metric_ricci_scalar(metric_obj.covariant, (time, space))


def flowing_christoffel_lab_frame(
    index: Any,
    flow_velocity: Any,
    signal_speed: Any,
    time: sp.Symbol,
    space: sp.Symbol,
) -> dict:
    """Return the six non-trivial Christoffel symbols of the flowing substrate
    metric in the laboratory ``(t, x)`` frame as exact expressions in
    ``n``, ``V``, ``c0``, and the first derivatives ``partial_t n``,
    ``partial_x n``, ``partial_t V``, ``partial_x V`` evaluated along the
    world-line.

    Delegates to
    :func:`substrate_framework.pseudo_riemannian.metric_christoffel_from_derivatives`
    so the same Levi-Civita connection machinery used by the accepted
    covariant sine-Gordon and Gordon modules is the single source of truth.
    """

    n_value = _positive(index, "index")
    v_value = _real(flow_velocity, "flow_velocity")
    c0 = _positive(signal_speed, "signal_speed")

    g_tt = -c0**2 / n_value + n_value * v_value**2
    g_tX = n_value * v_value
    g_XX = n_value
    covariant = sp.Matrix([[g_tt, g_tX], [g_tX, g_XX]])
    inverse = metric_inverse(exact_metric_matrix(covariant, dimension=2))

    dtn = sp.diff(n_value, time)
    dxn = sp.diff(n_value, space)
    dtv = sp.diff(v_value, time)
    dxv = sp.diff(v_value, space)

    # Metric-component derivatives via the chain rule on (n, V).
    # dg_tt = (c0**2/n**2 + V**2) dn + 2 n V dV
    # dg_tX = V dn + n dV
    # dg_XX = dn
    dg_tt_dt = (c0**2 / n_value**2 + v_value**2) * dtn + 2 * n_value * v_value * dtv
    dg_tt_dx = (c0**2 / n_value**2 + v_value**2) * dxn + 2 * n_value * v_value * dxv
    dg_tX_dt = v_value * dtn + n_value * dtv
    dg_tX_dx = v_value * dxn + n_value * dxv
    dg_XX_dt = dtn
    dg_XX_dx = dxn

    # derivatives tensor has shape (rho, mu, nu) with (mu, nu) = (metric index, metric index)
    # coord 0 = t, coord 1 = x
    derivatives = [
        [
            [dg_tt_dt, dg_tX_dt],
            [dg_tX_dt, dg_XX_dt],
        ],
        [
            [dg_tt_dx, dg_tX_dx],
            [dg_tX_dx, dg_XX_dx],
        ],
    ]

    gamma = metric_christoffel_from_derivatives(inverse, derivatives)

    # gamma has shape (dim, dim, dim) with indices (upper, mu, nu) where
    # upper=0=t, upper=1=x; mu, nu are covariant indices (0=t, 1=x).
    return {
        "Gamma^t_tt": sp.simplify(gamma[0, 0, 0]),
        "Gamma^t_tX": sp.simplify(gamma[0, 0, 1]),
        "Gamma^t_XX": sp.simplify(gamma[0, 1, 1]),
        "Gamma^X_tt": sp.simplify(gamma[1, 0, 0]),
        "Gamma^X_tX": sp.simplify(gamma[1, 0, 1]),
        "Gamma^X_XX": sp.simplify(gamma[1, 1, 1]),
    }
