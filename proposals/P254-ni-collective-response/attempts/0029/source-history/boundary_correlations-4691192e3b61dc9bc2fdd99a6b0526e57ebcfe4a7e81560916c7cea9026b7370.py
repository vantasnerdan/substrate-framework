"""Exact boundary sign correlations and half-line winding conversion.

The sign correlation defined here is a boundary observable, not a topological
charge.  The two are kept as separately named APIs because no algebraic
implication relates them without additional boundary dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


@dataclass(frozen=True)
class ScalarBoundaryParityLedger:
    """Exact parity decomposition of one scalar boundary residual."""

    temporal_trace: sp.Expr
    coordinate_trace: sp.Expr
    temporal_coefficient: sp.Expr
    coordinate_coefficient: sp.Expr
    source: sp.Expr
    residual: sp.Expr
    parity_even_component: sp.Expr
    parity_odd_component: sp.Expr
    parity_image_residual: sp.Expr
    reflected_coefficient_residual: sp.Expr
    fixed_parameter_parity_defect: sp.Expr


@dataclass(frozen=True)
class OrientedHalfLineParityLedger:
    """Parity map from a right half-line to its left-half-line image."""

    temporal_trace: sp.Expr
    right_coordinate_trace: sp.Expr
    right_outward_trace: sp.Expr
    left_parity_coordinate_trace: sp.Expr
    left_outward_trace: sp.Expr
    right_residual: sp.Expr
    left_parity_residual: sp.Expr


@dataclass(frozen=True)
class ScalarBoundaryTraceFamily:
    """Solution family of one linear residual for two boundary traces."""

    temporal_trace_parameter: sp.Expr
    coordinate_trace_solution: sp.Expr | None
    temporal_only_constraint: sp.Expr | None
    coordinate_trace_free: bool


def _real_scalar(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_real is not True:
        raise ValueError(f"{name} must be real")
    return expression


def _exact_real_scalar(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise TypeError(f"{name} must be exact")
    if expression.is_real is False:
        raise ValueError(f"{name} must be real")
    return expression


def scalar_boundary_parity_ledger(
    temporal_trace: Any,
    coordinate_trace: Any,
    temporal_coefficient: Any,
    coordinate_coefficient: Any,
    source: Any,
) -> ScalarBoundaryParityLedger:
    """Return the exact scalar-parity ledger for ``a*u + beta*v - J``.

    The supplied traces represent ``u=phi_t(t,b)`` and ``v=phi_x(t,b)``.
    At the reflected point, scalar spatial parity leaves ``u`` unchanged and
    changes ``v`` to ``-v``.  A scalar source is pulled back without an
    intrinsic sign.  The result therefore maps the coefficient family
    ``beta`` to ``-beta``.  A mixed residual is not itself a parity-odd
    eigenobject: it contains the even part ``a*u-J`` and odd part ``beta*v``.
    """

    u = _exact_real_scalar(temporal_trace, "temporal_trace")
    v = _exact_real_scalar(coordinate_trace, "coordinate_trace")
    a = _exact_real_scalar(temporal_coefficient, "temporal_coefficient")
    beta = _exact_real_scalar(coordinate_coefficient, "coordinate_coefficient")
    current = _exact_real_scalar(source, "source")
    even = sp.simplify(a * u - current)
    odd = sp.simplify(beta * v)
    residual = sp.simplify(even + odd)
    parity_image = sp.simplify(even - odd)
    reflected_coefficient = sp.simplify(a * u - beta * v - current)
    return ScalarBoundaryParityLedger(
        temporal_trace=u,
        coordinate_trace=v,
        temporal_coefficient=a,
        coordinate_coefficient=beta,
        source=current,
        residual=residual,
        parity_even_component=even,
        parity_odd_component=odd,
        parity_image_residual=parity_image,
        reflected_coefficient_residual=reflected_coefficient,
        fixed_parameter_parity_defect=sp.simplify(parity_image - residual),
    )


def oriented_half_line_parity_ledger(
    temporal_trace: Any,
    right_coordinate_trace: Any,
    temporal_coefficient: Any,
    normal_coefficient: Any,
    source: Any,
) -> OrientedHalfLineParityLedger:
    """Return the right-to-left half-line parity map with outward normals.

    For a right half-line ``x>=b``, the left-boundary outward derivative is
    ``-phi_x``.  Parity maps it to a left half-line ``x<=-b`` whose
    right-boundary outward derivative is ``+partial_x``; the transformed
    scalar field has coordinate trace ``-phi_x``.  The two outward-normal
    traces, and hence the normal coefficient, are unchanged by the combined
    field-and-domain map.
    """

    u = _exact_real_scalar(temporal_trace, "temporal_trace")
    v = _exact_real_scalar(right_coordinate_trace, "right_coordinate_trace")
    a = _exact_real_scalar(temporal_coefficient, "temporal_coefficient")
    eta = _exact_real_scalar(normal_coefficient, "normal_coefficient")
    current = _exact_real_scalar(source, "source")
    right_outward = sp.simplify(-v)
    left_coordinate = sp.simplify(-v)
    left_outward = left_coordinate
    right_residual = sp.simplify(a * u + eta * right_outward - current)
    left_residual = sp.simplify(a * u + eta * left_outward - current)
    return OrientedHalfLineParityLedger(
        temporal_trace=u,
        right_coordinate_trace=v,
        right_outward_trace=right_outward,
        left_parity_coordinate_trace=left_coordinate,
        left_outward_trace=left_outward,
        right_residual=right_residual,
        left_parity_residual=left_residual,
    )


def scalar_boundary_trace_family(
    temporal_trace_parameter: Any,
    temporal_coefficient: Any,
    coordinate_coefficient: Any,
    source: Any,
) -> ScalarBoundaryTraceFamily:
    """Solve one exact linear boundary residual without hiding its freedom.

    For declared nonzero ``beta``, the relation ``a*u+beta*v=J`` gives
    ``v=(J-a*u)/beta`` while ``u`` remains a free family parameter.  For
    ``beta=0``, the coordinate trace is wholly unconstrained and only
    ``a*u-J=0`` remains.  A symbolic coefficient must be declared zero or
    nonzero so the branch is auditable.
    """

    u = _exact_real_scalar(temporal_trace_parameter, "temporal_trace_parameter")
    a = _exact_real_scalar(temporal_coefficient, "temporal_coefficient")
    beta = _exact_real_scalar(coordinate_coefficient, "coordinate_coefficient")
    current = _exact_real_scalar(source, "source")
    if beta.is_zero is None:
        raise ValueError("coordinate_coefficient must be declared zero or nonzero")
    if beta.is_zero:
        return ScalarBoundaryTraceFamily(
            temporal_trace_parameter=u,
            coordinate_trace_solution=None,
            temporal_only_constraint=sp.simplify(a * u - current),
            coordinate_trace_free=True,
        )
    return ScalarBoundaryTraceFamily(
        temporal_trace_parameter=u,
        coordinate_trace_solution=sp.simplify((current - a * u) / beta),
        temporal_only_constraint=None,
        coordinate_trace_free=False,
    )


def _positive_frequency(value: Any) -> sp.Expr:
    frequency = _real_scalar(value, "angular_frequency")
    if frequency.is_number and not float(frequency) > 0.0:
        raise ValueError("angular_frequency must be positive")
    return frequency


def boundary_sign_correlation_density(
    boundary_time_derivative: Any,
    boundary_coordinate_derivative: Any,
) -> sp.Expr:
    """Return ``sign(phi_t)*phi_x`` for declared boundary traces.

    The second argument is the derivative in a fixed coordinate direction. It
    is not automatically an outward-normal derivative of a physical domain.
    """

    time_derivative = _real_scalar(
        boundary_time_derivative,
        "boundary_time_derivative",
    )
    coordinate_derivative = _real_scalar(
        boundary_coordinate_derivative,
        "boundary_coordinate_derivative",
    )
    return sp.sign(time_derivative) * coordinate_derivative


def sinusoidal_boundary_sign_correlation(
    time_amplitude: Any,
    coordinate_amplitude: Any,
    angular_frequency: Any,
    relative_phase: Any,
) -> sp.Expr:
    """Return the exact full-period correlation for two sinusoidal traces.

    The convention is ``phi_t=A*sin(theta)`` and
    ``phi_x=B*sin(theta+delta)`` over one period, where ``delta`` is the
    supplied relative phase.  The result is
    ``4*sign(A)*B*cos(delta)/omega``.  A cosine convention for the second
    trace requires the corresponding explicit phase shift.
    """

    temporal_amplitude = _real_scalar(time_amplitude, "time_amplitude")
    spatial_amplitude = _real_scalar(
        coordinate_amplitude,
        "coordinate_amplitude",
    )
    frequency = _positive_frequency(angular_frequency)
    phase = _real_scalar(relative_phase, "relative_phase")
    return sp.simplify(
        4 * sp.sign(temporal_amplitude) * spatial_amplitude * sp.cos(phase)
        / frequency
    )


def right_half_line_topological_charge_change(
    boundary_field_change: Any,
) -> sp.Expr:
    """Return ``-Delta(phi_boundary)/(2*pi)`` on a right half-line.

    This conversion uses the orientation ``epsilon**01=+1`` and domain
    ``x >= b``.  It additionally assumes the field at positive infinity is
    time independent, so the boundary field change is the complete integrated
    topological flux.
    """

    field_change = _real_scalar(boundary_field_change, "boundary_field_change")
    return sp.simplify(-field_change / (2 * sp.pi))
