"""Exact ledgers for a sine-Gordon equation in mixed physical coordinates.

The conditional equation is ``theta_z_tau = g*sin(theta)``, where ``z`` is a
length coordinate, ``tau`` is a time coordinate, and ``g`` therefore has
dimensions ``1/(length*time)``. This module derives the coordinate map and
linear characteristic. It does not derive ``g`` from an optical medium,
identify an absorption coefficient with a frequency, or select a material.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True or expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive and real")
    return expression


@dataclass(frozen=True)
class MixedSineGordonScaleChoice:
    """One member of the normalization family ``g*L*T=1``."""

    coefficient: sp.Expr
    length_scale: sp.Expr
    time_scale: sp.Expr
    inverse_length_scale: sp.Expr
    inverse_time_scale: sp.Expr


@dataclass(frozen=True)
class MixedSineGordonLinearSpectrum:
    """Positive-``k`` branch of the mixed-coordinate characteristic."""

    wavenumber: sp.Expr
    angular_frequency: sp.Expr
    phase_velocity: sp.Expr
    group_velocity: sp.Expr


def mixed_sine_gordon_residual(
    field: Any,
    coordinate: sp.Symbol,
    retarded_time: sp.Symbol,
    coefficient: Any,
) -> sp.Expr:
    """Return ``theta_z_tau-g*sin(theta)``."""

    expression = sp.sympify(field)
    g = _positive_exact(coefficient, "coefficient")
    return sp.diff(expression, coordinate, retarded_time) - g * sp.sin(expression)


def mixed_sine_gordon_linearized_residual(
    field: Any,
    coordinate: sp.Symbol,
    retarded_time: sp.Symbol,
    coefficient: Any,
) -> sp.Expr:
    """Return the exact vacuum linearization ``theta_z_tau-g*theta``."""

    expression = sp.sympify(field)
    g = _positive_exact(coefficient, "coefficient")
    return sp.diff(expression, coordinate, retarded_time) - g * expression


def mixed_sine_gordon_linear_spectrum(
    wavenumber: Any,
    coefficient: Any,
) -> MixedSineGordonLinearSpectrum:
    r"""Return the branch ``Omega=g/k`` derived from ``k*Omega=g``.

    This characteristic is not ``Omega**2=omega_0**2+c**2*k**2`` and has no
    finite, wavenumber-independent angular-frequency floor.
    """

    k = _positive_exact(wavenumber, "wavenumber")
    g = _positive_exact(coefficient, "coefficient")
    angular = g / k
    return MixedSineGordonLinearSpectrum(
        wavenumber=k,
        angular_frequency=angular,
        phase_velocity=angular / k,
        group_velocity=-g / k**2,
    )


def mixed_sine_gordon_dimensionless_coupling(
    coefficient: Any,
    length_scale: Any,
    time_scale: Any,
) -> sp.Expr:
    """Return the normalized coefficient ``g*L*T``."""

    g = _positive_exact(coefficient, "coefficient")
    length = _positive_exact(length_scale, "length_scale")
    time = _positive_exact(time_scale, "time_scale")
    return g * length * time


def mixed_sine_gordon_scale_choice(
    coefficient: Any,
    length_scale: Any,
) -> MixedSineGordonScaleChoice:
    """Choose one length and derive the time required by ``g*L*T=1``.

    A fixed ``g`` leaves ``L`` arbitrary and hence does not select either an
    inverse-time scale or an inverse-length scale.
    """

    g = _positive_exact(coefficient, "coefficient")
    length = _positive_exact(length_scale, "length_scale")
    time = 1 / (g * length)
    return MixedSineGordonScaleChoice(
        coefficient=g,
        length_scale=length,
        time_scale=time,
        inverse_length_scale=1 / length,
        inverse_time_scale=1 / time,
    )


def mixed_sine_gordon_log_scale_jacobian() -> sp.ImmutableMatrix:
    r"""Return ``d log(g*L*T)/d(log L,log T) = (1,1)``.

    At fixed ``g`` this rank-one map has nullspace spanned by ``(1,-1)``, the
    reciprocal rescaling ``L -> rho*L``, ``T -> T/rho``.
    """

    return sp.ImmutableMatrix([[1, 1]])


def mixed_sine_gordon_hyperbolic_coordinates(
    coordinate: Any,
    retarded_time: Any,
    length_scale: Any,
    time_scale: Any,
) -> tuple[sp.Expr, sp.Expr]:
    r"""Return ``X=z/L+tau/T`` and ``S=z/L-tau/T``.

    When ``g*L*T=1``, ``theta_z_tau=g*sin(theta)`` becomes
    ``theta_SS-theta_XX+sin(theta)=0``.
    """

    length = _positive_exact(length_scale, "length_scale")
    time = _positive_exact(time_scale, "time_scale")
    xi = sp.sympify(coordinate) / length
    eta = sp.sympify(retarded_time) / time
    return sp.simplify(xi + eta), sp.simplify(xi - eta)


def mixed_sine_gordon_physical_coordinates(
    hyperbolic_space: Any,
    hyperbolic_time: Any,
    length_scale: Any,
    time_scale: Any,
) -> tuple[sp.Expr, sp.Expr]:
    """Invert :func:`mixed_sine_gordon_hyperbolic_coordinates`."""

    length = _positive_exact(length_scale, "length_scale")
    time = _positive_exact(time_scale, "time_scale")
    space = sp.sympify(hyperbolic_space)
    clock = sp.sympify(hyperbolic_time)
    return (
        sp.simplify(length * (space + clock) / 2),
        sp.simplify(time * (space - clock) / 2),
    )


def normalized_hyperbolic_sine_gordon_residual(
    field: Any,
    hyperbolic_space: sp.Symbol,
    hyperbolic_time: sp.Symbol,
) -> sp.Expr:
    """Return ``theta_SS-theta_XX+sin(theta)``."""

    expression = sp.sympify(field)
    return (
        sp.diff(expression, hyperbolic_time, 2)
        - sp.diff(expression, hyperbolic_space, 2)
        + sp.sin(expression)
    )


def mixed_coefficient_from_absorption_rate(
    absorption_coefficient: Any,
    inverse_time_scale: Any,
    prefactor: Any = 1,
) -> sp.Expr:
    """Dimensionally complete ``g=factor*alpha*rate``.

    This helper records only the required dimensions. It does not derive the
    rate or dimensionless prefactor from Maxwell-Bloch dynamics.
    """

    alpha = _positive_exact(absorption_coefficient, "absorption_coefficient")
    rate = _positive_exact(inverse_time_scale, "inverse_time_scale")
    factor = _positive_exact(prefactor, "prefactor")
    return factor * alpha * rate


def mixed_sine_gordon_dimension_matrix() -> sp.ImmutableMatrix:
    r"""Return dimensions over rows ``(length,time)``.

    Columns are ``(g,alpha,rate,omega_squared)``. They distinguish the mixed
    coefficient ``1/(L*T)``, absorption ``1/L``, rate ``1/T``, and a
    laboratory angular-frequency squared ``1/T**2``.
    """

    return sp.ImmutableMatrix(
        [
            [-1, -1, 0, 0],
            [-1, 0, -1, -2],
        ]
    )
