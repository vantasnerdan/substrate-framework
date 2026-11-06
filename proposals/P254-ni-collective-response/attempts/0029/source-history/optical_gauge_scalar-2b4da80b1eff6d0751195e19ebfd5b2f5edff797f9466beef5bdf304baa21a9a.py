"""Exact charged-scalar identities on the accepted 1+1 optical metric.

The module composes the declared optical metric and local-U(1) convention with
a separately supplied charged-scalar action.  It distinguishes local
gauge-coordinate labels from circle holonomy and boundary data.  No material
action, breather realization, electromagnetic field, analog-gravity
observation, or real-gravity coupling is inferred.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .optical_geometry import optical_metric_1d


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be declared real")
    return sp.simplify(expression)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be declared positive")
    return expression


def _exact_expression(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    return expression


def _coordinates(
    coordinates: tuple[sp.Symbol, sp.Symbol],
) -> tuple[sp.Symbol, sp.Symbol]:
    if len(coordinates) != 2:
        raise ValueError("coordinates must contain (time, space)")
    time, space = coordinates
    if (
        not isinstance(time, sp.Symbol)
        or not isinstance(space, sp.Symbol)
        or time.is_real is not True
        or space.is_real is not True
        or time == space
    ):
        raise ValueError("coordinates must be distinct real SymPy Symbols")
    return time, space


def _connections(
    connections: tuple[Any, Any],
) -> tuple[sp.Expr, sp.Expr]:
    if len(connections) != 2:
        raise ValueError("connections must contain (A_t, A_x)")
    return (
        _exact_real(connections[0], "temporal connection"),
        _exact_real(connections[1], "spatial connection"),
    )


def _charged_derivative(
    field: sp.Expr,
    connection: sp.Expr,
    coordinate: sp.Symbol,
    coupling: sp.Expr,
) -> sp.Expr:
    return sp.diff(field, coordinate) - sp.I * coupling * connection * field


@dataclass(frozen=True)
class ConstantOpticalGaugeDispersion:
    """Exact invariant momenta and mass shell for a constant background."""

    metric: sp.ImmutableMatrix
    inverse_metric: sp.ImmutableMatrix
    volume_density: sp.Expr
    invariant_frequency: sp.Expr
    invariant_wavenumber: sp.Expr
    mass_shell_lhs: sp.Expr
    invariant_frequency_squared: sp.Expr


@dataclass(frozen=True)
class AffineGaugePlaneWave:
    """Plane-wave labels after ``chi=a*t+b*x`` and their invariants."""

    frequency: sp.Expr
    wavenumber: sp.Expr
    temporal_connection: sp.Expr
    spatial_connection: sp.Expr
    invariant_frequency: sp.Expr
    invariant_wavenumber: sp.Expr


@dataclass(frozen=True)
class CircleOpticalGaugeMode:
    """One exact quasi-periodic circle mode in a constant flat connection."""

    mode: sp.Expr
    boundary_phase: sp.Expr
    canonical_wavenumber: sp.Expr
    invariant_wavenumber: sp.Expr
    invariant_frequency_squared: sp.Expr
    connection_holonomy: sp.Expr


def charged_optical_scalar_lagrangian_density(
    field: Any,
    conjugate_field: Any,
    index: Any,
    signal_speed: Any,
    mass: Any,
    coupling: Any,
    connections: tuple[Any, Any],
    coordinates: tuple[sp.Symbol, sp.Symbol],
) -> sp.Expr:
    """Return the exact coordinate density of the declared scalar action.

    The convention is
    ``-sqrt(-g)*(g^munu*(D_mu Psi)^*D_nu Psi + m^2*Psi^*Psi)`` with
    ``D=partial-i*e*A`` and ``g=diag(-1/n,n/c0^2)``.  ``conjugate_field`` is
    supplied independently so symbolic differentiation does not assume a
    particular complex-field representation.
    """

    time, space = _coordinates(coordinates)
    temporal, spatial = _connections(connections)
    psi = _exact_expression(field, "field")
    psi_bar = _exact_expression(conjugate_field, "conjugate field")
    n_value = _positive_exact(index, "index")
    c0 = _positive_exact(signal_speed, "signal_speed")
    mass_value = _positive_exact(mass, "mass")
    charge = _positive_exact(coupling, "coupling")
    if sp.diff(n_value, time) != 0:
        raise ValueError("index must be static")

    d_psi_t = _charged_derivative(psi, temporal, time, charge)
    d_psi_x = _charged_derivative(psi, spatial, space, charge)
    d_bar_t = sp.diff(psi_bar, time) + sp.I * charge * temporal * psi_bar
    d_bar_x = sp.diff(psi_bar, space) + sp.I * charge * spatial * psi_bar
    kinetic = -n_value * d_bar_t * d_psi_t + c0**2 * d_bar_x * d_psi_x / n_value
    return sp.simplify(-(kinetic + mass_value**2 * psi_bar * psi) / c0)


def charged_optical_scalar_euler_operator(
    field: Any,
    index: Any,
    signal_speed: Any,
    mass: Any,
    coupling: Any,
    connections: tuple[Any, Any],
    coordinates: tuple[sp.Symbol, sp.Symbol],
) -> sp.Expr:
    """Return ``Box_A Psi-m^2*Psi`` from the declared optical action.

    Because ``sqrt(-g)=1/c0`` is constant, the exact operator is
    ``D_t[-n D_t Psi]+D_x[(c0^2/n)D_x Psi]-m^2 Psi``.  This divergence form
    remains valid for a static varying index and arbitrary real connection.
    """

    time, space = _coordinates(coordinates)
    temporal, spatial = _connections(connections)
    psi = _exact_expression(field, "field")
    n_value = _positive_exact(index, "index")
    c0 = _positive_exact(signal_speed, "signal_speed")
    mass_value = _positive_exact(mass, "mass")
    charge = _positive_exact(coupling, "coupling")
    if sp.diff(n_value, time) != 0:
        raise ValueError("index must be static")

    d_psi_t = _charged_derivative(psi, temporal, time, charge)
    d_psi_x = _charged_derivative(psi, spatial, space, charge)
    temporal_flux = -n_value * d_psi_t
    spatial_flux = c0**2 * d_psi_x / n_value
    return sp.simplify(
        _charged_derivative(temporal_flux, temporal, time, charge)
        + _charged_derivative(spatial_flux, spatial, space, charge)
        - mass_value**2 * psi
    )


def constant_optical_gauge_dispersion(
    index: Any,
    signal_speed: Any,
    mass: Any,
    coupling: Any,
    frequency: Any,
    wavenumber: Any,
    temporal_connection: Any = 0,
    spatial_connection: Any = 0,
) -> ConstantOpticalGaugeDispersion:
    """Return the exact constant-background gauge-invariant mass shell."""

    n_value = _positive_exact(index, "index")
    c0 = _positive_exact(signal_speed, "signal_speed")
    mass_value = _positive_exact(mass, "mass")
    charge = _positive_exact(coupling, "coupling")
    omega = _exact_real(frequency, "frequency")
    wave_number = _exact_real(wavenumber, "wavenumber")
    temporal = _exact_real(temporal_connection, "temporal connection")
    spatial = _exact_real(spatial_connection, "spatial connection")
    metric = sp.ImmutableMatrix(optical_metric_1d(n_value, c0))
    inverse = sp.ImmutableMatrix(metric.inv().applyfunc(sp.simplify))
    invariant_frequency = sp.simplify(omega + charge * temporal)
    invariant_wavenumber = sp.simplify(wave_number - charge * spatial)
    mass_shell_lhs = sp.simplify(
        n_value * invariant_frequency**2
        - c0**2 * invariant_wavenumber**2 / n_value
    )
    return ConstantOpticalGaugeDispersion(
        metric=metric,
        inverse_metric=inverse,
        volume_density=sp.simplify(sp.sqrt(-metric.det())),
        invariant_frequency=invariant_frequency,
        invariant_wavenumber=invariant_wavenumber,
        mass_shell_lhs=mass_shell_lhs,
        invariant_frequency_squared=sp.simplify(
            mass_value**2 / n_value
            + c0**2 * invariant_wavenumber**2 / n_value**2
        ),
    )


def affine_gauge_plane_wave(
    frequency: Any,
    wavenumber: Any,
    temporal_connection: Any,
    spatial_connection: Any,
    coupling: Any,
    time_slope: Any,
    space_slope: Any,
) -> AffineGaugePlaneWave:
    """Transform labels for ``chi=time_slope*t+space_slope*x`` exactly."""

    omega = _exact_real(frequency, "frequency")
    wave_number = _exact_real(wavenumber, "wavenumber")
    temporal = _exact_real(temporal_connection, "temporal connection")
    spatial = _exact_real(spatial_connection, "spatial connection")
    charge = _positive_exact(coupling, "coupling")
    time_gradient = _exact_real(time_slope, "time_slope")
    space_gradient = _exact_real(space_slope, "space_slope")
    transformed_frequency = sp.simplify(omega - charge * time_gradient)
    transformed_wavenumber = sp.simplify(wave_number + charge * space_gradient)
    transformed_temporal = sp.simplify(temporal + time_gradient)
    transformed_spatial = sp.simplify(spatial + space_gradient)
    return AffineGaugePlaneWave(
        frequency=transformed_frequency,
        wavenumber=transformed_wavenumber,
        temporal_connection=transformed_temporal,
        spatial_connection=transformed_spatial,
        invariant_frequency=sp.simplify(
            transformed_frequency + charge * transformed_temporal
        ),
        invariant_wavenumber=sp.simplify(
            transformed_wavenumber - charge * transformed_spatial
        ),
    )


def circle_optical_gauge_mode(
    index: Any,
    signal_speed: Any,
    mass: Any,
    coupling: Any,
    spatial_connection: Any,
    circumference: Any,
    mode: Any,
    boundary_phase: Any = 0,
) -> CircleOpticalGaugeMode:
    """Return one mode for ``Psi(x+L)=exp(i*theta)Psi(x)``."""

    n_value = _positive_exact(index, "index")
    c0 = _positive_exact(signal_speed, "signal_speed")
    mass_value = _positive_exact(mass, "mass")
    charge = _positive_exact(coupling, "coupling")
    spatial = _exact_real(spatial_connection, "spatial connection")
    length = _positive_exact(circumference, "circumference")
    phase = _exact_real(boundary_phase, "boundary_phase")
    integer_mode = sp.sympify(mode)
    if integer_mode.has(sp.Float) or integer_mode.is_integer is not True:
        raise ValueError("mode must be an exact integer")
    canonical_wavenumber = sp.simplify((2 * sp.pi * integer_mode + phase) / length)
    invariant_wavenumber = sp.simplify(canonical_wavenumber - charge * spatial)
    return CircleOpticalGaugeMode(
        mode=integer_mode,
        boundary_phase=phase,
        canonical_wavenumber=canonical_wavenumber,
        invariant_wavenumber=invariant_wavenumber,
        invariant_frequency_squared=sp.simplify(
            mass_value**2 / n_value
            + c0**2 * invariant_wavenumber**2 / n_value**2
        ),
        connection_holonomy=sp.simplify(sp.exp(sp.I * charge * spatial * length)),
    )


def berry_one_form_to_u1_connection(
    berry_component: Any,
    texture_coordinate: Any,
    spatial_coordinate: sp.Symbol,
    coupling: Any,
) -> sp.Expr:
    """Return the same-phase U(1) component from a Berry one-form pullback.

    For ``B=i*psi^dagger*d psi``, ``psi'=exp(i*chi)psi`` gives
    ``B'=B-d chi``.  Matching the same phase to
    ``Psi'=exp(i*e*lambda)Psi`` therefore requires ``e*A=-B``.  Pulling
    ``B_phi*dphi`` back along a separately declared texture ``phi(x)`` gives
    ``A_x=-B_phi*partial_x(phi)/e``.
    """

    if (
        not isinstance(spatial_coordinate, sp.Symbol)
        or spatial_coordinate.is_real is not True
    ):
        raise ValueError("spatial_coordinate must be a real SymPy Symbol")
    component = _exact_real(berry_component, "berry_component")
    texture = _exact_real(texture_coordinate, "texture_coordinate")
    charge = _positive_exact(coupling, "coupling")
    return sp.simplify(-component * sp.diff(texture, spatial_coordinate) / charge)
