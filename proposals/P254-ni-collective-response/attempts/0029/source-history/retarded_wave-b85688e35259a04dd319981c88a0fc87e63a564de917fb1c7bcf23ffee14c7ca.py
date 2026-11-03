"""Exact retarded radiation from a point source in one spatial dimension.

The declared scalar Lagrangian density is

``A*(phi_t**2-c**2*phi_x**2)/2 + B*phi*q(t)*delta(x)``.

This module derives the equation, jump, outgoing flux, and source-work ledger
from that one normalization. The retarded result additionally assumes no
incoming field and a source primitive with vanishing past boundary term. A
static countermodel is returned separately because the local sourced equation
alone does not select a radiating history.

Nothing here derives a dilaton action, a gravitational degree of freedom, a
breather source, a coupling value, or a radiation-reaction law.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _exact_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if result.is_real is not True:
        raise ValueError(f"{name} must be declared real")
    return sp.simplify(result)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    result = _exact_real(value, name)
    if result.is_positive is not True:
        raise ValueError(f"{name} must be declared positive")
    return result


@dataclass(frozen=True)
class RetardedPointSourceRadiation:
    """Exact outgoing ledger for the declared canonical scalar action.

    ``field_primitive_coefficient`` multiplies a primitive ``I`` satisfying
    ``I'(u)=q(u)`` in ``phi(t,x)=coefficient*I(t-|x|/c)``. The derivative and
    flux fields are evaluated at the appropriate retarded time.
    """

    kinetic_coefficient: sp.Expr
    source_coupling: sp.Expr
    wave_speed: sp.Expr
    source_amplitude: sp.Expr
    equation_delta_coefficient: sp.Expr
    field_primitive_coefficient: sp.Expr
    time_derivative: sp.Expr
    right_space_derivative: sp.Expr
    left_space_derivative: sp.Expr
    derivative_jump: sp.Expr
    right_outward_flux: sp.Expr
    left_outward_flux: sp.Expr
    total_outward_power: sp.Expr
    source_work_rate: sp.Expr


@dataclass(frozen=True)
class StaticPointSourceCountermodel:
    """Static solution of the same local point-source wave equation.

    The field is ``coefficient*abs(x)``. It has a nonzero asymptotic slope and
    therefore is not the retarded finite-past solution, but its identical
    derivative jump and zero flux expose the missing boundary/history premise.
    """

    kinetic_coefficient: sp.Expr
    source_coupling: sp.Expr
    wave_speed: sp.Expr
    source_amplitude: sp.Expr
    equation_delta_coefficient: sp.Expr
    absolute_value_coefficient: sp.Expr
    right_space_derivative: sp.Expr
    left_space_derivative: sp.Expr
    derivative_jump: sp.Expr
    total_outward_power: sp.Expr


def retarded_point_source_radiation(
    kinetic_coefficient: Any,
    source_coupling: Any,
    wave_speed: Any,
    source_amplitude: Any,
) -> RetardedPointSourceRadiation:
    r"""Derive the exact no-incoming retarded point-source flux ledger.

    Variation of the declared action gives
    ``phi_tt-c**2*phi_xx=(B/A)*q(t)*delta(x)``. With the retarded Green
    function and no past boundary contribution,
    ``phi=B*I(t-|x|/c)/(2*A*c)`` where ``I'=q``. The canonical energy flux is
    ``S=-A*c**2*phi_t*phi_x``. The two outgoing sides carry equal flux and
    their sum equals the local work ``B*q*phi_t(t,0)``.
    """

    coefficient = _positive_exact(kinetic_coefficient, "kinetic_coefficient")
    coupling = _exact_real(source_coupling, "source_coupling")
    speed = _positive_exact(wave_speed, "wave_speed")
    source = _exact_real(source_amplitude, "source_amplitude")

    equation_delta = sp.simplify(coupling * source / coefficient)
    primitive_coefficient = sp.simplify(coupling / (2 * coefficient * speed))
    time_derivative = sp.simplify(primitive_coefficient * source)
    right_space_derivative = sp.simplify(-time_derivative / speed)
    left_space_derivative = sp.simplify(time_derivative / speed)
    derivative_jump = sp.simplify(
        right_space_derivative - left_space_derivative
    )
    one_side_flux = sp.simplify(
        coefficient * speed**2 * time_derivative**2 / speed
    )
    total_power = sp.simplify(2 * one_side_flux)
    source_work = sp.simplify(coupling * source * time_derivative)

    if sp.simplify(-speed**2 * derivative_jump - equation_delta) != 0:
        raise AssertionError("retarded derivative jump does not reproduce source")
    if sp.simplify(total_power - source_work) != 0:
        raise AssertionError("outgoing flux and source work disagree")

    return RetardedPointSourceRadiation(
        kinetic_coefficient=coefficient,
        source_coupling=coupling,
        wave_speed=speed,
        source_amplitude=source,
        equation_delta_coefficient=equation_delta,
        field_primitive_coefficient=primitive_coefficient,
        time_derivative=time_derivative,
        right_space_derivative=right_space_derivative,
        left_space_derivative=left_space_derivative,
        derivative_jump=derivative_jump,
        right_outward_flux=one_side_flux,
        left_outward_flux=one_side_flux,
        total_outward_power=total_power,
        source_work_rate=source_work,
    )


def static_point_source_countermodel(
    kinetic_coefficient: Any,
    source_coupling: Any,
    wave_speed: Any,
    source_amplitude: Any,
) -> StaticPointSourceCountermodel:
    r"""Return a zero-flux static solution of the same sourced equation.

    For constant ``q``, ``phi=-B*q*abs(x)/(2*A*c**2)`` obeys the same delta
    source equation. It demonstrates that the local field equation and jump
    condition do not imply the retarded radiation result without boundary and
    history data.
    """

    coefficient = _positive_exact(kinetic_coefficient, "kinetic_coefficient")
    coupling = _exact_real(source_coupling, "source_coupling")
    speed = _positive_exact(wave_speed, "wave_speed")
    source = _exact_real(source_amplitude, "source_amplitude")

    equation_delta = sp.simplify(coupling * source / coefficient)
    absolute_value_coefficient = sp.simplify(
        -coupling * source / (2 * coefficient * speed**2)
    )
    right_space_derivative = absolute_value_coefficient
    left_space_derivative = -absolute_value_coefficient
    derivative_jump = sp.simplify(
        right_space_derivative - left_space_derivative
    )
    if sp.simplify(-speed**2 * derivative_jump - equation_delta) != 0:
        raise AssertionError("static derivative jump does not reproduce source")

    return StaticPointSourceCountermodel(
        kinetic_coefficient=coefficient,
        source_coupling=coupling,
        wave_speed=speed,
        source_amplitude=source,
        equation_delta_coefficient=equation_delta,
        absolute_value_coefficient=absolute_value_coefficient,
        right_space_derivative=right_space_derivative,
        left_space_derivative=left_space_derivative,
        derivative_jump=derivative_jump,
        total_outward_power=sp.Integer(0),
    )
