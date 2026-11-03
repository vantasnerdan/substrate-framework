"""Exact Berry ledgers for closed rank-one projector paths.

The local one-form returned here depends on the chosen normalized section.
The closed-ray holonomy also includes that section's endpoint transition and
is gauge invariant.  These definitions supply no physical vector potential,
curvature source, defect dynamics, electromagnetic dictionary, material, or
observation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _integer(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_integer is not True:
        raise ValueError(f"{name} must be an integer")
    return expression


def _real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


def _column(section: Any) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(section)
    if matrix.cols != 1 or matrix.rows < 2:
        raise ValueError("section must be a column with dimension at least two")
    return matrix


def _zero(expression: Any) -> bool:
    if isinstance(expression, sp.MatrixBase):
        return all(sp.simplify(entry) == 0 for entry in expression)
    return sp.simplify(expression) == 0


@dataclass(frozen=True)
class ClosedRayBerryLedger:
    """Exact local and endpoint data for one normalized closed ray path."""

    coordinate: sp.Symbol
    start: sp.Expr
    end: sp.Expr
    section: sp.ImmutableMatrix
    projector: sp.ImmutableMatrix
    projector_derivative: sp.ImmutableMatrix
    projector_is_constant: bool
    berry_connection: sp.Expr
    endpoint_transition: sp.Expr
    connection_integral: sp.Expr
    bare_integral_phase: sp.Expr
    holonomy: sp.Expr


def closed_ray_berry_ledger(
    section: Any,
    coordinate: sp.Symbol,
    *,
    start: Any = 0,
    end: Any = 2 * sp.pi,
) -> ClosedRayBerryLedger:
    """Return exact Berry data using ``A=i*psi^dagger*d_psi``.

    The section must be exactly normalized and its rank-one projector must
    close between ``start`` and ``end``.  If
    ``psi(end)=tau*psi(start)``, the gauge-invariant answer is
    ``tau*exp(i*integral(A))``.  A bare integral phase is exposed separately
    so callers cannot confuse it with the closed-ray holonomy.
    """

    if not isinstance(coordinate, sp.Symbol):
        raise ValueError("coordinate must be a Symbol")
    path = _column(section)
    lower = _real(start, "start")
    upper = _real(end, "end")
    if sp.simplify(upper - lower).is_positive is not True:
        raise ValueError("end must be explicitly greater than start")

    norm = sp.simplify((path.conjugate().T * path)[0, 0])
    if norm != 1:
        raise ValueError("section must be exactly normalized")

    projector = sp.ImmutableMatrix(sp.simplify(path * path.conjugate().T))
    projector_start = projector.subs(coordinate, lower)
    projector_end = projector.subs(coordinate, upper)
    if not _zero(projector_end - projector_start):
        raise ValueError("section must define a closed rank-one projector path")

    section_start = path.subs(coordinate, lower)
    section_end = path.subs(coordinate, upper)
    transition = sp.simplify(
        (section_start.conjugate().T * section_end)[0, 0]
    )
    if not _zero(section_end - transition * section_start):
        raise ValueError("endpoint section must close by one scalar transition")
    if sp.simplify(sp.conjugate(transition) * transition) != 1:
        raise ValueError("endpoint transition must have unit modulus")

    connection = sp.simplify(
        sp.I * (path.conjugate().T * sp.diff(path, coordinate))[0, 0]
    )
    if sp.simplify(connection - sp.conjugate(connection)) != 0:
        raise ValueError("Berry connection must simplify to an exact real expression")
    connection_integral = sp.simplify(
        sp.integrate(connection, (coordinate, lower, upper))
    )
    if connection_integral.has(sp.Integral):
        raise ValueError("Berry connection integral must evaluate exactly")
    bare_phase = sp.simplify(sp.exp(sp.I * connection_integral))
    derivative = sp.ImmutableMatrix(sp.simplify(sp.diff(projector, coordinate)))
    return ClosedRayBerryLedger(
        coordinate=coordinate,
        start=lower,
        end=upper,
        section=path,
        projector=projector,
        projector_derivative=derivative,
        projector_is_constant=_zero(derivative),
        berry_connection=connection,
        endpoint_transition=transition,
        connection_integral=connection_integral,
        bare_integral_phase=bare_phase,
        holonomy=sp.simplify(transition * bare_phase),
    )


def phase_transform_section(
    section: Any,
    phase: Any,
) -> sp.ImmutableMatrix:
    """Return the normalized-section phase change ``exp(i*chi)*psi``."""

    chi = _real(phase, "phase")
    return sp.ImmutableMatrix(sp.exp(sp.I * chi) * _column(section))


def projective_loop_section(
    winding: Any,
    coordinate: sp.Symbol,
    *,
    periodic: bool = False,
) -> sp.ImmutableMatrix:
    """Return the real lift, or its periodic complex gauge, for winding ``k``."""

    if not isinstance(coordinate, sp.Symbol) or coordinate.is_real is not True:
        raise ValueError("coordinate must be an explicitly real Symbol")
    step = _integer(winding, "winding")
    lift = sp.ImmutableMatrix(
        [sp.cos(step * coordinate / 2), sp.sin(step * coordinate / 2)]
    )
    if periodic:
        return phase_transform_section(lift, -step * coordinate / 2)
    return lift


@dataclass(frozen=True)
class ProjectiveLoopBerryLedger:
    """Two exact gauges of the same integer-winding projective loop."""

    winding: sp.Expr
    parity_character: sp.Expr
    real_lift: ClosedRayBerryLedger
    periodic_section: ClosedRayBerryLedger


def projective_loop_berry_ledger(
    winding: Any,
    coordinate: sp.Symbol,
) -> ProjectiveLoopBerryLedger:
    """Return real and periodic-section ledgers with holonomy ``(-1)**k``."""

    step = _integer(winding, "winding")
    real_lift = closed_ray_berry_ledger(
        projective_loop_section(step, coordinate), coordinate
    )
    periodic_section = closed_ray_berry_ledger(
        projective_loop_section(step, coordinate, periodic=True), coordinate
    )
    return ProjectiveLoopBerryLedger(
        winding=step,
        parity_character=sp.simplify((-1) ** step),
        real_lift=real_lift,
        periodic_section=periodic_section,
    )
