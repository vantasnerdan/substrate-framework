"""Exact ledgers for compact-gauge charge and carrier-action selection.

These helpers separate a dimensionful electromagnetic action scale from a
law selecting charge or mode action.  They describe conditional classical
Euler--Maxwell extensions and do not quantize a particle.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


def _positive(name: str, value: sp.Expr) -> None:
    if value.is_positive is False or value.is_zero is True:
        raise ValueError(f"{name} must be positive")


def _nonzero(name: str, value: sp.Expr) -> None:
    if value.is_zero is True:
        raise ValueError(f"{name} must be nonzero")


def electromagnetic_action_scale(charge, permittivity, maxwell_speed):
    """Return ``Q**2/(4*pi*epsilon*c)`` with dimensions of action."""

    q = sp.sympify(charge)
    epsilon = sp.sympify(permittivity)
    speed = sp.sympify(maxwell_speed)
    _positive("permittivity", epsilon)
    _positive("maxwell_speed", speed)
    return sp.factor(q**2 / (4 * sp.pi * epsilon * speed))


@dataclass(frozen=True)
class GaugeRescalingLedger:
    """Coefficients and invariant combinations under a field rescaling."""

    permittivity: sp.Expr
    permeability: sp.Expr
    coupling: sp.Expr
    maxwell_speed_squared: sp.Expr
    charge_action_numerator: sp.Expr


def gauge_field_rescaling(permittivity, permeability, coupling, field_scale):
    """Apply ``A' = a*A`` and return the transformed coefficient ledger.

    ``field_scale`` is a positive normalization magnitude here.  A possible
    orientation/sign reversal is tracked separately and is not a charge-sign
    restriction.
    """

    epsilon = sp.sympify(permittivity)
    mu = sp.sympify(permeability)
    g = sp.sympify(coupling)
    scale = sp.sympify(field_scale)
    for name, value in (
        ("permittivity", epsilon),
        ("permeability", mu),
        ("field_scale", scale),
    ):
        _positive(name, value)
    epsilon_new = epsilon / scale**2
    mu_new = mu * scale**2
    g_new = g / scale
    return GaugeRescalingLedger(
        permittivity=sp.factor(epsilon_new),
        permeability=sp.factor(mu_new),
        coupling=sp.factor(g_new),
        maxwell_speed_squared=sp.factor(1 / (epsilon_new * mu_new)),
        charge_action_numerator=sp.factor(g_new**2 / epsilon_new),
    )


def transported_tag_charge(tag_coupling, tag_integral, tag_scale=1):
    """Return continuous tag charge ``Q=g_tag*lambda*C_chi``."""

    return sp.factor(
        sp.sympify(tag_coupling)
        * sp.sympify(tag_scale)
        * sp.sympify(tag_integral)
    )


def compact_phase_moment_charge(charge_unit, phase_number, character_weight=1):
    """Return ``Q=g_0*m*N`` for the classical compact-phase moment map.

    Here ``N`` is the integral of a dimensionless number density.  With
    phase-action coefficient ``S_0``, the cotangent action is
    ``J_phase=S_0*N`` and the charge-per-action coefficient is
    ``kappa_m=g_0*m/S_0``.  The phase coordinate is periodic, while ``N`` and
    ``J_phase`` remain real and continuous.  Choosing integer character ``m``
    is additional representation data and does not discretize either
    classical variable.
    """

    weight = sp.sympify(character_weight)
    if weight.is_integer is False:
        raise ValueError("character_weight must be an integer")
    return sp.factor(
        sp.sympify(charge_unit) * weight * sp.sympify(phase_number)
    )


def gauge_action_ratio(charge, permittivity, maxwell_speed, mode_action):
    """Return ``S_EM(Q)/J``; no equation forces this ratio to one."""

    action = sp.sympify(mode_action)
    _positive("mode_action", action)
    return sp.factor(
        electromagnetic_action_scale(charge, permittivity, maxwell_speed)
        / action
    )


def conditional_self_consistent_action(
    charge_per_action, permittivity, maxwell_speed
):
    """Return the root of imposed ``Q=kappa_m*J``, ``J=S_EM``.

    ``charge_per_action`` has units charge/action and is distinct from the
    transported-tag coupling unless a separate physical map identifies them.
    This is algebra for an explicitly added constraint, not a selection law
    of the classical Euler--Maxwell action.
    """

    kappa_m = sp.sympify(charge_per_action)
    epsilon = sp.sympify(permittivity)
    speed = sp.sympify(maxwell_speed)
    _nonzero("charge_per_action", kappa_m)
    _positive("permittivity", epsilon)
    _positive("maxwell_speed", speed)
    return sp.factor(4 * sp.pi * epsilon * speed / kappa_m**2)


@dataclass(frozen=True)
class CompactU1TopologyLedger:
    """Low-degree topology relevant to smooth U(1) fields in three-space."""

    pi_3_u1: int
    h2_s3: int
    h2_s2: str
    whole_space_magnetic_flux: int


@dataclass(frozen=True)
class BFRouteLedger:
    """Topological data for the candidate compact 3+1-dimensional BF route."""

    euler_vorticity_is_compact_connection: bool
    euler_vorticity_periods_are_integral: bool
    compact_two_form_is_added_structure: bool
    level_is_euler_selected: bool


def compact_u1_topology_ledger() -> CompactU1TopologyLedger:
    """Return ``pi3(U1)=0``, ``H2(S3)=0``, and ``H2(S2)=Z``.

    A smooth global source-free Maxwell potential on all of R3 has zero
    magnetic flux through the sphere at infinity.
    """

    return CompactU1TopologyLedger(
        pi_3_u1=0,
        h2_s3=0,
        h2_s2="Z",
        whole_space_magnetic_flux=0,
    )


def compact_bf_route_ledger() -> BFRouteLedger:
    """Classify the direct Euler-vorticity candidate for compact BF coupling.

    A smooth Euler vorticity two-form is a real advected exact form on
    whole-space ``R**3``.  It is not thereby a compact two-form gauge
    connection, and its continuous circulation/flux data do not supply an
    integral BF level.
    """

    return BFRouteLedger(
        euler_vorticity_is_compact_connection=False,
        euler_vorticity_periods_are_integral=False,
        compact_two_form_is_added_structure=True,
        level_is_euler_selected=False,
    )


def charge_action_map_jacobian(
    tag_coupling,
    tag_integral,
    tag_scale,
    mode_action,
    permittivity,
    maxwell_speed,
):
    """Return the Jacobian of ``(Q,J,S_EM(Q))`` in ``(lambda,J)``.

    The second column is exactly ``(0, 1, 0)**T``.  Consequently fixing the
    classical charge and its electromagnetic action scale cannot fix an
    independently supplied carrier-mode action.
    """

    g_tag = sp.sympify(tag_coupling)
    c_tag = sp.sympify(tag_integral)
    lam = sp.sympify(tag_scale)
    action = sp.sympify(mode_action)
    epsilon = sp.sympify(permittivity)
    speed = sp.sympify(maxwell_speed)
    _positive("permittivity", epsilon)
    _positive("maxwell_speed", speed)
    charge = transported_tag_charge(g_tag, c_tag, lam)
    outputs = sp.Matrix(
        [
            charge,
            action,
            electromagnetic_action_scale(charge, epsilon, speed),
        ]
    )
    return outputs.jacobian(sp.Matrix([lam, action]))
