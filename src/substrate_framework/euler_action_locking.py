"""Exact ledgers for classical charge/action locking candidates.

The helpers in this module expose what a reduced nonlinear energy, a compact
phase character, and elementary carrier-complement topology can determine.
They do not construct the missing Euler--Maxwell center manifold, nonlinear
coefficients, compact two-form field, or quantum representation.
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


def leading_nonzero_action(linear_coefficient, quadratic_coefficient):
    """Return ``-nu/(2*a2)`` for ``E=E0+nu*J+a2*J**2+...``.

    A positive result is only a leading candidate.  Establishing a true local
    minimum also requires a controlled analytic remainder and a reduction on
    an actual invariant or constrained same-carrier manifold.
    """

    nu = sp.sympify(linear_coefficient)
    a2 = sp.sympify(quadratic_coefficient)
    _nonzero("quadratic_coefficient", a2)
    return sp.factor(-nu / (2 * a2))


@dataclass(frozen=True)
class ReducedEnergyLocalLedger:
    """Local sign and remainder ledger for a reduced action energy."""

    leading_action: sp.Expr
    leading_second_derivative: sp.Expr
    positive_interior_candidate: bool | None
    zero_is_boundary_selected: bool | None
    remainder_bound_required: bool
    invariant_reduction_required: bool


def reduced_energy_local_ledger(linear_coefficient, quadratic_coefficient):
    """Classify the leading local energy geometry on the half-line ``J>=0``.

    For decidable signs, ``nu<0<a2`` supplies a positive interior candidate,
    while ``nu>0`` makes ``J=0`` the sufficiently local boundary minimum.
    Symbolically undecidable signs are reported as ``None``.
    """

    nu = sp.sympify(linear_coefficient)
    a2 = sp.sympify(quadratic_coefficient)
    _nonzero("quadratic_coefficient", a2)
    positive_candidate = None
    if nu.is_negative is True and a2.is_positive is True:
        positive_candidate = True
    elif nu.is_nonnegative is True or a2.is_nonpositive is True:
        positive_candidate = False
    zero_selected = None
    if nu.is_positive is True:
        zero_selected = True
    elif nu.is_negative is True:
        zero_selected = False
    return ReducedEnergyLocalLedger(
        leading_action=leading_nonzero_action(nu, a2),
        leading_second_derivative=sp.factor(2 * a2),
        positive_interior_candidate=positive_candidate,
        zero_is_boundary_selected=zero_selected,
        remainder_bound_required=True,
        invariant_reduction_required=True,
    )


def compact_phase_charge_per_action(
    charge_unit, character_weight, phase_action_unit
):
    """Return ``kappa_m=g0*m/S0`` for a chosen compact character.

    ``m`` is integral representation data.  ``S0`` is an explicit positive
    action coefficient and remains imported unless a same-carrier law derives
    it.
    """

    g0 = sp.sympify(charge_unit)
    m = sp.sympify(character_weight)
    s0 = sp.sympify(phase_action_unit)
    if m.is_integer is False:
        raise ValueError("character_weight must be an integer")
    _positive("phase_action_unit", s0)
    return sp.factor(g0 * m / s0)


def compact_phase_charge(
    charge_unit, character_weight, phase_number
):
    """Return ``Q=g0*m*N``; the classical number ``N`` is continuous."""

    g0 = sp.sympify(charge_unit)
    m = sp.sympify(character_weight)
    number = sp.sympify(phase_number)
    if m.is_integer is False:
        raise ValueError("character_weight must be an integer")
    return sp.factor(g0 * m * number)


def action_locked_by_fixed_charge(
    total_charge, charge_unit, character_weight, phase_action_unit
):
    """Return the action imposed by ``Q=(g0*m/S0)J``.

    This is a kinematic moment-map relation.  It does not select ``Q``, ``S0``,
    or an energetic minimum.
    """

    kappa_m = compact_phase_charge_per_action(
        charge_unit, character_weight, phase_action_unit
    )
    _nonzero("charge_per_action", kappa_m)
    return sp.factor(sp.sympify(total_charge) / kappa_m)


def electromagnetic_candidate_action_unit(
    charge_unit, permittivity, maxwell_speed, dimensionless_factor=1
):
    """Return ``C*g0**2/(4*pi*epsilon*c)``.

    Dimensional and field-normalization invariance do not determine the
    dimensionless factor ``C``.
    """

    g0 = sp.sympify(charge_unit)
    epsilon = sp.sympify(permittivity)
    speed = sp.sympify(maxwell_speed)
    coefficient = sp.sympify(dimensionless_factor)
    _positive("permittivity", epsilon)
    _positive("maxwell_speed", speed)
    return sp.factor(coefficient * g0**2 / (4 * sp.pi * epsilon * speed))


@dataclass(frozen=True)
class ComplementCohomologyLedger:
    """Integral low-degree cohomology of point and ring complements in R3."""

    point_h1_rank: int
    point_h2_rank: int
    ring_h1_rank: int
    ring_h2_rank: int


def carrier_complement_cohomology() -> ComplementCohomologyLedger:
    """Return the Alexander-duality ranks for point and circle complements."""

    return ComplementCohomologyLedger(
        point_h1_rank=0,
        point_h2_rank=1,
        ring_h1_rank=1,
        ring_h2_rank=0,
    )


@dataclass(frozen=True)
class LabelPullbackBFRoute:
    """Exact local status of ``B=d alpha wedge d beta`` in the BF equation."""

    exterior_derivative_vanishes_away_from_defects: bool
    supplies_smooth_bulk_electric_source: bool
    defect_or_extra_compact_field_required: bool


def smooth_two_label_bf_ledger() -> LabelPullbackBFRoute:
    """Record ``d(d alpha wedge d beta)=0`` away from label defects."""

    return LabelPullbackBFRoute(
        exterior_derivative_vanishes_away_from_defects=True,
        supplies_smooth_bulk_electric_source=False,
        defect_or_extra_compact_field_required=True,
    )


def euler_similarity_action_factor(velocity_scale, inverse_length_scale):
    """Return the exact bare-Euler action factor ``A/B**4``."""

    scale_a = sp.sympify(velocity_scale)
    scale_b = sp.sympify(inverse_length_scale)
    _nonzero("inverse_length_scale", scale_b)
    return sp.factor(scale_a / scale_b**4)


@dataclass(frozen=True)
class LeadingRingScaleBalance:
    """Leading fixed-aspect-ratio fluid/Coulomb ring balance."""

    radius: sp.Expr
    axial_impulse: sp.Expr
    second_derivative: sp.Expr


def leading_ring_scale_balance(
    circulation, charge, fluid_density, permittivity, log_factor
):
    """Minimize the leading thin-ring energy ``L*(A*R+B/R)``.

    The coefficients are ``A=rho*kappa**2/2`` and
    ``B=Q**2/(8*pi**2*epsilon)``. Inputs are positive, selecting the
    positive-radius branch. The result is a fixed-aspect-ratio asymptotic
    candidate rather than a theorem for the full charged Cao energy.
    """

    kappa = sp.sympify(circulation)
    charge_q = sp.sympify(charge)
    rho = sp.sympify(fluid_density)
    epsilon = sp.sympify(permittivity)
    logarithm = sp.sympify(log_factor)
    for name, value in (
        ("circulation", kappa),
        ("charge", charge_q),
        ("fluid_density", rho),
        ("permittivity", epsilon),
        ("log_factor", logarithm),
    ):
        _positive(name, value)
    radius = sp.factor(
        charge_q / (2 * sp.pi * kappa * sp.sqrt(rho * epsilon))
    )
    coefficient_b = charge_q**2 / (8 * sp.pi**2 * epsilon)
    second = sp.factor(2 * logarithm * coefficient_b / radius**3)
    impulse = sp.factor(sp.pi * rho * kappa * radius**2)
    return LeadingRingScaleBalance(
        radius=radius,
        axial_impulse=impulse,
        second_derivative=second,
    )
