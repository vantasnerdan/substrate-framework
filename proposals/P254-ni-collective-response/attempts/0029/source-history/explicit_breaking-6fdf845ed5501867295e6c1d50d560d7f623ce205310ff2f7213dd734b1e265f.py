"""Exact evidence for declared periodic explicit-breaking coordinates.

The functions in this module keep potential curvature, scalar kinetic
normalization, and SU(2) trace conventions in one calculation.  They do not
derive a chiral sector, a GMOR relation, a physical pion dictionary, or any
parameter value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .symmetry_breaking import leading_exponential_kinetic_metric


def _positive(expression: Any, *, name: str) -> sp.Expr:
    value = sp.sympify(expression)
    if value.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return value


@dataclass(frozen=True)
class PeriodicPotentialEvidence:
    """Exact local and quadratic data for ``A*(1-cos(q*field/F))``."""

    field: sp.Symbol
    amplitude: sp.Expr
    coordinate_scale: sp.Expr
    angle_multiplier: sp.Expr
    kinetic_coefficient: sp.Expr
    potential: sp.Expr
    period: sp.Expr
    value_at_origin: sp.Expr
    slope_at_origin: sp.Expr
    curvature_at_origin: sp.Expr
    fourth_derivative_at_origin: sp.Expr
    sixth_order_series: sp.Expr
    generalized_mass_squared: sp.Expr


def periodic_potential_evidence(
    field: sp.Symbol,
    amplitude: Any,
    coordinate_scale: Any,
    kinetic_coefficient: Any,
    *,
    angle_multiplier: Any = 1,
) -> PeriodicPotentialEvidence:
    """Derive a periodic potential and its mass relative to supplied kinetics.

    The scalar quadratic convention is
    ``L_kinetic=kinetic_coefficient*(d field)**2/2``.  The generalized mass
    squared is therefore ``V''(0)/kinetic_coefficient``.  The amplitude is
    deliberately not required to be positive so sign mutations remain
    observable, while the coordinate scale, angle multiplier, and kinetic
    coefficient must be provably positive.
    """

    if not isinstance(field, sp.Symbol):
        raise ValueError("field must be a SymPy symbol")
    scale = _positive(coordinate_scale, name="coordinate scale")
    multiplier = _positive(angle_multiplier, name="angle multiplier")
    kinetic = _positive(kinetic_coefficient, name="kinetic coefficient")
    amplitude_expression = sp.sympify(amplitude)
    angle = multiplier * field / scale
    potential = sp.simplify(amplitude_expression * (1 - sp.cos(angle)))
    curvature = sp.simplify(sp.diff(potential, field, 2).subs(field, 0))
    return PeriodicPotentialEvidence(
        field=field,
        amplitude=amplitude_expression,
        coordinate_scale=scale,
        angle_multiplier=multiplier,
        kinetic_coefficient=kinetic,
        potential=potential,
        period=sp.simplify(2 * sp.pi * scale / multiplier),
        value_at_origin=sp.simplify(potential.subs(field, 0)),
        slope_at_origin=sp.simplify(sp.diff(potential, field).subs(field, 0)),
        curvature_at_origin=curvature,
        fourth_derivative_at_origin=sp.simplify(
            sp.diff(potential, field, 4).subs(field, 0)
        ),
        sixth_order_series=sp.expand(sp.series(potential, field, 0, 8).removeO()),
        generalized_mass_squared=sp.simplify(curvature / kinetic),
    )


@dataclass(frozen=True)
class LocalCurvatureNonuniquenessEvidence:
    """Two globally distinct potentials with the same local quadratic data."""

    field: sp.Symbol
    curvature: sp.Expr
    coordinate_scale: sp.Expr
    periodic_potential: sp.Expr
    quadratic_potential: sp.Expr
    periodic_shift_residual: sp.Expr
    quadratic_shift_residual: sp.Expr
    hessian_difference_at_origin: sp.Expr
    fourth_derivative_difference_at_origin: sp.Expr


def matched_local_curvature_potentials(
    field: sp.Symbol,
    curvature: Any,
    coordinate_scale: Any,
) -> LocalCurvatureNonuniquenessEvidence:
    """Construct cosine and quadratic potentials with one matched Hessian.

    The equal curvature proves only local quadratic equivalence.  Periodicity
    and fourth derivatives distinguish the global functions exactly.
    """

    if not isinstance(field, sp.Symbol):
        raise ValueError("field must be a SymPy symbol")
    scale = _positive(coordinate_scale, name="coordinate scale")
    curvature_expression = sp.sympify(curvature)
    periodic = sp.simplify(
        curvature_expression * scale**2 * (1 - sp.cos(field / scale))
    )
    quadratic = sp.simplify(curvature_expression * field**2 / 2)
    period = 2 * sp.pi * scale
    return LocalCurvatureNonuniquenessEvidence(
        field=field,
        curvature=curvature_expression,
        coordinate_scale=scale,
        periodic_potential=periodic,
        quadratic_potential=quadratic,
        periodic_shift_residual=sp.trigsimp(
            periodic.subs(field, field + period) - periodic
        ),
        quadratic_shift_residual=sp.expand(
            quadratic.subs(field, field + period) - quadratic
        ),
        hessian_difference_at_origin=sp.simplify(
            (sp.diff(periodic, field, 2) - sp.diff(quadratic, field, 2)).subs(field, 0)
        ),
        fourth_derivative_difference_at_origin=sp.simplify(
            (sp.diff(periodic, field, 4) - sp.diff(quadratic, field, 4)).subs(field, 0)
        ),
    )


@dataclass(frozen=True)
class SU2TraceBreakingEvidence:
    """Kinetic and trace-potential data in one SU(2) coordinate convention."""

    field: sp.Symbol
    coordinate_scale: sp.Expr
    angle_multiplier: sp.Expr
    kinetic_prefactor: sp.Expr
    lagrangian_trace_prefactor: sp.Expr
    group_element: sp.ImmutableMatrix
    trace_u_minus_identity: sp.Expr
    lagrangian_breaking_term: sp.Expr
    potential: sp.Expr
    kinetic_coefficient: sp.Expr
    potential_curvature: sp.Expr
    generalized_mass_squared: sp.Expr
    generalized_mass_coordinate_residual: sp.Expr


def su2_trace_breaking_evidence(
    field: sp.Symbol,
    coordinate_scale: Any,
    angle_multiplier: Any,
    kinetic_prefactor: Any,
    lagrangian_trace_prefactor: Any,
) -> SU2TraceBreakingEvidence:
    """Derive the generalized mass from a matched SU(2) kinetic/trace pair.

    The declared coordinate is ``U=exp(i*q*tau3*field/F)`` and the action
    terms are ``Z*Tr(dU*dU.H) + C*Tr(U-I)``.  With potential ``-C*Tr(U-I)``,
    the scalar kinetic coefficient is ``4*Z*q**2/F**2`` and the curvature is
    ``2*C*q**2/F**2``.  Their ratio is ``C/(2*Z)``, independent of the field
    coordinate multiplier.  This covariance does not select ``Z`` or ``C``.
    """

    if not isinstance(field, sp.Symbol):
        raise ValueError("field must be a SymPy symbol")
    scale = _positive(coordinate_scale, name="coordinate scale")
    multiplier = _positive(angle_multiplier, name="angle multiplier")
    kinetic_prefactor_expression = _positive(
        kinetic_prefactor, name="kinetic prefactor"
    )
    trace_prefactor_expression = sp.sympify(lagrangian_trace_prefactor)
    angle = sp.simplify(multiplier * field / scale)
    group_element = sp.ImmutableMatrix(
        [[sp.exp(sp.I * angle), 0], [0, sp.exp(-sp.I * angle)]]
    )
    trace_shift = sp.simplify(sp.trace(group_element - sp.eye(2)).rewrite(sp.cos))
    lagrangian_term = sp.simplify(trace_prefactor_expression * trace_shift)
    potential = sp.simplify(-lagrangian_term)
    tau3 = sp.ImmutableMatrix([[1, 0], [0, -1]])
    kinetic = leading_exponential_kinetic_metric(
        (multiplier * tau3,),
        scale,
        kinetic_prefactor_expression,
    ).kinetic_metric[0, 0]
    curvature = sp.simplify(sp.diff(potential, field, 2).subs(field, 0))
    generalized_mass = sp.simplify(curvature / kinetic)
    coordinate_free_mass = sp.simplify(
        trace_prefactor_expression / (2 * kinetic_prefactor_expression)
    )
    return SU2TraceBreakingEvidence(
        field=field,
        coordinate_scale=scale,
        angle_multiplier=multiplier,
        kinetic_prefactor=kinetic_prefactor_expression,
        lagrangian_trace_prefactor=trace_prefactor_expression,
        group_element=group_element,
        trace_u_minus_identity=trace_shift,
        lagrangian_breaking_term=lagrangian_term,
        potential=potential,
        kinetic_coefficient=sp.simplify(kinetic),
        potential_curvature=curvature,
        generalized_mass_squared=generalized_mass,
        generalized_mass_coordinate_residual=sp.simplify(
            generalized_mass - coordinate_free_mass
        ),
    )


@dataclass(frozen=True)
class ConditionalGMOREvidence:
    """Algebraic consequences of a separately declared GMOR convention."""

    quark_mass_sum: sp.Expr
    condensate: sp.Expr
    decay_scale: sp.Expr
    convention_factor: sp.Expr
    mass_squared: sp.Expr
    relation_residual: sp.Expr
    quark_mass_log_exponent: sp.Expr
    condensate_log_exponent: sp.Expr
    decay_scale_log_exponent: sp.Expr
    convention_factor_log_exponent: sp.Expr
    zero_quark_mass_limit: sp.Expr
    scale_condensate_degeneracy_residual: sp.Expr


def conditional_gmor_evidence(
    quark_mass_sum: Any,
    condensate: Any,
    decay_scale: Any,
    *,
    convention_factor: Any = 1,
) -> ConditionalGMOREvidence:
    """Solve a declared ``M^2 F^2=-c*m_q*condensate`` relation.

    ``c`` keeps condensate and decay-constant conventions explicit.  The
    relation is an input: this helper proves its algebraic scaling and exposes
    a continuous rescaling degeneracy, but it does not derive GMOR from a QFT
    action or determine any of its parameters.
    """

    quark_mass = _positive(quark_mass_sum, name="quark mass sum")
    scale = _positive(decay_scale, name="decay scale")
    factor = _positive(convention_factor, name="convention factor")
    condensate_expression = sp.sympify(condensate)
    mass_squared = sp.simplify(-factor * quark_mass * condensate_expression / scale**2)

    mass_slot, condensate_slot, scale_slot, factor_slot = sp.symbols(
        "m_slot Sigma_slot F_slot c_slot", nonzero=True
    )
    generic_mass = -factor_slot * mass_slot * condensate_slot / scale_slot**2

    def generic_log_exponent(variable: sp.Symbol) -> sp.Expr:
        return sp.simplify(variable * sp.diff(generic_mass, variable) / generic_mass)

    rescaling = sp.Symbol("rho", positive=True)
    rescaled_mass = sp.simplify(
        mass_squared.subs(
            {
                scale: rescaling * scale,
                condensate_expression: rescaling**2 * condensate_expression,
            },
            simultaneous=True,
        )
    )
    return ConditionalGMOREvidence(
        quark_mass_sum=quark_mass,
        condensate=condensate_expression,
        decay_scale=scale,
        convention_factor=factor,
        mass_squared=mass_squared,
        relation_residual=sp.simplify(
            mass_squared * scale**2 + factor * quark_mass * condensate_expression
        ),
        quark_mass_log_exponent=generic_log_exponent(mass_slot),
        condensate_log_exponent=generic_log_exponent(condensate_slot),
        decay_scale_log_exponent=generic_log_exponent(scale_slot),
        convention_factor_log_exponent=generic_log_exponent(factor_slot),
        zero_quark_mass_limit=sp.simplify(mass_squared.subs(quark_mass, 0)),
        scale_condensate_degeneracy_residual=sp.simplify(rescaled_mass - mass_squared),
    )
