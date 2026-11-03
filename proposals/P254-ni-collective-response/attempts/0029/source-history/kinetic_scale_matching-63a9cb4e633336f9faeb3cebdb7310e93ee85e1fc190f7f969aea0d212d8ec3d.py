"""Exact composition of inverse-length scale maps and affine kinetic flow.

The ledgers in this module compose :mod:`scale_transmutation` with
:mod:`vacuum_polarization`.  They do not identify either length with a
substrate cutoff or soliton, select a conversion factor or matching boundary,
or turn an inverse kinetic coordinate into a physical gauge coupling.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .scale_transmutation import (
    OneLoopLengthRatioLedger,
    one_loop_inverse_energy_length_ledger,
)
from .vacuum_polarization import (
    MatterInducedKineticEvidence,
    matter_induced_kinetic_evidence,
)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


@dataclass(frozen=True)
class InverseLengthKineticEvidence:
    """An affine kinetic family evaluated on two declared length scales."""

    reference_length: sp.Expr
    evaluation_length: sp.Expr
    reference_conversion: sp.Expr
    evaluation_conversion: sp.Expr
    reference_energy: sp.Expr
    evaluation_energy: sp.Expr
    evaluation_to_reference_length_ratio: sp.Expr
    evaluation_to_reference_conversion_ratio: sp.Expr
    reference_to_evaluation_energy_ratio: sp.Expr
    energy_length_ratio_residual: sp.Expr
    scale_logarithm: sp.Expr
    kinetic: MatterInducedKineticEvidence
    affine_composition_residual: sp.Expr
    common_length_rescaling: sp.Symbol
    rescaled_reference_energy: sp.Expr
    rescaled_evaluation_energy: sp.Expr
    common_rescaling_log_residual: sp.Expr
    physical_scale_identification_is_separate_premise: bool
    positive_inverse_interpretation_is_separate_premise: bool


def inverse_length_scale_kinetic_evidence(
    reference_length: Any,
    evaluation_length: Any,
    reference_conversion: Any,
    evaluation_conversion: Any,
    renormalized_local_coefficient: Any,
    finite_matching_offset: Any,
    scalar_weight: Any,
    dirac_weight: Any,
) -> InverseLengthKineticEvidence:
    r"""Compose ``E=K/ell`` with the complete affine kinetic family.

    For independently supplied positive lengths and conversions, define
    ``E0=K0/ell0`` and ``E1=K1/ell1``.  If
    ``R_ell=ell1/ell0`` and ``R_K=K1/K0``, then exactly

    ``E0/E1=R_ell/R_K``

    and the kinetic logarithm is ``log(R_ell/R_K)``.  Thus unequal
    conversions cannot be discarded while the lengths are held fixed.  The
    returned kinetic object retains the independent affine reference value
    from :func:`matter_induced_kinetic_evidence`.
    """

    length0 = _positive_exact(reference_length, "reference length")
    length1 = _positive_exact(evaluation_length, "evaluation length")
    conversion0 = _positive_exact(reference_conversion, "reference conversion")
    conversion1 = _positive_exact(
        evaluation_conversion,
        "evaluation conversion",
    )
    energy0 = sp.simplify(conversion0 / length0)
    energy1 = sp.simplify(conversion1 / length1)
    length_ratio = sp.simplify(length1 / length0)
    conversion_ratio = sp.simplify(conversion1 / conversion0)
    energy_ratio = sp.simplify(energy0 / energy1)
    ratio_residual = sp.simplify(energy_ratio - length_ratio / conversion_ratio)
    scale_logarithm = sp.simplify(sp.log(energy_ratio))

    kinetic = matter_induced_kinetic_evidence(
        energy1,
        energy0,
        renormalized_local_coefficient,
        finite_matching_offset,
        scalar_weight,
        dirac_weight,
    )
    composition_residual = sp.simplify(
        kinetic.kinetic_coefficient
        - kinetic.reference_value
        - kinetic.running_coefficient * scale_logarithm
    )

    rho = sp.Symbol("rho_common_length", positive=True)
    rescaled_energy0 = sp.simplify(conversion0 / (rho * length0))
    rescaled_energy1 = sp.simplify(conversion1 / (rho * length1))
    rescaled_logarithm = sp.simplify(
        sp.log(sp.simplify(rescaled_energy0 / rescaled_energy1))
    )

    return InverseLengthKineticEvidence(
        reference_length=length0,
        evaluation_length=length1,
        reference_conversion=conversion0,
        evaluation_conversion=conversion1,
        reference_energy=energy0,
        evaluation_energy=energy1,
        evaluation_to_reference_length_ratio=length_ratio,
        evaluation_to_reference_conversion_ratio=conversion_ratio,
        reference_to_evaluation_energy_ratio=energy_ratio,
        energy_length_ratio_residual=ratio_residual,
        scale_logarithm=scale_logarithm,
        kinetic=kinetic,
        affine_composition_residual=composition_residual,
        common_length_rescaling=rho,
        rescaled_reference_energy=rescaled_energy0,
        rescaled_evaluation_energy=rescaled_energy1,
        common_rescaling_log_residual=sp.simplify(
            rescaled_logarithm - scale_logarithm
        ),
        physical_scale_identification_is_separate_premise=True,
        positive_inverse_interpretation_is_separate_premise=True,
    )


@dataclass(frozen=True)
class OneLoopScaleMatchedKineticEvidence:
    """A consistent one-loop length map composed with affine kinetic flow."""

    transmutation: OneLoopLengthRatioLedger
    matched: InverseLengthKineticEvidence
    exponent: sp.Expr
    scale_log_cancellation_residual: sp.Expr
    general_kinetic_coefficient: sp.Expr
    expected_general_kinetic_coefficient: sp.Expr
    general_kinetic_residual: sp.Expr
    zero_matching_kinetic_coefficient: sp.Expr
    expected_zero_matching_kinetic_coefficient: sp.Expr
    zero_matching_residual: sp.Expr
    zero_matching_inverse_kinetic_coordinate: sp.Expr
    expected_zero_matching_inverse_kinetic_coordinate: sp.Expr
    zero_matching_inverse_residual: sp.Expr
    conversion_factors_need_not_be_equal: bool
    zero_matching_is_separate_premise: bool
    physical_coupling_interpretation_is_separate_premise: bool


def one_loop_scale_matched_kinetic_evidence(
    reference_energy: Any,
    coupling_squared: Any,
    beta_coefficient: Any,
    *,
    reference_conversion: Any,
    transmuted_conversion: Any,
    renormalized_local_coefficient: Any,
    finite_matching_offset: Any,
    scalar_weight: Any,
    dirac_weight: Any,
) -> OneLoopScaleMatchedKineticEvidence:
    r"""Compose the accepted formal one-loop hierarchy with kinetic flow.

    Let ``X=8*pi**2/(b0*g2)`` and ``E1=E0*exp(-X)``.  The paired inverse-
    energy lengths obey
    ``ell1/ell0=(K1/K0)*exp(X)``.  Substitution into the general scale map
    cancels ``K1/K0`` without requiring ``K1=K0`` and gives

    ``Z(E1)=Z_ref+b/(b0*g2)``.

    This API additionally requires a provably positive total matter
    coefficient so the separately declared zero-matching branch has the
    positive inverse coordinate ``b0*g2/b``.  It does not infer that branch or
    authorize a physical gauge-coupling interpretation.
    """

    energy0 = _positive_exact(reference_energy, "reference energy")
    g2 = _positive_exact(coupling_squared, "coupling squared")
    b0 = _positive_exact(beta_coefficient, "beta coefficient")
    conversion0 = _positive_exact(reference_conversion, "reference conversion")
    conversion1 = _positive_exact(
        transmuted_conversion,
        "transmuted conversion",
    )
    transmutation = one_loop_inverse_energy_length_ledger(
        energy0,
        g2,
        b0,
        reference_conversion=conversion0,
        transmuted_conversion=conversion1,
    )
    matched = inverse_length_scale_kinetic_evidence(
        transmutation.reference_length,
        transmutation.transmuted_length,
        conversion0,
        conversion1,
        renormalized_local_coefficient,
        finite_matching_offset,
        scalar_weight,
        dirac_weight,
    )
    coefficient = matched.kinetic.one_loop_coefficient
    if coefficient.is_positive is not True:
        raise ValueError(
            "total matter coefficient must be provably positive for the "
            "zero-matching inverse coordinate"
        )

    exponent = sp.simplify(8 * sp.pi**2 / (b0 * g2))
    expected_general = sp.simplify(
        matched.kinetic.reference_value + coefficient / (b0 * g2)
    )
    zero_expected = sp.simplify(coefficient / (b0 * g2))
    inverse_expected = sp.simplify(b0 * g2 / coefficient)
    zero_actual = matched.kinetic.zero_matching_kinetic_coefficient
    inverse_actual = sp.simplify(1 / zero_actual)

    return OneLoopScaleMatchedKineticEvidence(
        transmutation=transmutation,
        matched=matched,
        exponent=exponent,
        scale_log_cancellation_residual=sp.simplify(
            matched.scale_logarithm - exponent
        ),
        general_kinetic_coefficient=matched.kinetic.kinetic_coefficient,
        expected_general_kinetic_coefficient=expected_general,
        general_kinetic_residual=sp.simplify(
            matched.kinetic.kinetic_coefficient - expected_general
        ),
        zero_matching_kinetic_coefficient=zero_actual,
        expected_zero_matching_kinetic_coefficient=zero_expected,
        zero_matching_residual=sp.simplify(zero_actual - zero_expected),
        zero_matching_inverse_kinetic_coordinate=inverse_actual,
        expected_zero_matching_inverse_kinetic_coordinate=inverse_expected,
        zero_matching_inverse_residual=sp.simplify(
            inverse_actual - inverse_expected
        ),
        conversion_factors_need_not_be_equal=True,
        zero_matching_is_separate_premise=True,
        physical_coupling_interpretation_is_separate_premise=True,
    )
