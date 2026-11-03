"""Exact conditional exponential-fixation and intensity-coordinate ledgers.

The functions here analyze a supplied two-boundary exponential family.  They
do not derive a Wright--Fisher or Moran process, quantum state space, Born
postulate, measurement or actualization rule, microscopic noise law, physical
amplitude normalization, observed deviation, action quantum, cutoff, or
granularity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _unit_interval(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_number and (
        bool(expression < 0) or bool(expression > 1)
    ):
        raise ValueError(f"{name} must lie in the closed unit interval")
    return expression


@dataclass(frozen=True)
class ExponentialFixationLedger:
    """Exact nonzero-selection branch and its continuous neutral extension."""

    initial_frequency: sp.Expr
    selection_ratio: sp.Expr
    probability: sp.Expr
    neutral_limit: sp.Expr
    continuous_probability: sp.Expr
    boundary_at_zero: sp.Expr
    boundary_at_one: sp.Expr
    frequency_derivative: sp.Expr
    selection_derivative: sp.Expr
    selection_derivative_convexity_gap: sp.Expr
    selection_derivative_factorization_residual: sp.Expr
    complement_probability: sp.Expr
    complement_symmetry_residual: sp.Expr
    small_selection_series_through_cubic: sp.Expr
    bvp_boundary_matrix: sp.ImmutableMatrix
    bvp_boundary_determinant: sp.Expr
    bvp_constant_offset: sp.Expr
    bvp_exponential_coefficient: sp.Expr
    bvp_residual: sp.Expr


@dataclass(frozen=True)
class IntensitySelectionLedger:
    """Normalized-frequency and raw-contrast coordinates for two intensities."""

    first_intensity: sp.Expr
    second_intensity: sp.Expr
    total_intensity: sp.Expr
    initial_frequency: sp.Expr
    normalized_contrast: sp.Expr
    raw_contrast: sp.Expr
    selection_coefficient: sp.Expr
    selection_ratio: sp.Expr
    amplitude_scale: sp.Expr
    scaled_first_intensity: sp.Expr
    scaled_second_intensity: sp.Expr
    scaled_total_intensity: sp.Expr
    scaled_initial_frequency: sp.Expr
    scaled_raw_contrast: sp.Expr
    fixed_coefficient_scaled_selection: sp.Expr
    covariant_scaled_coefficient: sp.Expr
    covariant_scaled_selection: sp.Expr
    unit_normalized_first_intensity: sp.Expr
    unit_normalized_second_intensity: sp.Expr
    unit_normalized_contrast: sp.Expr
    unit_normalized_coefficient: sp.Expr
    unit_normalized_selection: sp.Expr


def continuous_exponential_fixation_probability(
    initial_frequency: Any,
    selection_ratio: Any,
) -> sp.Expr:
    r"""Return the continuous extension of the supplied fixation expression.

    For ``S != 0`` the expression is
    ``(1-exp(-S*x))/(1-exp(-S))``. At ``S=0`` its unique continuous value is
    ``x``. Numeric ``x`` is checked to lie in ``[0,1]``; a symbolic real ``x``
    carries that domain as a caller obligation.
    """

    frequency = _unit_interval(initial_frequency, "initial_frequency")
    selection = _exact_real(selection_ratio, "selection_ratio")
    if selection.is_zero is True:
        return frequency
    branch = sp.simplify(
        (1 - sp.exp(-selection * frequency))
        / (1 - sp.exp(-selection))
    )
    if selection.is_zero is False:
        return branch
    return sp.Piecewise(
        (frequency, sp.Eq(selection, 0)),
        (branch, True),
    )


def exponential_fixation_ledger(
    initial_frequency: Any,
    selection_ratio: Any,
) -> ExponentialFixationLedger:
    r"""Return exact diagnostics for a separately supplied nonzero ``S``.

    The nonzero branch uniquely solves ``u_xx+S*u_x=0`` with
    ``u(0)=0,u(1)=1`` because the two-constant boundary matrix has determinant
    ``exp(-S)-1 != 0``.  Calling this equation a diffusion generator, and
    fixing the physical convention absorbed into ``S``, require separate
    provenance.

    For ``0<x<1``, strict convexity of ``exp`` makes the returned gap
    ``(1-x)+x*exp(S)-exp(S*x)`` positive for every real nonzero ``S``. Since
    the remaining derivative factor is positive, the family is strictly
    increasing in ``S``. This inequality interpretation retains the stated
    open-domain obligation rather than encoding it as a sampled boolean.
    """

    frequency = _unit_interval(initial_frequency, "initial_frequency")
    selection = _exact_real(selection_ratio, "selection_ratio")
    if selection.is_zero is not False:
        raise ValueError("selection_ratio must be provably nonzero")

    x = sp.Symbol("_fixation_x", real=True)
    s = sp.Symbol("_fixation_S", real=True, nonzero=True)
    branch = (1 - sp.exp(-s * x)) / (1 - sp.exp(-s))
    probability = sp.simplify(branch.subs({x: frequency, s: selection}))
    frequency_derivative = sp.simplify(
        sp.diff(branch, x).subs({x: frequency, s: selection})
    )
    selection_derivative = sp.simplify(
        sp.diff(branch, s).subs({x: frequency, s: selection})
    )
    convexity_gap = sp.simplify(
        (1 - frequency)
        + frequency * sp.exp(selection)
        - sp.exp(selection * frequency)
    )
    derivative_factorized = sp.simplify(
        convexity_gap
        * sp.exp(selection)
        * sp.exp(-selection * frequency)
        / (sp.exp(selection) - 1) ** 2
    )

    boundary_matrix = sp.ImmutableMatrix(
        [[1, 1], [1, sp.exp(-selection)]]
    )
    determinant = sp.simplify(boundary_matrix.det())
    constant_offset = sp.simplify(1 / (1 - sp.exp(-selection)))
    exponential_coefficient = sp.simplify(-constant_offset)
    bvp_candidate = constant_offset + exponential_coefficient * sp.exp(-s * x)
    bvp_residual = sp.simplify(
        (sp.diff(bvp_candidate, x, 2) + s * sp.diff(bvp_candidate, x))
        .subs(s, selection)
        .subs(x, frequency)
    )

    neutral_parameter = sp.Symbol("_neutral_parameter", real=True)
    neutral_branch = (
        1 - sp.exp(-neutral_parameter * frequency)
    ) / (1 - sp.exp(-neutral_parameter))
    neutral_limit = sp.simplify(sp.limit(neutral_branch, neutral_parameter, 0))
    series = sp.series(neutral_branch, neutral_parameter, 0, 4).removeO()
    complement = sp.simplify(
        continuous_exponential_fixation_probability(
            1 - frequency,
            -selection,
        )
    )
    return ExponentialFixationLedger(
        initial_frequency=frequency,
        selection_ratio=selection,
        probability=probability,
        neutral_limit=neutral_limit,
        continuous_probability=continuous_exponential_fixation_probability(
            frequency, selection
        ),
        boundary_at_zero=sp.simplify(branch.subs({x: 0, s: selection})),
        boundary_at_one=sp.simplify(branch.subs({x: 1, s: selection})),
        frequency_derivative=frequency_derivative,
        selection_derivative=selection_derivative,
        selection_derivative_convexity_gap=convexity_gap,
        selection_derivative_factorization_residual=sp.simplify(
            selection_derivative - derivative_factorized
        ),
        complement_probability=complement,
        complement_symmetry_residual=sp.simplify(
            probability + complement - 1
        ),
        small_selection_series_through_cubic=sp.simplify(
            series.subs(neutral_parameter, selection)
        ),
        bvp_boundary_matrix=boundary_matrix,
        bvp_boundary_determinant=determinant,
        bvp_constant_offset=constant_offset,
        bvp_exponential_coefficient=exponential_coefficient,
        bvp_residual=bvp_residual,
    )


def two_intensity_selection_ledger(
    first_intensity: Any,
    second_intensity: Any,
    selection_coefficient: Any,
    amplitude_scale: Any,
) -> IntensitySelectionLedger:
    r"""Return exact raw, normalized, and rescaled two-intensity coordinates.

    Let ``N=I1+I2``, ``x=I1/N``, and ``S=kappa*(I1-I2)``. A common amplitude
    rescaling by ``lambda`` sends both intensities to ``lambda**2`` times their
    old values. It preserves ``x`` but sends ``S`` to ``lambda**2*S`` only if
    ``kappa`` is held fixed. The coordinate change
    ``kappa->kappa/lambda**2`` preserves ``S``. Unit-normalizing the intensities
    similarly requires ``kappa_unit=kappa*N`` to represent the same ``S``.
    Whether raw intensity or ``kappa`` is a physical observable is outside this
    algebra and must be established separately.
    """

    intensity1 = _exact_positive(first_intensity, "first_intensity")
    intensity2 = _exact_positive(second_intensity, "second_intensity")
    coefficient = _exact_real(selection_coefficient, "selection_coefficient")
    scale = _exact_positive(amplitude_scale, "amplitude_scale")
    total = sp.simplify(intensity1 + intensity2)
    frequency = sp.simplify(intensity1 / total)
    normalized_contrast = sp.simplify((intensity1 - intensity2) / total)
    raw_contrast = sp.simplify(intensity1 - intensity2)
    selection = sp.simplify(coefficient * raw_contrast)
    intensity_scale = sp.simplify(scale**2)
    scaled1 = sp.simplify(intensity_scale * intensity1)
    scaled2 = sp.simplify(intensity_scale * intensity2)
    scaled_total = sp.simplify(scaled1 + scaled2)
    scaled_frequency = sp.simplify(scaled1 / scaled_total)
    scaled_contrast = sp.simplify(scaled1 - scaled2)
    covariant_coefficient = sp.simplify(coefficient / intensity_scale)
    unit1 = frequency
    unit2 = sp.simplify(1 - frequency)
    unit_contrast = sp.simplify(unit1 - unit2)
    unit_coefficient = sp.simplify(coefficient * total)
    return IntensitySelectionLedger(
        first_intensity=intensity1,
        second_intensity=intensity2,
        total_intensity=total,
        initial_frequency=frequency,
        normalized_contrast=normalized_contrast,
        raw_contrast=raw_contrast,
        selection_coefficient=coefficient,
        selection_ratio=selection,
        amplitude_scale=scale,
        scaled_first_intensity=scaled1,
        scaled_second_intensity=scaled2,
        scaled_total_intensity=scaled_total,
        scaled_initial_frequency=scaled_frequency,
        scaled_raw_contrast=scaled_contrast,
        fixed_coefficient_scaled_selection=sp.simplify(
            coefficient * scaled_contrast
        ),
        covariant_scaled_coefficient=covariant_coefficient,
        covariant_scaled_selection=sp.simplify(
            covariant_coefficient * scaled_contrast
        ),
        unit_normalized_first_intensity=unit1,
        unit_normalized_second_intensity=unit2,
        unit_normalized_contrast=unit_contrast,
        unit_normalized_coefficient=unit_coefficient,
        unit_normalized_selection=sp.simplify(
            unit_coefficient * unit_contrast
        ),
    )
