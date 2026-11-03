"""Exact provenance diagnostics for dimensionful one-loop scale maps.

The helpers in this module separate four questions that are often conflated:
whether declared primitives span a target dimension, how a formal one-loop
scale depends on its dimensionful reference, whether a desired target can be
reconstructed by choosing that reference, and how a fixed quantity's numeric
coordinate changes with the unit standard.  They establish no physical beta
function, absolute scale, preferred unit system, or empirical identification.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_analysis import monomial_exponents
from .linear_systems import LinearSystemDiagnostics, diagnose_linear_system
from .scale_transmutation import one_loop_inverse_energy_length_ledger


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


@dataclass(frozen=True)
class ScaleDimensionLedger:
    """M,L,T target-span diagnostics for ``c``, ``hbar``, and a length."""

    base_dimensions: tuple[str, ...]
    primitive_names: tuple[str, ...]
    primitive_dimension_matrix: sp.ImmutableMatrix
    length_target: sp.ImmutableMatrix
    without_length: LinearSystemDiagnostics
    with_length_primitive_names: tuple[str, ...]
    with_length_dimension_matrix: sp.ImmutableMatrix
    with_length: LinearSystemDiagnostics
    with_length_exponents: sp.ImmutableMatrix


@dataclass(frozen=True)
class OneLoopScaleProvenanceLedger:
    """Reference dependence of a formal scale and its inverse-energy length."""

    reference_energy: sp.Expr
    coupling_squared: sp.Expr
    beta_coefficient: sp.Expr
    conversion: sp.Expr
    reference_rescaling: sp.Expr
    exponent: sp.Expr
    transmuted_energy: sp.Expr
    inverse_energy_length: sp.Expr
    transmuted_to_reference_energy_ratio: sp.Expr
    rescaled_reference_energy: sp.Expr
    rescaled_transmuted_energy: sp.Expr
    rescaled_inverse_energy_length: sp.Expr
    transmuted_energy_rescaling_ratio: sp.Expr
    inverse_length_rescaling_ratio: sp.Expr


@dataclass(frozen=True)
class UnitCoordinateLedger:
    """Coordinates of one fixed quantity in two rescaled unit standards."""

    quantity: sp.Expr
    unit_standard: sp.Expr
    unit_rescaling: sp.Expr
    coordinate: sp.Expr
    rescaled_unit_standard: sp.Expr
    rescaled_coordinate: sp.Expr
    coordinate_rescaling_ratio: sp.Expr


def speed_action_length_dimension_ledger() -> ScaleDimensionLedger:
    """Diagnose a length target from ``c`` and ``hbar`` in M,L,T order.

    With columns ``c=(0,1,-1)`` and ``hbar=(1,2,-1)``, the pure-length target
    ``(0,1,0)`` is outside the column span.  Appending that same target as a
    primitive named ``a`` makes the solve uniquely ``a**1``.  The latter is a
    supplied target primitive, not a derivation of its magnitude.
    """

    primitives = sp.Matrix([[0, 1], [1, 2], [-1, -1]])
    target = sp.Matrix([0, 1, 0])
    without_length = diagnose_linear_system(primitives, target)
    with_length_matrix = primitives.row_join(target)
    with_length = diagnose_linear_system(with_length_matrix, target)
    exponents = monomial_exponents(with_length_matrix, target)
    return ScaleDimensionLedger(
        base_dimensions=("M", "L", "T"),
        primitive_names=("c", "hbar"),
        primitive_dimension_matrix=sp.ImmutableMatrix(primitives),
        length_target=sp.ImmutableMatrix(target),
        without_length=without_length,
        with_length_primitive_names=("c", "hbar", "a"),
        with_length_dimension_matrix=sp.ImmutableMatrix(with_length_matrix),
        with_length=with_length,
        with_length_exponents=sp.ImmutableMatrix(exponents),
    )


def one_loop_scale_provenance_ledger(
    reference_energy: Any,
    coupling_squared: Any,
    beta_coefficient: Any,
    *,
    conversion: Any,
    reference_rescaling: Any,
) -> OneLoopScaleProvenanceLedger:
    """Return exact fixed-input covariance under ``mu0 -> rho*mu0``.

    The formal relation is ``Lambda=mu0*exp(-X)`` with
    ``X=8*pi**2/(b0*g2)`` and the length is ``conversion/Lambda``.  Holding
    ``g2``, ``b0``, and ``conversion`` fixed while rescaling ``mu0`` by
    ``rho`` rescales ``Lambda`` by ``rho`` and the length by ``1/rho``.
    """

    mu0 = _positive(reference_energy, "reference_energy")
    g2 = _positive(coupling_squared, "coupling_squared")
    b0 = _positive(beta_coefficient, "beta_coefficient")
    conversion_value = _positive(conversion, "conversion")
    rho = _positive(reference_rescaling, "reference_rescaling")

    original = one_loop_inverse_energy_length_ledger(
        mu0,
        g2,
        b0,
        reference_conversion=conversion_value,
        transmuted_conversion=conversion_value,
    )
    rescaled = one_loop_inverse_energy_length_ledger(
        rho * mu0,
        g2,
        b0,
        reference_conversion=conversion_value,
        transmuted_conversion=conversion_value,
    )
    return OneLoopScaleProvenanceLedger(
        reference_energy=mu0,
        coupling_squared=g2,
        beta_coefficient=b0,
        conversion=conversion_value,
        reference_rescaling=rho,
        exponent=original.exponent,
        transmuted_energy=original.transmuted_energy,
        inverse_energy_length=original.transmuted_length,
        transmuted_to_reference_energy_ratio=(
            original.transmuted_to_reference_energy_ratio
        ),
        rescaled_reference_energy=sp.simplify(rho * mu0),
        rescaled_transmuted_energy=rescaled.transmuted_energy,
        rescaled_inverse_energy_length=rescaled.transmuted_length,
        transmuted_energy_rescaling_ratio=sp.simplify(
            rescaled.transmuted_energy / original.transmuted_energy
        ),
        inverse_length_rescaling_ratio=sp.simplify(
            rescaled.transmuted_length / original.transmuted_length
        ),
    )


def reference_energy_for_target_transmuted_energy(
    target_energy: Any,
    coupling_squared: Any,
    beta_coefficient: Any,
) -> sp.Expr:
    """Choose the reference energy that reproduces a supplied target scale."""

    target = _positive(target_energy, "target_energy")
    g2 = _positive(coupling_squared, "coupling_squared")
    b0 = _positive(beta_coefficient, "beta_coefficient")
    exponent = 8 * sp.pi**2 / (b0 * g2)
    return sp.simplify(target * sp.exp(exponent))


def reference_energy_for_target_inverse_length(
    target_length: Any,
    coupling_squared: Any,
    beta_coefficient: Any,
    *,
    conversion: Any,
) -> sp.Expr:
    """Choose the reference energy that reproduces a supplied inverse length."""

    target = _positive(target_length, "target_length")
    g2 = _positive(coupling_squared, "coupling_squared")
    b0 = _positive(beta_coefficient, "beta_coefficient")
    conversion_value = _positive(conversion, "conversion")
    exponent = 8 * sp.pi**2 / (b0 * g2)
    return sp.simplify(conversion_value * sp.exp(exponent) / target)


def unit_coordinate_ledger(
    quantity: Any,
    unit_standard: Any,
    unit_rescaling: Any,
) -> UnitCoordinateLedger:
    """Keep a quantity fixed while rescaling its positive unit standard.

    If ``q=N*u`` and the unit is replaced by ``u'=rho*u``, the numeric
    coordinate becomes ``N'=N/rho``.  This covariance changes neither the
    physical quantity nor its provenance and cannot derive either one.
    """

    quantity_value = _positive(quantity, "quantity")
    unit_value = _positive(unit_standard, "unit_standard")
    rho = _positive(unit_rescaling, "unit_rescaling")
    coordinate = sp.simplify(quantity_value / unit_value)
    rescaled_unit = sp.simplify(rho * unit_value)
    rescaled_coordinate = sp.simplify(quantity_value / rescaled_unit)
    return UnitCoordinateLedger(
        quantity=quantity_value,
        unit_standard=unit_value,
        unit_rescaling=rho,
        coordinate=coordinate,
        rescaled_unit_standard=rescaled_unit,
        rescaled_coordinate=rescaled_coordinate,
        coordinate_rescaling_ratio=sp.simplify(
            rescaled_coordinate / coordinate
        ),
    )
