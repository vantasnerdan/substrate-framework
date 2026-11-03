"""Exact finite charge traces with explicit Abelian normalization provenance.

The APIs in this module operate on a separately supplied finite state table.
They derive no representation, anomaly condition, gauge action, kinetic
coefficient, unification boundary, or physical mixing angle.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp


@dataclass(frozen=True)
class WeightedChargeState:
    """One exact generator-eigenvalue pair with a positive multiplicity."""

    label: str
    multiplicity: int
    t3: Any
    abelian_charge: Any


@dataclass(frozen=True)
class FiniteChargeTraceLedger:
    """Weighted traces for ``Q=t3+c*y`` over a declared finite table."""

    states: tuple[WeightedChargeState, ...]
    electric_coefficient: sp.Expr
    state_count: int
    trace_t3_squared: sp.Expr
    trace_abelian_squared: sp.Expr
    trace_cross: sp.Expr
    trace_electric_squared: sp.Expr
    expanded_trace_electric_squared: sp.Expr
    decomposition_residual: sp.Expr
    trace_ratio: sp.Expr | None


@dataclass(frozen=True)
class AbelianNormalizationLedger:
    """One Abelian generator-coordinate change and its covariant companions."""

    base: FiniteChargeTraceLedger
    generator_rescaling: sp.Expr
    base_abelian_coupling: sp.Expr
    rescaled_states: tuple[WeightedChargeState, ...]
    rescaled_abelian_coupling: sp.Expr
    rescaled_electric_coefficient: sp.Expr
    fixed_coefficient: FiniteChargeTraceLedger
    covariant: FiniteChargeTraceLedger
    charge_product_residuals: tuple[sp.Expr, ...]
    coupled_trace_norm: sp.Expr
    rescaled_coupled_trace_norm: sp.Expr
    coupled_trace_norm_residual: sp.Expr


@dataclass(frozen=True)
class ChargeCouplingAngleLedger:
    """Compare a supplied two-coupling angle with a supplied trace ratio."""

    su2_trace: sp.Expr
    abelian_trace: sp.Expr
    su2_coupling: sp.Expr
    abelian_coupling: sp.Expr
    coupling_angle: sp.Expr
    trace_angle: sp.Expr
    angle_residual: sp.Expr
    coupling_squared_ratio: sp.Expr
    required_coupling_squared_ratio: sp.Expr
    equality_numerator: sp.Expr
    su2_inverse_trace_coefficient: sp.Expr
    abelian_inverse_trace_coefficient: sp.Expr
    common_coefficient_residual: sp.Expr


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _normalize_states(
    states: Iterable[WeightedChargeState],
) -> tuple[WeightedChargeState, ...]:
    normalized: list[WeightedChargeState] = []
    labels: set[str] = set()
    for index, state in enumerate(states):
        if not isinstance(state, WeightedChargeState):
            raise TypeError("states must contain WeightedChargeState records")
        if not isinstance(state.label, str) or not state.label.strip():
            raise ValueError("state labels must be non-empty strings")
        if state.label in labels:
            raise ValueError("state labels must be unique provenance keys")
        labels.add(state.label)
        multiplicity = sp.sympify(state.multiplicity)
        if (
            multiplicity.is_number is not True
            or multiplicity.is_integer is not True
            or multiplicity.is_positive is not True
        ):
            raise ValueError("state multiplicities must be positive integers")
        normalized.append(
            WeightedChargeState(
                label=state.label,
                multiplicity=int(multiplicity),
                t3=_exact_real(state.t3, f"states[{index}].t3"),
                abelian_charge=_exact_real(
                    state.abelian_charge,
                    f"states[{index}].abelian_charge",
                ),
            )
        )
    if not normalized:
        raise ValueError("at least one weighted state is required")
    return tuple(normalized)


def finite_charge_trace_ledger(
    states: Iterable[WeightedChargeState],
    *,
    electric_coefficient: Any = 1,
) -> FiniteChargeTraceLedger:
    """Return exact weighted traces for a declared finite state table.

    The declared electric coordinate is ``Q=t3+electric_coefficient*y``.
    ``trace_ratio`` is ``Tr(t3^2)/Tr(Q^2)`` when the denominator is provably
    nonzero and otherwise ``None``. No physical angle or coupling relation is
    inferred from that ratio.
    """

    table = _normalize_states(states)
    coefficient = _exact_real(electric_coefficient, "electric_coefficient")

    def weighted_sum(function: Any) -> sp.Expr:
        return sp.simplify(
            sum(
                state.multiplicity
                * function(sp.sympify(state.t3), sp.sympify(state.abelian_charge))
                for state in table
            )
        )

    trace_t3 = weighted_sum(lambda t3, _y: t3**2)
    trace_y = weighted_sum(lambda _t3, y: y**2)
    trace_cross = weighted_sum(lambda t3, y: t3 * y)
    trace_q = weighted_sum(lambda t3, y: (t3 + coefficient * y) ** 2)
    expanded = sp.simplify(
        trace_t3 + 2 * coefficient * trace_cross + coefficient**2 * trace_y
    )
    ratio = (
        sp.simplify(trace_t3 / trace_q)
        if trace_q.is_zero is False
        else None
    )
    return FiniteChargeTraceLedger(
        states=table,
        electric_coefficient=coefficient,
        state_count=sum(state.multiplicity for state in table),
        trace_t3_squared=trace_t3,
        trace_abelian_squared=trace_y,
        trace_cross=trace_cross,
        trace_electric_squared=trace_q,
        expanded_trace_electric_squared=expanded,
        decomposition_residual=sp.simplify(trace_q - expanded),
        trace_ratio=ratio,
    )


def weighted_abelian_moment(
    states: Iterable[WeightedChargeState],
    power: int,
) -> sp.Expr:
    """Return ``sum multiplicity*y**power`` for a nonnegative integer power."""

    exponent = sp.sympify(power)
    if (
        exponent.is_number is not True
        or exponent.is_integer is not True
        or exponent.is_nonnegative is not True
    ):
        raise ValueError("power must be a nonnegative integer")
    table = _normalize_states(states)
    return sp.simplify(
        sum(
            state.multiplicity * sp.sympify(state.abelian_charge) ** int(exponent)
            for state in table
        )
    )


def abelian_normalization_ledger(
    states: Iterable[WeightedChargeState],
    generator_rescaling: Any,
    abelian_coupling: Any,
    *,
    electric_coefficient: Any = 1,
) -> AbelianNormalizationLedger:
    """Rescale ``y'=rho*y`` with ``g'=g/rho`` and ``c'=c/rho``.

    The inverse coupling and electric coefficient transformations preserve
    both ``g*y`` and ``Q=t3+c*y``. Holding ``c`` fixed instead changes the
    electric coordinate and generally changes its trace ratio.
    """

    table = _normalize_states(states)
    rho = _positive_exact(generator_rescaling, "generator_rescaling")
    coupling = _positive_exact(abelian_coupling, "abelian_coupling")
    coefficient = _exact_real(electric_coefficient, "electric_coefficient")
    base = finite_charge_trace_ledger(
        table,
        electric_coefficient=coefficient,
    )
    rescaled_states = tuple(
        WeightedChargeState(
            label=state.label,
            multiplicity=state.multiplicity,
            t3=state.t3,
            abelian_charge=sp.simplify(rho * sp.sympify(state.abelian_charge)),
        )
        for state in table
    )
    rescaled_coupling = sp.simplify(coupling / rho)
    rescaled_coefficient = sp.simplify(coefficient / rho)
    fixed = finite_charge_trace_ledger(
        rescaled_states,
        electric_coefficient=coefficient,
    )
    covariant = finite_charge_trace_ledger(
        rescaled_states,
        electric_coefficient=rescaled_coefficient,
    )
    product_residuals = tuple(
        sp.simplify(
            rescaled_coupling * sp.sympify(rescaled.abelian_charge)
            - coupling * sp.sympify(original.abelian_charge)
        )
        for original, rescaled in zip(table, rescaled_states, strict=True)
    )
    coupled_norm = sp.simplify(coupling**2 * base.trace_abelian_squared)
    rescaled_coupled_norm = sp.simplify(
        rescaled_coupling**2 * covariant.trace_abelian_squared
    )
    return AbelianNormalizationLedger(
        base=base,
        generator_rescaling=rho,
        base_abelian_coupling=coupling,
        rescaled_states=rescaled_states,
        rescaled_abelian_coupling=rescaled_coupling,
        rescaled_electric_coefficient=rescaled_coefficient,
        fixed_coefficient=fixed,
        covariant=covariant,
        charge_product_residuals=product_residuals,
        coupled_trace_norm=coupled_norm,
        rescaled_coupled_trace_norm=rescaled_coupled_norm,
        coupled_trace_norm_residual=sp.simplify(
            rescaled_coupled_norm - coupled_norm
        ),
    )


def charge_coupling_angle_ledger(
    su2_trace: Any,
    abelian_trace: Any,
    su2_coupling: Any,
    abelian_coupling: Any,
) -> ChargeCouplingAngleLedger:
    """Compare a declared two-coupling angle with a declared trace ratio.

    Equality holds exactly when
    ``g_y**2/g_2**2 = su2_trace/abelian_trace``. Equivalently, the separately
    defined inverse-trace kinetic coefficients are equal. The function checks
    the algebra; it does not derive that common-coefficient premise.
    """

    trace_2 = _positive_exact(su2_trace, "su2_trace")
    trace_y = _positive_exact(abelian_trace, "abelian_trace")
    coupling_2 = _positive_exact(su2_coupling, "su2_coupling")
    coupling_y = _positive_exact(abelian_coupling, "abelian_coupling")
    coupling_angle = sp.simplify(
        coupling_y**2 / (coupling_2**2 + coupling_y**2)
    )
    trace_angle = sp.simplify(trace_2 / (trace_2 + trace_y))
    equality_numerator = sp.simplify(
        coupling_y**2 * trace_y - coupling_2**2 * trace_2
    )
    coefficient_2 = sp.simplify(1 / (coupling_2**2 * trace_2))
    coefficient_y = sp.simplify(1 / (coupling_y**2 * trace_y))
    return ChargeCouplingAngleLedger(
        su2_trace=trace_2,
        abelian_trace=trace_y,
        su2_coupling=coupling_2,
        abelian_coupling=coupling_y,
        coupling_angle=coupling_angle,
        trace_angle=trace_angle,
        angle_residual=sp.simplify(coupling_angle - trace_angle),
        coupling_squared_ratio=sp.simplify(coupling_y**2 / coupling_2**2),
        required_coupling_squared_ratio=sp.simplify(trace_2 / trace_y),
        equality_numerator=equality_numerator,
        su2_inverse_trace_coefficient=coefficient_2,
        abelian_inverse_trace_coefficient=coefficient_y,
        common_coefficient_residual=sp.simplify(coefficient_y - coefficient_2),
    )


def common_trace_normalized_coupling_angle(
    su2_trace: Any,
    abelian_trace: Any,
    common_inverse_coefficient: Any,
) -> ChargeCouplingAngleLedger:
    """Instantiate a separately supplied common inverse-trace coefficient.

    The premises are ``1/g_2^2=C*S_2`` and ``1/g_y^2=C*S_y``. They imply the
    trace ratio conditionally, but this helper does not derive ``C`` or justify
    why the same coefficient applies to both sectors.
    """

    trace_2 = _positive_exact(su2_trace, "su2_trace")
    trace_y = _positive_exact(abelian_trace, "abelian_trace")
    coefficient = _positive_exact(
        common_inverse_coefficient,
        "common_inverse_coefficient",
    )
    coupling_2 = sp.sqrt(1 / (coefficient * trace_2))
    coupling_y = sp.sqrt(1 / (coefficient * trace_y))
    return charge_coupling_angle_ledger(
        trace_2,
        trace_y,
        coupling_2,
        coupling_y,
    )
