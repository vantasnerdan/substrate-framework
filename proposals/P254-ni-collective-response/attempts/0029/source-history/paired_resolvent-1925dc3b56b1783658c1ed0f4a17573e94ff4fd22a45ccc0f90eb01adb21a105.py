"""Exact finite paired-resolvent and Schur-complement identities.

These helpers evaluate declared finite complex matrices. A common imaginary
shift is phenomenological input; the functions do not derive a Lindblad
equation, physical loss mechanism, transition probability, or rate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Sequence

import sympy as sp


PairScaling = Literal["fixed_per_pair", "fixed_sum"]


def _real_expression(value: Any, *, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


def _nonzero_real(value: Any, *, name: str) -> sp.Expr:
    expression = _real_expression(value, name=name)
    if expression.is_zero is not False:
        raise ValueError(f"{name} must be explicitly nonzero")
    return expression


def _nonnegative_real(value: Any, *, name: str) -> sp.Expr:
    expression = _real_expression(value, name=name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be explicitly nonnegative")
    return expression


def _positive_integer(value: Any, *, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_positive is not True
    ):
        raise ValueError(f"{name} must be a positive integer")
    return int(expression)


@dataclass(frozen=True)
class SymmetricPairLossLedger:
    """Exact zero-energy loss scales for one symmetric detuning pair."""

    detuning: sp.Expr
    coupling_product: sp.Expr
    small_loss_linear_coefficient: sp.Expr
    large_loss_inverse_coefficient: sp.Expr
    stationary_positive_loss: sp.Expr
    peak_magnitude: sp.Expr


def asymmetric_pair_resolvent(
    detuning: Any,
    loss: Any,
    positive_energy_product: Any,
    negative_energy_product: Any,
    *,
    spectral_energy: Any = 0,
) -> sp.Expr:
    """Return a declared two-state intermediate resolvent contribution.

    Intermediate energies are ``+Delta-i*Gamma/2`` and
    ``-Delta-i*Gamma/2``. The coupling products multiply their respective
    denominators; they may be complex and need not be conjugates. The
    returned expression is an off-diagonal finite-matrix element, not a rate.
    """

    delta = _nonzero_real(detuning, name="detuning")
    gamma = _nonnegative_real(loss, name="loss")
    energy = _real_expression(spectral_energy, name="spectral_energy")
    positive_product = sp.sympify(positive_energy_product)
    negative_product = sp.sympify(negative_energy_product)
    half_width = sp.I * gamma / 2
    return sp.factor(
        positive_product / (energy - delta + half_width)
        + negative_product / (energy + delta + half_width)
    )


def symmetric_pair_resolvent(
    detuning: Any,
    loss: Any,
    coupling_product: Any = 1,
    *,
    spectral_energy: Any = 0,
) -> sp.Expr:
    """Specialize :func:`asymmetric_pair_resolvent` to equal products."""

    return asymmetric_pair_resolvent(
        detuning,
        loss,
        coupling_product,
        coupling_product,
        spectral_energy=spectral_energy,
    )


def symmetric_pair_loss_ledger(
    detuning: Any,
    coupling_product: Any = 1,
) -> SymmetricPairLossLedger:
    """Return exact zero-energy small/large-loss and one-pair peak data.

    For real positive ``Gamma``, the zero-energy contribution is
    ``-i*c*Gamma/(Delta**2+Gamma**2/4)``. Its magnitude has a unique positive
    stationary point at ``Gamma=2*abs(Delta)`` when ``c`` is nonzero.
    """

    delta = _nonzero_real(detuning, name="detuning")
    product = sp.sympify(coupling_product)
    absolute_delta = sp.Abs(delta)
    return SymmetricPairLossLedger(
        detuning=delta,
        coupling_product=product,
        small_loss_linear_coefficient=sp.simplify(-sp.I * product / delta**2),
        large_loss_inverse_coefficient=sp.simplify(-4 * sp.I * product),
        stationary_positive_loss=2 * absolute_delta,
        peak_magnitude=sp.simplify(sp.Abs(product) / absolute_delta),
    )


def finite_resolvent_effective_block(
    endpoint_block: Sequence[Sequence[Any]] | sp.MatrixBase,
    endpoint_to_intermediate: Sequence[Sequence[Any]] | sp.MatrixBase,
    intermediate_block: Sequence[Sequence[Any]] | sp.MatrixBase,
    intermediate_to_endpoint: Sequence[Sequence[Any]] | sp.MatrixBase,
    *,
    spectral_energy: Any = 0,
) -> sp.ImmutableMatrix:
    """Return ``H_PP+H_PQ*(E*I-H_QQ)^-1*H_QP`` exactly.

    The caller supplies both coupling orientations because a non-Hermitian or
    complex-symmetric model must not silently substitute an adjoint. Matrix
    invertibility at the declared spectral energy is required by SymPy's exact
    inverse operation.
    """

    endpoint = sp.ImmutableMatrix(endpoint_block)
    to_intermediate = sp.ImmutableMatrix(endpoint_to_intermediate)
    intermediate = sp.ImmutableMatrix(intermediate_block)
    to_endpoint = sp.ImmutableMatrix(intermediate_to_endpoint)
    if endpoint.rows != endpoint.cols or intermediate.rows != intermediate.cols:
        raise ValueError("endpoint and intermediate blocks must be square")
    if to_intermediate.shape != (endpoint.rows, intermediate.rows):
        raise ValueError("endpoint_to_intermediate has incompatible shape")
    if to_endpoint.shape != (intermediate.rows, endpoint.rows):
        raise ValueError("intermediate_to_endpoint has incompatible shape")
    energy = sp.sympify(spectral_energy)
    resolvent = (energy * sp.eye(intermediate.rows) - intermediate).inv()
    return sp.ImmutableMatrix(
        sp.simplify(endpoint + to_intermediate * resolvent * to_endpoint)
    )


def equal_pair_resolvent_sum(
    pair_count: Any,
    detuning: Any,
    loss: Any,
    coupling_product: Any = 1,
    *,
    scaling: PairScaling = "fixed_per_pair",
    spectral_energy: Any = 0,
) -> sp.Expr:
    """Return an equal-pair sum under one declared size convention.

    ``fixed_per_pair`` keeps each pair product fixed, so the sum is extensive.
    ``fixed_sum`` divides each pair product by the pair count, keeping the
    total product weight fixed. Changing pair count is model enlargement, not
    numerical mesh refinement.
    """

    count = _positive_integer(pair_count, name="pair_count")
    product = sp.sympify(coupling_product)
    if scaling == "fixed_per_pair":
        per_pair = product
    elif scaling == "fixed_sum":
        per_pair = product / count
    else:
        raise ValueError("scaling must be 'fixed_per_pair' or 'fixed_sum'")
    return sp.simplify(
        count
        * symmetric_pair_resolvent(
            detuning,
            loss,
            per_pair,
            spectral_energy=spectral_energy,
        )
    )
