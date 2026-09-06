"""Exact probability algebra for a declared first-event detector hypothesis.

The functions do not derive stochastic clocks, detector reset, or quantum
measurement from Euler.  They expose what follows once those ingredients are
supplied and keep physical channel intensities explicit.
"""

from dataclasses import dataclass
from typing import Sequence

import sympy as sp


@dataclass(frozen=True)
class FirstEventLaw:
    """Competing-clock probabilities and their total event rate."""

    probabilities: tuple[sp.Expr, ...]
    total_rate: sp.Expr


def _is_zero_matrix(matrix):
    return all(sp.trigsimp(sp.simplify(value)) == 0 for value in matrix)


def projective_channel_intensities(state, analyzer):
    """Return ``|U psi|^2`` for an exactly normalized state and unitary ``U``."""
    psi = sp.Matrix(state)
    unitary = sp.Matrix(analyzer)
    if psi.cols != 1 or psi.rows == 0:
        raise ValueError("state must be a nonempty column")
    if unitary.shape != (psi.rows, psi.rows):
        raise ValueError("analyzer must be square on the state space")
    identity = sp.eye(psi.rows)
    if not _is_zero_matrix(unitary.H * unitary - identity):
        raise ValueError("analyzer must be exactly unitary")
    norm = sp.simplify((psi.H * psi)[0])
    if norm != 1:
        raise ValueError("state must be exactly normalized")
    output = unitary * psi
    return tuple(sp.simplify(sp.conjugate(value) * value) for value in output)


def first_event_probabilities(intensities: Sequence, rate_scale):
    """Return the winner law for independent exponential clocks ``kappa*I_i``."""
    values = tuple(sp.sympify(value) for value in intensities)
    scale = sp.sympify(rate_scale)
    if not values:
        raise ValueError("at least one intensity is required")
    if scale.is_positive is not True:
        raise ValueError("rate_scale must be positive")
    if any(value.is_nonnegative is not True for value in values):
        raise ValueError("intensities must be explicitly nonnegative")
    total_intensity = sp.simplify(sum(values))
    if total_intensity.is_positive is not True:
        raise ValueError("total intensity must be positive")
    return FirstEventLaw(
        probabilities=tuple(sp.simplify(value / total_intensity) for value in values),
        total_rate=sp.simplify(scale * total_intensity),
    )


def post_reset_transition_probabilities(analyzer):
    """Return ``|U_ji|^2`` after an assumed reset to the prior outcome basis."""
    unitary = sp.Matrix(analyzer)
    if unitary.rows == 0 or unitary.rows != unitary.cols:
        raise ValueError("analyzer must be a nonempty square matrix")
    if not _is_zero_matrix(unitary.H * unitary - sp.eye(unitary.rows)):
        raise ValueError("analyzer must be exactly unitary")
    return sp.ImmutableMatrix(
        unitary.rows,
        unitary.cols,
        lambda row, col: sp.simplify(
            sp.conjugate(unitary[row, col]) * unitary[row, col]
        ),
    )
