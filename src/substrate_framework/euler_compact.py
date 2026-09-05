"""Conditional exact compact-isovortical operator algebra (P251/0085, 0103).

For a scalar test function f, write a vector differential operator in RIGHT
normal form, S_j f = sum_alpha d^alpha(s_{j,alpha} f).  This module forms
the coefficient matrix of D S, where D xi=(div xi, div(xi cross omega)).
It also forms the distinct adjoint rows for integral r cross xi and
integral r cross (xi cross omega). Density, when wanted, multiplies those
rows externally.

The API performs no rank test and assumes no Beltrami relation. A null
vector at one point is not itself a compact field: constructing S requires
coefficient SECTIONS satisfying the matrix identity on an open set. If
such sections exist, applying S to a smooth compact scalar gives the two
divergence constraints exactly. Background existence, rank persistence,
positive Hessians, moment controllability and continuum closure remain
separate conditional theorems, not import-time conclusions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import product
from math import comb
from operator import index

import sympy as sp

MultiIndex = tuple[int, int, int]


def derivative_indices(order: int) -> tuple[MultiIndex, ...]:
    """Return three-dimensional multi-indices through order, graded lexically."""
    if isinstance(order, bool):
        raise ValueError("order must be a nonnegative integer")
    try:
        order = index(order)
    except TypeError as exc:
        raise ValueError("order must be a nonnegative integer") from exc
    if order < 0:
        raise ValueError("order must be a nonnegative integer")
    return tuple(sorted((a for a in product(range(order + 1), repeat=3)
                         if sum(a) <= order), key=lambda a: (sum(a), a)))


def _exact_vector(values, name: str) -> sp.ImmutableMatrix:
    vector = sp.ImmutableMatrix(values)
    if vector.shape != (3, 1):
        raise ValueError(f"{name} must have three components")
    if any(entry.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo) for entry in vector):
        raise ValueError(f"{name} must contain finite exact symbolic expressions")
    return vector


@dataclass(frozen=True)
class CompactJetSystem:
    """An exact right-normal matrix and its physical angular adjoint rows.

    Columns are (component, alpha), component-major, with alpha taken from
    ``operator_indices``. The first output block is div S and the second
    div(S cross omega); each is indexed by ``output_indices`` (orders one
    through N+1). An outer divergence leaves no zeroth-order output.

    The three angular rows are ordered by Cartesian axis. Both angular
    matrices are adjoint coefficients evaluated at ``radius`` = point
    minus the chosen moment origin. They act on the SAME coefficient
    column as ``constraints``; neither includes a density factor.
    """

    operator_indices: tuple[MultiIndex, ...]
    output_indices: tuple[MultiIndex, ...]
    constraints: sp.ImmutableMatrix
    generator_angular: sp.ImmutableMatrix
    velocity_angular: sp.ImmutableMatrix


def compact_isovortical_jet_system(
    vorticity_jets: Mapping[MultiIndex, Sequence],
    order: int = 6,
    *,
    radius: Sequence = (0, 0, 0),
) -> CompactJetSystem:
    """Construct D S and both angular rows from actual vorticity derivatives.

    ``vorticity_jets[alpha]`` is d^alpha omega, NOT a Taylor coefficient
    divided by alpha factorial. Every jet through ``order`` is required;
    missing derivatives are not silently treated as zero. Values may be
    exact point data or symbolic expressions of position. Floating inputs
    are rejected: use explicit SymPy rationals for an exact calculation.

    The commutation identity used for the second divergence block is
    m_c d^alpha = sum_{gamma<=alpha} (-1)^|gamma| binom(alpha,gamma)
                 d^(alpha-gamma) m_(d^gamma c).
    It differentiates omega, not the unknown right-normal coefficients.
    """
    indices = derivative_indices(order)
    order = index(order)
    outputs = derivative_indices(order + 1)[1:]
    missing = [alpha for alpha in indices if alpha not in vorticity_jets]
    if missing:
        raise ValueError(f"missing vorticity derivative {missing[0]}")
    jets = {alpha: _exact_vector(vorticity_jets[alpha], f"jet {alpha}")
            for alpha in indices}
    r = _exact_vector(radius, "radius")
    out_index = {alpha: row for row, alpha in enumerate(outputs)}
    width, height = 3 * len(indices), len(outputs)
    constraints = sp.MutableSparseMatrix(2 * height, width, {})
    generator = sp.MutableSparseMatrix(3, width, {})
    velocity = sp.MutableSparseMatrix(3, width, {})

    for component in range(3):
        basis = sp.eye(3)[:, component]
        for column, alpha in enumerate(indices):
            col = component * len(indices) + column
            beta = tuple(alpha[k] + int(k == component) for k in range(3))
            constraints[out_index[beta], col] += 1
            for gamma in product(*(range(a + 1) for a in alpha)):
                factor = (-1) ** sum(gamma)
                for a, g in zip(alpha, gamma, strict=True):
                    factor *= comb(a, g)
                force = basis.cross(jets[gamma])
                for axis in range(3):
                    beta = tuple(alpha[k] - gamma[k] + int(k == axis)
                                 for k in range(3))
                    constraints[height + out_index[beta], col] += factor * force[axis]

            degree, sign = sum(alpha), (-1) ** sum(alpha)
            for axis in range(3):
                unit = sp.eye(3)[:, axis]
                if degree == 0:
                    generator[axis, col] = unit.cross(r)[component]
                elif degree == 1:
                    generator[axis, col] = -unit.cross(sp.Matrix(alpha))[component]

                # d^alpha [omega cross (e_axis cross r)] at the actual radius.
                value = -r[component] * jets[alpha][axis]
                if component == axis:
                    value += jets[alpha].dot(r)
                    for k in range(3):
                        if alpha[k]:
                            lower = tuple(alpha[i] - int(i == k) for i in range(3))
                            value += alpha[k] * jets[lower][k]
                if alpha[component]:
                    lower = tuple(alpha[i] - int(i == component) for i in range(3))
                    value -= alpha[component] * jets[lower][axis]
                velocity[axis, col] = sign * value

    return CompactJetSystem(indices, outputs, sp.ImmutableMatrix(constraints),
                            sp.ImmutableMatrix(generator), sp.ImmutableMatrix(velocity))
