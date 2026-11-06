"""Dimension-aware exact coordinate geometry.

Authority status: conditional, unpromoted infrastructure linked to issue #40.
Callers supply an exact symmetric invertible metric and a coordinate chart.
These helpers derive its Levi-Civita connection and Ricci contraction; they do
not select a physical metric, signature, field equation, or source.
"""

from __future__ import annotations

from typing import Any, Iterable

import sympy as sp


def coordinate_symbols(coordinates: Iterable[Any]) -> tuple[sp.Symbol, ...]:
    """Return at least two distinct SymPy coordinate symbols."""

    result = tuple(coordinates)
    if len(result) < 2:
        raise ValueError("coordinates must contain at least two symbols")
    if any(not isinstance(coordinate, sp.Symbol) for coordinate in result):
        raise ValueError("each coordinate must be a SymPy Symbol")
    if len(set(result)) != len(result):
        raise ValueError("coordinates must be distinct")
    return result


def exact_metric_matrix(metric: Any, *, dimension: int | None = None) -> sp.Matrix:
    """Return an exact symmetric invertible square metric matrix."""

    result = sp.Matrix(metric)
    rows, columns = result.shape
    if rows != columns or rows < 2:
        raise ValueError("metric must be a square matrix of dimension at least two")
    if dimension is not None and result.shape != (dimension, dimension):
        raise ValueError(f"metric must be {dimension} by {dimension}")
    if result.has(sp.Float):
        raise ValueError("metric entries must be exact rather than floating")
    if sp.simplify(result - result.T) != sp.zeros(rows):
        raise ValueError("metric must be symmetric")
    if sp.simplify(result.det()) == 0:
        raise ValueError("metric must be invertible")
    return result


def metric_inverse(metric: Any) -> sp.Matrix:
    """Return the simplified exact inverse of a supplied metric."""

    return exact_metric_matrix(metric).inv().applyfunc(sp.simplify)


def metric_christoffel_symbols(
    metric: Any,
    coordinates: Iterable[Any],
) -> tuple[tuple[tuple[sp.Expr, ...], ...], ...]:
    r"""Return ``Gamma^a_bc`` for the supplied exact coordinate metric."""

    coords = coordinate_symbols(coordinates)
    metric_matrix = exact_metric_matrix(metric, dimension=len(coords))
    inverse = metric_inverse(metric_matrix)
    dimension = len(coords)
    return tuple(
        tuple(
            tuple(
                sp.simplify(
                    sum(
                        inverse[upper, sigma]
                        * (
                            sp.diff(metric_matrix[nu, sigma], coords[mu])
                            + sp.diff(metric_matrix[mu, sigma], coords[nu])
                            - sp.diff(metric_matrix[mu, nu], coords[sigma])
                        )
                        for sigma in range(dimension)
                    )
                    / 2
                )
                for nu in range(dimension)
            )
            for mu in range(dimension)
        )
        for upper in range(dimension)
    )


def metric_christoffel_from_derivatives(
    inverse_metric: Any,
    metric_derivatives: Any,
) -> sp.ImmutableDenseNDimArray:
    r"""Return ``Gamma^a_bc`` from local ``g^ab`` and ``d_r g_mn`` data.

    The derivative array uses order ``[rho, mu, nu]`` and must be symmetric in
    its last two indices. This local form is useful when no coordinate-function
    representation of the metric is available.
    """

    inverse = exact_metric_matrix(inverse_metric)
    dimension = inverse.rows
    derivatives = sp.ImmutableDenseNDimArray(metric_derivatives)
    if derivatives.shape != (dimension, dimension, dimension):
        raise ValueError("metric derivatives must have shape (n, n, n)")
    flattened = sp.flatten(derivatives.tolist())
    if any(sp.sympify(value).has(sp.Float) for value in flattened):
        raise ValueError("metric derivatives must be exact rather than floating")
    for rho in range(dimension):
        for mu in range(dimension):
            for nu in range(dimension):
                if sp.simplify(
                    derivatives[rho, mu, nu] - derivatives[rho, nu, mu]
                ) != 0:
                    raise ValueError(
                        "metric derivatives must be symmetric in metric indices"
                    )
    entries = [
        sp.simplify(
            sum(
                inverse[upper, sigma]
                * (
                    derivatives[mu, sigma, nu]
                    + derivatives[nu, sigma, mu]
                    - derivatives[sigma, mu, nu]
                )
                for sigma in range(dimension)
            )
            / 2
        )
        for upper in range(dimension)
        for mu in range(dimension)
        for nu in range(dimension)
    ]
    return sp.ImmutableDenseNDimArray(entries, (dimension, dimension, dimension))


def metric_ricci_tensor(metric: Any, coordinates: Iterable[Any]) -> sp.Matrix:
    """Return the exact Ricci tensor in the displayed connection convention."""

    coords = coordinate_symbols(coordinates)
    metric_matrix = exact_metric_matrix(metric, dimension=len(coords))
    gamma = metric_christoffel_symbols(metric_matrix, coords)
    dimension = len(coords)
    tensor = sp.zeros(dimension, dimension)
    for mu in range(dimension):
        for nu in range(dimension):
            entry = sp.Integer(0)
            for upper in range(dimension):
                entry += sp.diff(gamma[upper][mu][nu], coords[upper])
                entry -= sp.diff(gamma[upper][mu][upper], coords[nu])
                for sigma in range(dimension):
                    entry += gamma[upper][upper][sigma] * gamma[sigma][mu][nu]
                    entry -= gamma[upper][nu][sigma] * gamma[sigma][mu][upper]
            tensor[mu, nu] = sp.simplify(entry)
    return tensor


def metric_ricci_scalar(metric: Any, coordinates: Iterable[Any]) -> sp.Expr:
    """Return the exact scalar ``g^ab R_ab`` for a supplied metric."""

    coords = coordinate_symbols(coordinates)
    metric_matrix = exact_metric_matrix(metric, dimension=len(coords))
    inverse = metric_inverse(metric_matrix)
    tensor = metric_ricci_tensor(metric_matrix, coords)
    dimension = len(coords)
    return sp.simplify(
        sum(
            inverse[mu, nu] * tensor[mu, nu]
            for mu in range(dimension)
            for nu in range(dimension)
        )
    )
