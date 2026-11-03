"""Exact finite local non-Abelian gauge algebra.

The convention is ``D_mu = partial_mu - i*g*W_mu``.  Consequently the finite
connection law has the inhomogeneous sign
``W_mu' = U*W_mu*U.H - (i/g)*(partial_mu U)*U.H``.  This module checks that
law, curvature covariance, and the covariant-derivative commutator for exact
matrices.  It supplies no kinetic-action coefficient, source equation, matter
current, anomaly statement, mass mechanism, or physical gauge-sector map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp

from .su2_doublets import su2_chiral_factor_ledger


def _simplified_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.Matrix(value)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError(f"{name} must be nonempty")
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError(f"{name} must contain exact entries")
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


def _zero(matrix: sp.MatrixBase) -> bool:
    simplified = sp.Matrix(matrix).applyfunc(sp.simplify)
    return simplified == sp.zeros(*simplified.shape)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return sp.simplify(expression)


def _coordinate_pair(
    coordinates: Sequence[sp.Symbol],
) -> tuple[sp.Symbol, sp.Symbol]:
    if len(coordinates) != 2:
        raise ValueError("coordinates must contain exactly two symbols")
    first, second = coordinates
    if not isinstance(first, sp.Symbol) or not isinstance(second, sp.Symbol):
        raise ValueError("coordinates must be SymPy symbols")
    if first == second:
        raise ValueError("coordinates must be distinct")
    return first, second


def nonabelian_covariant_derivative(
    field: Any,
    connection: Any,
    coordinate: sp.Symbol,
    coupling: Any,
) -> sp.ImmutableMatrix:
    """Return ``(partial_mu-i*g*W_mu) field`` in the declared convention."""

    vector = _simplified_matrix(field, "field")
    if vector.cols != 1:
        raise ValueError("field must be a column matrix")
    gauge = _simplified_matrix(connection, "connection")
    if gauge.shape != (vector.rows, vector.rows):
        raise ValueError("connection must be square and match the field")
    if not isinstance(coordinate, sp.Symbol):
        raise ValueError("coordinate must be a SymPy symbol")
    strength = _positive_exact(coupling, "coupling")
    return _simplified_matrix(
        vector.diff(coordinate) - sp.I * strength * gauge * vector,
        "covariant_derivative",
    )


def _field_strength_from_matrices(
    connection_mu: sp.ImmutableMatrix,
    connection_nu: sp.ImmutableMatrix,
    coordinate_mu: sp.Symbol,
    coordinate_nu: sp.Symbol,
    coupling: sp.Expr,
) -> sp.ImmutableMatrix:
    commutator = connection_mu * connection_nu - connection_nu * connection_mu
    return _simplified_matrix(
        connection_nu.diff(coordinate_mu)
        - connection_mu.diff(coordinate_nu)
        - sp.I * coupling * commutator,
        "field_strength",
    )


def nonabelian_field_strength(
    connections: Sequence[Any],
    coordinates: Sequence[sp.Symbol],
    coupling: Any,
) -> sp.ImmutableMatrix:
    """Return ``F_mu_nu=d_mu W_nu-d_nu W_mu-i*g[W_mu,W_nu]``."""

    if len(connections) != 2:
        raise ValueError("connections must contain exactly two matrices")
    coordinate_mu, coordinate_nu = _coordinate_pair(coordinates)
    connection_mu = _simplified_matrix(connections[0], "connection_mu")
    connection_nu = _simplified_matrix(connections[1], "connection_nu")
    if connection_mu.rows != connection_mu.cols:
        raise ValueError("connection_mu must be square")
    if connection_nu.shape != connection_mu.shape:
        raise ValueError("connections must have the same square shape")
    if not _zero(connection_mu - connection_mu.H):
        raise ValueError("connection_mu must be Hermitian")
    if not _zero(connection_nu - connection_nu.H):
        raise ValueError("connection_nu must be Hermitian")
    strength = _positive_exact(coupling, "coupling")
    return _field_strength_from_matrices(
        connection_mu,
        connection_nu,
        coordinate_mu,
        coordinate_nu,
        strength,
    )


@dataclass(frozen=True)
class NonAbelianGaugeLedger:
    """Exact finite transformation, covariance, curvature, and commutator data."""

    coupling: sp.Expr
    coordinates: tuple[sp.Symbol, sp.Symbol]
    field: sp.ImmutableMatrix
    transformed_field: sp.ImmutableMatrix
    connections: tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]
    transformed_connections: tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]
    covariant_derivatives: tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]
    transformed_covariant_derivatives: tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]
    covariance_residuals: tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]
    curvature: sp.ImmutableMatrix
    transformed_curvature: sp.ImmutableMatrix
    curvature_covariance_residual: sp.ImmutableMatrix
    covariant_commutator: sp.ImmutableMatrix
    commutator_curvature_residual: sp.ImmutableMatrix
    curvature_trace_square: sp.Expr
    transformed_curvature_trace_square: sp.Expr
    trace_invariance_residual: sp.Expr


def local_nonabelian_gauge_ledger(
    field: Any,
    connections: Sequence[Any],
    unitary: Any,
    coordinates: Sequence[sp.Symbol],
    coupling: Any,
) -> NonAbelianGaugeLedger:
    """Return exact finite local-gauge identities for two coordinates.

    ``unitary`` is an explicitly unitary matrix in the same carrier as the
    field.  This function verifies mathematical covariance only.  It does not
    assert that the inputs are physical fields or derive their dynamics.
    """

    coordinate_pair = _coordinate_pair(coordinates)
    strength = _positive_exact(coupling, "coupling")
    vector = _simplified_matrix(field, "field")
    if vector.cols != 1:
        raise ValueError("field must be a column matrix")
    transformation = _simplified_matrix(unitary, "unitary")
    if transformation.shape != (vector.rows, vector.rows):
        raise ValueError("unitary must be square and match the field")
    identity = sp.eye(vector.rows)
    if not _zero(transformation.H * transformation - identity):
        raise ValueError("unitary must be exactly unitary")
    if len(connections) != 2:
        raise ValueError("connections must contain exactly two matrices")
    original_connections = tuple(
        _simplified_matrix(value, f"connection_{index}")
        for index, value in enumerate(connections)
    )
    for connection in original_connections:
        if connection.shape != transformation.shape:
            raise ValueError("connections must match the transformation")
        if not _zero(connection - connection.H):
            raise ValueError("connections must be Hermitian")

    transformed_field = _simplified_matrix(
        transformation * vector,
        "transformed_field",
    )
    transformed_connections = tuple(
        _simplified_matrix(
            transformation * connection * transformation.H
            - sp.I
            / strength
            * transformation.diff(coordinate)
            * transformation.H,
            f"transformed_connection_{index}",
        )
        for index, (connection, coordinate) in enumerate(
            zip(original_connections, coordinate_pair, strict=True)
        )
    )
    original_derivatives = tuple(
        nonabelian_covariant_derivative(
            vector,
            connection,
            coordinate,
            strength,
        )
        for connection, coordinate in zip(
            original_connections,
            coordinate_pair,
            strict=True,
        )
    )
    transformed_derivatives = tuple(
        nonabelian_covariant_derivative(
            transformed_field,
            connection,
            coordinate,
            strength,
        )
        for connection, coordinate in zip(
            transformed_connections,
            coordinate_pair,
            strict=True,
        )
    )
    covariance_residuals = tuple(
        _simplified_matrix(
            transformed_derivative - transformation * original_derivative,
            f"covariance_residual_{index}",
        )
        for index, (transformed_derivative, original_derivative) in enumerate(
            zip(transformed_derivatives, original_derivatives, strict=True)
        )
    )

    curvature = nonabelian_field_strength(
        original_connections,
        coordinate_pair,
        strength,
    )
    transformed_curvature = _field_strength_from_matrices(
        transformed_connections[0],
        transformed_connections[1],
        coordinate_pair[0],
        coordinate_pair[1],
        strength,
    )
    curvature_covariance_residual = _simplified_matrix(
        transformed_curvature
        - transformation * curvature * transformation.H,
        "curvature_covariance_residual",
    )

    derivative_mu_nu = nonabelian_covariant_derivative(
        nonabelian_covariant_derivative(
            vector,
            original_connections[1],
            coordinate_pair[1],
            strength,
        ),
        original_connections[0],
        coordinate_pair[0],
        strength,
    )
    derivative_nu_mu = nonabelian_covariant_derivative(
        nonabelian_covariant_derivative(
            vector,
            original_connections[0],
            coordinate_pair[0],
            strength,
        ),
        original_connections[1],
        coordinate_pair[1],
        strength,
    )
    commutator = _simplified_matrix(
        derivative_mu_nu - derivative_nu_mu,
        "covariant_commutator",
    )
    commutator_residual = _simplified_matrix(
        commutator + sp.I * strength * curvature * vector,
        "commutator_curvature_residual",
    )
    trace_square = sp.simplify(sp.trace(curvature * curvature))
    transformed_trace_square = sp.simplify(
        sp.trace(transformed_curvature * transformed_curvature)
    )

    return NonAbelianGaugeLedger(
        coupling=strength,
        coordinates=coordinate_pair,
        field=vector,
        transformed_field=transformed_field,
        connections=original_connections,
        transformed_connections=transformed_connections,
        covariant_derivatives=original_derivatives,
        transformed_covariant_derivatives=transformed_derivatives,
        covariance_residuals=covariance_residuals,
        curvature=curvature,
        transformed_curvature=transformed_curvature,
        curvature_covariance_residual=curvature_covariance_residual,
        covariant_commutator=commutator,
        commutator_curvature_residual=commutator_residual,
        curvature_trace_square=trace_square,
        transformed_curvature_trace_square=transformed_trace_square,
        trace_invariance_residual=sp.simplify(
            transformed_trace_square - trace_square
        ),
    )


def su2_projected_connection(
    components: Sequence[Any],
    projector: Any,
) -> sp.ImmutableMatrix:
    """Return ``sum_a W^a (T_a tensor P)`` on C-REP-002's carrier."""

    if len(components) != 3:
        raise ValueError("components must contain exactly three values")
    carrier = su2_chiral_factor_ledger(projector)
    exact_components = []
    for index, value in enumerate(components):
        expression = sp.sympify(value)
        if expression.has(sp.Float):
            raise ValueError(f"component_{index} must be exact")
        if expression.is_real is not True:
            raise ValueError(f"component_{index} must be explicitly real")
        exact_components.append(expression)
    return _simplified_matrix(
        sum(
            (
                component * generator
                for component, generator in zip(
                    exact_components,
                    carrier.left_generators,
                    strict=True,
                )
            ),
            sp.zeros(*carrier.left_generators[0].shape),
        ),
        "projected_connection",
    )


def su2_projected_unitary(
    isospin_unitary: Any,
    projector: Any,
) -> sp.ImmutableMatrix:
    """Embed a declared SU(2) unitary on ``image(P)`` and identity on its kernel."""

    isospin = _simplified_matrix(isospin_unitary, "isospin_unitary")
    if isospin.shape != (2, 2):
        raise ValueError("isospin_unitary must be two by two")
    if not _zero(isospin.H * isospin - sp.eye(2)):
        raise ValueError("isospin_unitary must be exactly unitary")
    if sp.simplify(isospin.det() - 1) != 0:
        raise ValueError("isospin_unitary must have determinant one")
    carrier = su2_chiral_factor_ledger(projector)
    embedded = sp.kronecker_product(isospin, carrier.projector) + sp.kronecker_product(
        sp.eye(2),
        carrier.complementary_projector,
    )
    return _simplified_matrix(embedded, "projected_unitary")
