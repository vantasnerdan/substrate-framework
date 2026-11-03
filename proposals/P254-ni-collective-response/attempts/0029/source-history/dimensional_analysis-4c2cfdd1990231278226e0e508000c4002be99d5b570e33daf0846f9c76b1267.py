"""Exact linear-algebra helpers for Buckingham-style monomial groups."""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def dimensionless_monomial_basis(dimension_matrix: Any) -> tuple[sp.Matrix, ...]:
    """Return a basis for exponent vectors with zero total dimensions.

    Rows are base dimensions and columns are primitive quantities.
    """

    matrix = sp.Matrix(dimension_matrix)
    return tuple(matrix.nullspace())


def dimensionless_group_count(dimension_matrix: Any) -> int:
    """Return ``number_of_primitives - rank`` for a dimension matrix."""

    matrix = sp.Matrix(dimension_matrix)
    return matrix.cols - matrix.rank()


def monomial_exponents(
    dimension_matrix: Any, target_dimension: Any
) -> sp.Matrix:
    """Return the unique primitive exponents for a target dimension.

    Rows of ``dimension_matrix`` are base dimensions and columns are primitive
    quantities. A unique representation requires full column rank and the
    target to lie in the column span. Dimensionless coefficients are outside
    this exponent calculation and remain unconstrained.
    """

    matrix = sp.Matrix(dimension_matrix)
    target = sp.Matrix(target_dimension)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError("dimension_matrix must be non-empty")
    if target.cols != 1 or target.rows != matrix.rows:
        raise ValueError("target_dimension must be a column matching matrix rows")
    if matrix.rank() < matrix.cols:
        raise ValueError("target monomial exponents are not unique")
    if matrix.row_join(target).rank() > matrix.rank():
        raise ValueError("target dimension is outside the primitive span")
    solution, parameters = matrix.gauss_jordan_solve(target)
    if parameters.rows:
        raise ValueError("target monomial exponents are not unique")
    return sp.Matrix(solution)


def dimensionless_mass_coordinate(
    mass: Any,
    speed: Any,
    action: Any,
    length: Any,
) -> sp.Expr:
    """Return ``mass*speed*length/action`` for a declared positive basis.

    This is a lossless reparameterization at fixed basis values. It does not
    predict the coordinate or establish that the basis is physically complete.
    """

    mass_value = _positive(mass, "mass")
    speed_value = _positive(speed, "speed")
    action_value = _positive(action, "action")
    length_value = _positive(length, "length")
    return sp.simplify(mass_value * speed_value * length_value / action_value)


def mass_from_dimensionless_coordinate(
    coordinate: Any,
    speed: Any,
    action: Any,
    length: Any,
) -> sp.Expr:
    """Return ``coordinate*action/(speed*length)`` for a positive basis."""

    coordinate_value = _positive(coordinate, "coordinate")
    speed_value = _positive(speed, "speed")
    action_value = _positive(action, "action")
    length_value = _positive(length, "length")
    return sp.simplify(
        coordinate_value * action_value / (speed_value * length_value)
    )


def mass_coordinate_from_unit_product(
    coupling: Any,
    unit_product_coefficient: Any,
    mass_energy_coefficient: Any,
) -> sp.Expr:
    """Return a conditional mass coordinate from two declared unit equations.

    The premises are ``U*L=r*S*c/e**2`` and ``U=k*m*c**2``. Their composition
    gives ``m*c*L/S=r/(k*e**2)``. The coefficients and coupling are required
    inputs; this helper does not derive or select them.
    """

    coupling_value = _positive(coupling, "coupling")
    product_value = _positive(
        unit_product_coefficient, "unit_product_coefficient"
    )
    energy_value = _positive(
        mass_energy_coefficient, "mass_energy_coefficient"
    )
    return sp.simplify(product_value / (energy_value * coupling_value**2))
