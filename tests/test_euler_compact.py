"""Independent differential composition and angular-adjoint checks."""

import pytest
import sympy as sp

from substrate_framework.euler_compact import (
    compact_isovortical_jet_system,
    derivative_indices,
)

X = sp.symbols("x y z", real=True)
x, y, z = X
OMEGA = sp.Matrix([1 + y + z**2, x - y*z, 2 - x**2 + x*y])


def differentiate(expression, alpha):
    for coordinate, count in zip(X, alpha, strict=True):
        expression = sp.diff(expression, coordinate, count)
    return expression


def jets(field, order, point=None):
    result = {alpha: field.applyfunc(lambda e: differentiate(e, alpha))
              for alpha in derivative_indices(order)}
    if point is not None:
        result = {alpha: value.subs(dict(zip(X, point, strict=True)))
                  for alpha, value in result.items()}
    return result


def divergence(field):
    return sum(sp.diff(field[i], X[i]) for i in range(3))


@pytest.mark.parametrize("order", [0, 1, 2])
def test_matrix_composes_actual_right_normal_operators(order):
    system = compact_isovortical_jet_system(jets(OMEGA, order), order)
    indices = system.operator_indices
    f = (1 + x + y + z)**4 + x*y**2*z**3
    # Spatially VARIABLE coefficients catch accidental left-normal ordering.
    coefficients = [(j + 1)*(1 + x*y) + (sum(alpha) + 1)*z + alpha[1]*x
                    for j in range(3) for alpha in indices]
    xi = sp.Matrix([sum(differentiate(coefficients[j*len(indices) + n]*f, alpha)
                        for n, alpha in enumerate(indices)) for j in range(3)])
    direct = [sp.expand(divergence(xi)), sp.expand(divergence(xi.cross(OMEGA)))]
    right_coefficients = system.constraints * sp.Matrix(coefficients)
    count = len(system.output_indices)
    composed = [sp.expand(sum(differentiate(right_coefficients[b*count + n]*f, beta)
                              for n, beta in enumerate(system.output_indices)))
                for b in range(2)]
    assert all(sp.expand(a - b) == 0 for a, b in zip(direct, composed, strict=True))


@pytest.mark.parametrize("point,origin", [((0, 0, 0), (0, 0, 0)),
                                         ((1, -2, 3), (2, 1, -1))])
def test_angular_rows_by_independent_physical_adjoint(point, origin):
    order = 2
    radius = sp.Matrix(point) - sp.Matrix(origin)
    system = compact_isovortical_jet_system(jets(OMEGA, order, point), order,
                                           radius=radius)
    r = sp.Matrix(X) - sp.Matrix(origin)
    evaluation = dict(zip(X, point, strict=True))
    for axis in range(3):
        rotation = sp.eye(3)[:, axis].cross(r)
        force_test = OMEGA.cross(rotation)
        for component in range(3):
            for n, alpha in enumerate(system.operator_indices):
                column = component*len(system.operator_indices) + n
                for field, matrix in ((rotation, system.generator_angular),
                                      (force_test, system.velocity_angular)):
                    direct = (-1)**sum(alpha) * differentiate(field[component], alpha)
                    assert sp.expand(matrix[axis, column] - direct.subs(evaluation)) == 0
    # These are genuinely different physical rows, not aliases.
    assert system.generator_angular != system.velocity_angular


def test_actual_derivative_normalization_and_sign_are_exposed():
    order = 2
    field = sp.Matrix([0, 0, x**2])
    actual = jets(field, order)
    correct = compact_isovortical_jet_system(actual, order)
    factorial_mutation = dict(actual)
    factorial_mutation[(2, 0, 0)] = [0, 0, 1]
    wrong = compact_isovortical_jet_system(factorial_mutation, order)
    assert correct.constraints != wrong.constraints
    negated = compact_isovortical_jet_system({a: -w for a, w in actual.items()}, order)
    n = len(correct.output_indices)
    assert correct.constraints[:n, :] == negated.constraints[:n, :]
    assert correct.constraints[n:, :] == -negated.constraints[n:, :]
    assert correct.generator_angular == negated.generator_angular
    assert correct.velocity_angular == -negated.velocity_angular


def test_zero_order_constant_background_and_shifted_origin():
    omega = sp.Matrix(sp.symbols("w0:3"))
    r = sp.Matrix(sp.symbols("r0:3"))
    result = compact_isovortical_jet_system({(0, 0, 0): omega}, 0, radius=r)
    assert result.constraints.shape == (6, 3)
    assert result.generator_angular.shape == result.velocity_angular.shape == (3, 3)
    for axis in range(3):
        rotation = sp.eye(3)[:, axis].cross(r)
        assert sp.Matrix(result.generator_angular[axis, :]).T == rotation
        assert sp.Matrix(result.velocity_angular[axis, :]).T == omega.cross(rotation)


@pytest.mark.parametrize("order", [-1, 1.5, True])
def test_invalid_order(order):
    with pytest.raises(ValueError, match="nonnegative integer"):
        derivative_indices(order)


def test_missing_malformed_or_inexact_jet_is_not_silently_accepted():
    with pytest.raises(ValueError, match="missing vorticity derivative"):
        compact_isovortical_jet_system({(0, 0, 0): [0, 0, 1]}, 1)
    with pytest.raises(ValueError, match="three components"):
        compact_isovortical_jet_system({(0, 0, 0): [0, 1]}, 0)
    with pytest.raises(ValueError, match="finite exact"):
        compact_isovortical_jet_system({(0, 0, 0): [0, 0, 1.0]}, 0)
    with pytest.raises(ValueError, match="finite exact"):
        compact_isovortical_jet_system({(0, 0, 0): [0, 0, sp.oo]}, 0)
    with pytest.raises(ValueError, match="three components"):
        compact_isovortical_jet_system({(0, 0, 0): [0, 0, 1]}, 0, radius=(0, 0))


def test_returned_matrices_are_immutable_and_rank_is_not_assumed():
    result = compact_isovortical_jet_system({(0, 0, 0): [0, 0, 0]}, 0)
    assert result.constraints.rank() == 3
    assert result.velocity_angular == sp.zeros(3, 3)
    with pytest.raises(TypeError):
        result.constraints[0, 0] = 7
