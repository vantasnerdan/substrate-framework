import pytest
import sympy as sp

from substrate_framework.euler_retained_memory import (
    euler_pressure_poisson_source,
    linear_retained_memory,
)


def test_pressure_source_uses_euler_cross_contraction_not_frobenius_norm():
    rho = sp.symbols("rho", positive=True)
    gradient = sp.Matrix([[1, 2, 0], [3, -1, 4], [0, 5, 0]])
    source = euler_pressure_poisson_source(gradient, density=rho)
    assert sp.simplify(source - rho * sp.trace(gradient * gradient)) == 0
    assert sp.simplify(source - rho * sp.trace(gradient.T * gradient)) != 0


def test_two_block_memory_has_exact_markov_noise_and_kernel_signs():
    a, b, c, d, t = sp.symbols("a b c d t", real=True)
    operator = sp.Matrix([[a, b], [c, d]])
    projection = sp.diag(1, 0)
    memory = linear_retained_memory(operator, projection, t)
    assert memory.markov == sp.ImmutableMatrix([[a, 0], [0, 0]])
    assert memory.unresolved_propagator == sp.ImmutableMatrix(
        [[0, b * sp.exp(d * t)], [0, 0]])
    assert memory.memory_kernel == sp.ImmutableMatrix(
        [[b * c * sp.exp(d * t), 0], [0, 0]])


def test_invalid_projection_and_density_are_rejected():
    with pytest.raises(ValueError, match="idempotent"):
        linear_retained_memory(sp.eye(2), sp.diag(2, 0), sp.symbols("t", real=True))
    with pytest.raises(ValueError, match="positive"):
        euler_pressure_poisson_source(sp.eye(3), density=0)
