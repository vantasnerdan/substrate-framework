import pytest
import sympy as sp

from substrate_framework.euler_quantum_two_state import (
    affine_kks_poisson_bracket,
    euler_circulation_scale,
    finite_ccr_trace_obstruction,
    kks_plane_dynamics,
)


def test_positive_kks_plane_is_one_classical_oscillator():
    b, hq, hp = sp.symbols("b hq hp", positive=True)
    data = kks_plane_dynamics(b, ((hq, 0), (0, hp)))
    assert data.generator == sp.Matrix([[0, hp / b], [-hq / b, 0]])
    assert sp.simplify(data.generator**2) == sp.Matrix(
        [[-hq * hp / b**2, 0], [0, -hq * hp / b**2]]
    )
    assert data.frequency_squared == hq * hp / b**2


def test_angle_momentum_rows_form_a_central_heisenberg_bracket():
    b, gain = sp.symbols("b gain", nonzero=True)
    assert affine_kks_poisson_bracket((1, 0), (0, b), b) == 1
    assert affine_kks_poisson_bracket((gain, 0), (0, b), b) == gain
    assert affine_kks_poisson_bracket((1, 0), (gain, 0), b) == 0


def test_exact_ccr_has_no_finite_dimensional_matrix_representation():
    hbar = sp.symbols("hbar", positive=True)
    assert finite_ccr_trace_obstruction(2, hbar) == -2 * sp.I * hbar
    with pytest.raises(ValueError, match="positive integer"):
        finite_ccr_trace_obstruction(0, hbar)


def test_kelvin_circulation_is_continuously_rescaled_by_euler_similarity():
    gamma, a, b = sp.symbols("gamma a b", positive=True)
    assert euler_circulation_scale(gamma, a, b) == a * gamma / b


def test_domains_expose_singular_or_nonsymmetric_inputs():
    with pytest.raises(ValueError, match="nonzero"):
        kks_plane_dynamics(0, ((1, 0), (0, 1)))
    with pytest.raises(ValueError, match="symmetric"):
        kks_plane_dynamics(1, ((1, 1), (0, 1)))
