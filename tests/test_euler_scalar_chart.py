import pytest
import sympy as s

from substrate_framework.euler_phase import physical_scalar_chart


def test_moving_physical_chart_matches_differentiation_and_variation():
    t = s.Symbol("t", real=True)
    omega = s.Matrix([[0, 3], [-3, 0]])
    b = s.Matrix([[0, 1+t**2], [-2, 0]])
    c = s.Matrix([[1, t]])
    measured = s.Matrix([[2*t, 5]])
    result = physical_scalar_chart(
        omega, b, c, angle_rate=c.diff(t), angle_acceleration=c.diff(t, 2),
        generator_rate=b.diff(t), spin=measured)
    transform = c.col_join(c.diff(t)+c*b)
    expected = (transform.diff(t)+transform*b)*transform.inv()
    assert s.simplify(result.generator-expected) == s.zeros(2)
    assert s.simplify(result.mass.diff(t)-result.mass_rate) == 0
    assert s.simplify(result.generator[1, 1]+result.mass_rate/result.mass) == 0
    assert s.simplify(result.hamiltonian-s.Matrix([
        [result.stiffness, result.mass_rate/2],
        [result.mass_rate/2, result.mass]])) == s.zeros(2)
    theta = s.Function("theta")(t)
    action = result.mass*theta.diff(t)**2/2-result.stiffness*theta**2/2
    variation = s.diff(s.diff(action, theta.diff(t)), t)-s.diff(action, theta)
    dynamic = result.mass*(theta.diff(t, 2)-expected[1, 0]*theta
                          -expected[1, 1]*theta.diff(t))
    assert s.simplify(variation-dynamic) == 0
    assert s.simplify(result.mass_rate) != 0  # omitted moving term is exposed
    phase = s.Matrix(s.symbols("q p"))
    observation = result.coordinates*phase
    actual = (measured*phase)[0]
    reconstructed = result.spin_connection*observation[0]+result.spin_inertia*observation[1]
    assert s.simplify(actual-reconstructed) == 0
    assert s.simplify(result.angle_spin_bracket-result.spin_inertia/result.mass) == 0


def test_positive_physical_spin_is_not_automatically_canonical():
    omega = s.Matrix([[0, 2], [-2, 0]])
    b = s.Matrix([[0, 3], [-3, 0]])
    arguments = dict(angle_rate=s.zeros(1, 2), angle_acceleration=s.zeros(1, 2),
                     generator_rate=s.zeros(2))
    wrong = physical_scalar_chart(omega, b, [[1, 0]], spin=[[0, 1]], **arguments)
    matched = physical_scalar_chart(omega, b, [[1, 0]], spin=[[7, 2]], **arguments)
    assert wrong.mass == s.Rational(2, 3)
    assert wrong.spin_inertia == s.Rational(1, 3)
    assert wrong.angle_spin_bracket == s.Rational(1, 2)
    assert matched.mass == matched.spin_inertia
    assert matched.angle_spin_bracket == 1
    assert matched.spin_connection == 7  # normalized bracket retains current


def test_actual_observation_winding_changes_mass_and_can_destroy_chart():
    t = s.Symbol("t", real=True)
    omega = s.Matrix([[0, 1], [-1, 0]])
    b = omega
    def chart(frequency):
        row = s.Matrix([[s.cos(frequency*t), s.sin(frequency*t)]])
        return physical_scalar_chart(
            omega, b, row, angle_rate=row.diff(t), angle_acceleration=row.diff(t, 2),
            generator_rate=s.zeros(2), spin=[[0, 1]])
    assert chart(0).mass == 1
    assert chart(-2).mass == -1
    with pytest.raises(ValueError, match="Wronskian"):
        chart(-1)


def test_physical_chart_rejects_unlicensed_matrix_domains():
    omega = s.Matrix([[0, 1], [-1, 0]])
    arguments = dict(angle_rate=s.zeros(1, 2), angle_acceleration=s.zeros(1, 2),
                     generator_rate=s.zeros(2), spin=[[0, 1]])
    for invalid in (s.eye(2), s.zeros(2), s.zeros(3)):
        with pytest.raises(ValueError, match="symplectic"):
            physical_scalar_chart(invalid, omega, [[1, 0]], **arguments)
    with pytest.raises(ValueError, match="Hamiltonian"):
        physical_scalar_chart(omega, s.eye(2), [[1, 0]], **arguments)
    with pytest.raises(ValueError, match="rows"):
        physical_scalar_chart(omega, omega, [[1], [0]], **arguments)
    with pytest.raises(ValueError, match="Wronskian"):
        physical_scalar_chart(omega, s.zeros(2), [[1, 0]], **arguments)
