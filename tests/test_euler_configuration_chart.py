"""Physical coupled-observation chart versus direct action variation."""

import pytest
import sympy as s

from substrate_framework.euler_phase import (
    physical_configuration_chart,
    physical_scalar_chart,
)


def canonical_form(n):
    return s.zeros(n).row_join(s.eye(n)).col_join(
        (-s.eye(n)).row_join(s.zeros(n)))


def test_moving_two_coordinate_action_against_direct_variation():
    t = s.Symbol("t", real=True)
    omega = canonical_form(2)
    potential = s.Matrix([[2, 1], [1, 3]])
    h = s.diag(potential, s.eye(2))
    b = -omega.inv()*h
    transform = s.Matrix([[1, t], [0, 1]])
    c = transform.row_join(s.zeros(2))
    result = physical_configuration_chart(
        omega, b, c, configuration_rate=c.diff(t),
        configuration_acceleration=c.diff(t, 2), generator_rate=s.zeros(4),
        momentum=s.zeros(2).row_join(s.eye(2)),
    )
    q = s.Matrix([s.Function("q1")(t), s.Function("q2")(t)])
    original = transform.inv()*q
    direct = (original.diff(t).dot(original.diff(t))-(original.T*potential*original)[0])/2
    encoded = ((q.diff(t).T*result.mass*q.diff(t))[0]
               -(q.T*result.gyroscopic_form*q.diff(t))[0]
               -(q.T*result.stiffness*q)[0])/2

    def variation(action):
        return s.Matrix([s.diff(s.diff(action, row.diff(t)), t)-s.diff(action, row) for row in q])

    assert result.ordinary_action_condition is s.true
    assert result.configuration_bracket == s.zeros(2)
    assert s.simplify(variation(direct)-variation(encoded)) == s.zeros(2, 1)
    equation = (result.mass*q.diff(t, 2)
                +(result.mass.diff(t)+result.gyroscopic_form)*q.diff(t)
                +(result.stiffness+result.gyroscopic_form.diff(t)/2)*q)
    assert s.simplify(variation(encoded)-equation) == s.zeros(2, 1)
    assert result.mass == s.Matrix([[1, -t], [-t, t*t+1]])
    assert result.gyroscopic_form != s.zeros(2)
    assert result.momentum_difference != s.zeros(2, 4)
    assert s.simplify(result.generator[:2, :]-s.zeros(2).row_join(s.eye(2))) == s.zeros(2, 4)
    # The missing gyroscopic term is exposed by the actual variation.
    wrong = encoded+(q.T*result.gyroscopic_form*q.diff(t))[0]/2
    assert s.simplify(variation(direct)-variation(wrong)) != s.zeros(2, 1)


def test_scalar_reduction_matches_existing_complete_physical_chart():
    t = s.Symbol("t", real=True)
    omega = s.Matrix([[0, 2], [-2, 0]])
    b = s.Matrix([[0, 1], [-3, 0]])
    c = s.Matrix([[1+t, 0]])
    measured = s.Matrix([[0, 5]])
    scalar = physical_scalar_chart(
        omega, b, c, angle_rate=c.diff(t), angle_acceleration=c.diff(t, 2),
        generator_rate=s.zeros(2), spin=measured,
    )
    result = physical_configuration_chart(
        omega, b, c, configuration_rate=c.diff(t),
        configuration_acceleration=c.diff(t, 2), generator_rate=s.zeros(2), momentum=measured,
    )
    assert result.ordinary_action_condition is s.true
    assert result.mass[0, 0] == scalar.mass
    assert result.stiffness[0, 0] == scalar.stiffness
    assert result.generator == scalar.generator
    assert s.simplify(result.hamiltonian-scalar.hamiltonian) == s.zeros(2)
    assert result.measured_momentum == s.Matrix([[scalar.spin_connection, scalar.spin_inertia]])


def test_noncommuting_physical_positions_retain_full_phase_obstruction():
    omega = canonical_form(2)
    b = -omega.inv()*s.diag(2, 3, 1, 1)
    c = s.Matrix([[1, 0, 0, 1], [0, 1, 0, 0]])
    result = physical_configuration_chart(
        omega, b, c, configuration_rate=s.zeros(2, 4),
        configuration_acceleration=s.zeros(2, 4), generator_rate=s.zeros(4),
        momentum=s.zeros(2).row_join(s.eye(2)),
    )
    assert result.configuration_bracket != s.zeros(2)
    assert result.rate_rate_form != s.zeros(2)
    assert result.ordinary_action_condition is s.false
    assert s.simplify(result.symplectic*result.generator+result.hamiltonian) == s.zeros(4)


def test_initial_position_commutation_does_not_supply_an_all_time_action():
    omega = canonical_form(2)
    b = -omega.inv()*s.diag(2, 3, 1, 1)
    c = s.eye(2).row_join(s.zeros(2))
    cd = s.Matrix([[0, 0, 0, 1], [0, 0, 0, 0]])
    result = physical_configuration_chart(
        omega, b, c, configuration_rate=cd,
        configuration_acceleration=s.zeros(2, 4), generator_rate=s.zeros(4),
        momentum=s.zeros(2).row_join(s.eye(2)),
    )
    assert result.configuration_bracket == s.zeros(2)
    assert result.rate_rate_form == s.zeros(2)
    assert result.symplectic_rate[2:, 2:] != s.zeros(2)
    assert result.ordinary_action_condition is s.false


def test_invalid_phases_rows_and_degenerate_rate_map():
    omega = canonical_form(2)
    b = omega
    c = s.eye(2).row_join(s.zeros(2))
    kwargs = dict(configuration_rate=s.zeros(2, 4), configuration_acceleration=s.zeros(2, 4),
                  generator_rate=s.zeros(4), momentum=s.zeros(2, 4))
    with pytest.raises(ValueError, match="even-dimensional|skew"):
        physical_configuration_chart(s.eye(4), b, c, **kwargs)
    with pytest.raises(ValueError, match="nondegenerate"):
        physical_configuration_chart(s.zeros(4), b, c, **kwargs)
    with pytest.raises(ValueError, match="Hamiltonian"):
        physical_configuration_chart(omega, s.eye(4), c, **kwargs)
    with pytest.raises(ValueError, match="n by 2n"):
        physical_configuration_chart(omega, b, s.zeros(1, 4), **kwargs)
    with pytest.raises(ValueError, match="configuration/rate"):
        physical_configuration_chart(omega, s.zeros(4), c, **kwargs)
