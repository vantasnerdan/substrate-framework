import pytest
import sympy as s

from substrate_framework.euler_phase import moving_phase_pullback


def test_moving_pullback_against_direct_euler_lagrange_variation():
    t = s.Symbol("t", real=True)
    pair = s.Matrix([[0, 1], [-1, 0]])
    omega = s.diag(pair, 2 * pair)
    h = s.diag(3, 5, 7, 11)
    e = s.Matrix([[1, t], [0, 1], [t, 0], [0, t]])
    result = moving_phase_pullback(omega, h, e, e.diff(t))
    q, p = s.Function("q")(t), s.Function("p")(t)
    z = s.Matrix([q, p])
    x = e * z
    action = -(x.T * omega * x.diff(t))[0] / 2 - (x.T * h * x)[0] / 2
    el = s.Matrix([s.diff(s.diff(action, v.diff(t)), t) - s.diff(action, v)
                   for v in z])
    encoded = (result.symplectic * z.diff(t)
               + (result.hamiltonian + result.symplectic_rate / 2) * z)
    assert s.simplify(el - encoded) == s.zeros(2, 1)
    assert result.symplectic_rate == result.symplectic.diff(t)
    assert result.symplectic_rate != s.zeros(2)
    assert s.simplify(result.coordinates * e) == s.eye(2)
    assert s.simplify(e.T * omega * result.residual) == s.zeros(2)
    assert result.residual != s.zeros(4, 2)
    # Omitting the moving symplectic form changes the defining variation.
    wrong = result.symplectic * z.diff(t) + result.hamiltonian * z
    assert s.simplify(el - wrong) != s.zeros(2, 1)


def test_periodic_winding_changes_energy_sign_not_physical_motion():
    t = s.Symbol("t", real=True)
    omega = s.Matrix([[0, 1], [-1, 0]])
    h = s.eye(2)
    rotation = s.Matrix([[s.cos(2*t), s.sin(2*t)],
                         [-s.sin(2*t), s.cos(2*t)]])
    result = moving_phase_pullback(omega, h, rotation, rotation.diff(t))
    assert result.hamiltonian == -s.eye(2)
    assert result.symplectic == omega
    ambient_generator = -omega.inv() * h
    assert s.simplify(rotation.diff(t) + rotation * result.generator
                      - ambient_generator * rotation) == s.zeros(2)
    assert result.residual == s.zeros(2)
    # E is the physical observation, even though coordinate energy flipped.
    z = s.Matrix(s.symbols("q p"))
    physical_rate = rotation.diff(t)*z + rotation*result.generator*z
    assert s.simplify(physical_rate - ambient_generator*rotation*z) == s.zeros(2, 1)
    assert s.simplify(rotation.subs(t, t+s.pi)-rotation) == s.zeros(2)


def test_full_moving_complement_reconstructs_the_ambient_equation():
    t = s.Symbol("t", real=True)
    omega = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]]))
    h = s.diag(2, 3, 5, 7)
    e = s.Matrix([[1, 0], [0, 1], [t, 0], [0, t]])
    result = moving_phase_pullback(omega, h, e, e.diff(t))
    pi, projection = result.coordinates, result.projection
    a = -omega.inv()*h
    z = s.Matrix(s.symbols("q p"))
    r = (s.eye(4)-projection)*s.Matrix(s.symbols("r:4"))
    zdot = result.generator*z + (pi.diff(t)+pi*a)*r
    rdot = -result.residual*z + ((s.eye(4)-projection)*a-projection.diff(t))*r
    assert s.simplify(e.diff(t)*z+e*zdot+rdot-a*(e*z+r)) == s.zeros(4, 1)
    assert s.simplify(pi*rdot+pi.diff(t)*r) == s.zeros(2, 1)


def test_invalid_forms_and_static_restriction():
    omega = s.Matrix([[0, 2], [-2, 0]])
    result = moving_phase_pullback(omega, s.eye(2), s.eye(2), s.zeros(2))
    assert result.generator == -omega.inv()
    assert result.residual == s.zeros(2)
    for invalid in (s.eye(2), s.zeros(2)):
        with pytest.raises(ValueError):
            moving_phase_pullback(invalid, s.eye(2), s.eye(2), s.zeros(2))
    with pytest.raises(ValueError, match="restricted"):
        moving_phase_pullback(omega, s.eye(2), s.Matrix([[1], [0]]), s.zeros(2, 1))
    with pytest.raises(ValueError, match="symmetric"):
        moving_phase_pullback(omega, s.Matrix([[1, 1], [0, 1]]), s.eye(2), s.zeros(2))
