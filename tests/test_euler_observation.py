import pytest
import sympy as s

from substrate_framework.euler_observation import (
    material_tag_fourier_dipole,
    material_tag_moments,
)


def test_material_shape_rate_and_spin_from_moving_mass_points():
    t = s.Symbol("t", real=True)
    masses = [1, 2, 3]
    points = [s.Matrix([t, 1+t*t, 2]), s.Matrix([1+t, -t, 3*t]),
              s.Matrix([-t*t, 2, 1-t])]
    velocities = [x.diff(t) for x in points]
    result = material_tag_moments(masses, points, velocities)
    assert s.simplify(result.momentum-result.mass*result.centroid.diff(t)) == s.zeros(3, 1)
    assert s.simplify(result.shape_rate-result.second_mass.diff(t)) == s.zeros(3)
    total_angular = sum((m*x.cross(v) for m, x, v in zip(masses, points, velocities)),
                        s.zeros(3, 1))
    assert s.simplify(total_angular-result.centroid.cross(result.momentum)
                      - result.spin) == s.zeros(3, 1)
    # Translating/boosting the actual tag does not change intrinsic spin.
    shift = s.Matrix([2+t, t*t, 3])
    boosted = material_tag_moments(masses, [x+shift for x in points],
                                   [v+shift.diff(t) for v in velocities])
    assert s.simplify(boosted.spin-result.spin) == s.zeros(3, 1)


def test_dipole_against_defining_fourier_sum_with_shape_mutation():
    epsilon = s.Symbol("epsilon", real=True)
    k = s.Matrix([2, -1, 3])
    masses = [1, 2, 4]
    points = [s.Matrix([1, 2, 0]), s.Matrix([-1, 0, 3]), s.Matrix([0, -2, 1])]
    rates = [s.Matrix([2, 1, -1]), s.Matrix([0, 3, 2]), s.Matrix([-2, 1, 4])]
    result = material_tag_moments(masses, points, rates)
    resolved = sum((m*v*s.exp(-s.I*epsilon*k.dot(x))
                    for m, x, v in zip(masses, points, rates)), s.zeros(3, 1))
    collapsed = result.momentum*s.exp(-s.I*epsilon*k.dot(result.centroid))
    derivative = (resolved-collapsed).diff(epsilon).subs(epsilon, 0)
    dipole = material_tag_fourier_dipole(k, result.spin, result.shape_rate)
    assert s.simplify(derivative-dipole) == s.zeros(3, 1)
    assert s.simplify(dipole+s.I*result.first_momentum*k) == s.zeros(3, 1)
    # A spin-only representation is false for an actually deforming tag.
    assert s.simplify(derivative-s.I*k.cross(result.spin)/2) != s.zeros(3, 1)
    assert s.simplify(derivative+s.I*k.cross(result.spin)/2
                      +s.I*result.shape_rate*k/2) != s.zeros(3, 1)


def test_symmetric_euler_stress_has_no_isotropic_axial_zero_jet():
    # The sole isotropic rank-three tensor is epsilon_ijl; symmetrizing
    # the Euler momentum stress annihilates it, not its fluctuations.
    q = s.Matrix(s.symbols("q:3"))
    stress = s.Matrix(3, 3, lambda i, j: sum(s.LeviCivita(i, j, axis)*q[axis]
                                          for axis in range(3)))
    assert stress+stress.T == s.zeros(3)
    assert stress != s.zeros(3)
    # Exact SO(3) vector-response average: trace(ab^T)*I/3. A triad
    # rotation witness suffices for a diagonal rank-one response.
    rotations = [s.eye(3), s.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
                 s.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])]
    response = s.diag(1, 0, 0)
    assert sum((r*response*r.T for r in rotations), s.zeros(3))/3 == s.eye(3)/3


def test_invalid_material_inputs():
    with pytest.raises(ValueError):
        material_tag_moments([], [], [])
    with pytest.raises(ValueError, match="positive"):
        material_tag_moments([0], [[0, 0, 0]], [[0, 0, 0]])
    with pytest.raises(ValueError, match="equal"):
        material_tag_moments([1], [[0, 0, 0]], [])
    with pytest.raises(ValueError, match="symmetric"):
        material_tag_fourier_dipole([1, 0, 0], [0, 1, 0], [[0, 1, 0], [0, 0, 0], [0, 0, 0]])
