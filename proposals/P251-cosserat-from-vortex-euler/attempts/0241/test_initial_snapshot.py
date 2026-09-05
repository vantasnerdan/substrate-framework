import pytest
import sympy as s

from substrate_framework.euler_joint import (
    material_momentum_second_jet,
    passive_output_weights,
)
from substrate_framework.euler_observation import material_tag_moments
from substrate_framework.euler_passive_control import passive_packet


def test_full_second_momentum_jet_against_defining_moving_fourier_sum():
    e, t = s.symbols('epsilon t', real=True)
    k = s.Matrix([2, -1, 3])
    masses = [1, 2, 3]
    points = [s.Matrix([t, 1+t*t, 2]), s.Matrix([1+t, -t, 3*t]),
              s.Matrix([-t*t, 2, 1-t])]
    rates = [x.diff(t) for x in points]
    jet = material_momentum_second_jet(masses, points, rates, k)
    moments = material_tag_moments(masses, points, rates)
    exact = sum((m*v*s.exp(-s.I*e*k.dot(x))
                 for m, x, v in zip(masses, points, rates)), s.zeros(3, 1))
    exact -= moments.momentum*s.exp(-s.I*e*k.dot(moments.centroid))
    assert s.simplify(exact.subs(e, 0)) == s.zeros(3, 1)
    assert s.simplify(exact.diff(e).subs(e, 0)-jet.first) == s.zeros(3, 1)
    assert s.simplify(exact.diff(e, 2).subs(e, 0)/2-jet.second) == s.zeros(3, 1)
    # The acceleration correction requires time derivatives of the full jet.
    assert s.simplify(exact.diff(t).diff(e, 2).subs(e, 0)/2
                      -jet.second.diff(t)) == s.zeros(3, 1)
    assert s.simplify(jet.second-jet.local_second) != s.zeros(3, 1)
    assert jet.local_second != s.zeros(3, 1)
    assert s.simplify(jet.first-s.I*k.cross(moments.spin)/2) != s.zeros(3, 1)


def test_uniform_boost_keeps_finite_size_mass_quadrupole():
    # An exact uniform Euler flow with a transported two-parcel marker.
    t = s.Symbol('t', real=True)
    v = s.Matrix([1, 2, 0])
    points = [s.Matrix([-1, 0, 0])+t*v, s.Matrix([1, 0, 0])+t*v]
    result = material_momentum_second_jet([1, 1], points, [v, v], [1, 0, 0])
    assert result.first == s.zeros(3, 1)
    assert result.local_second == -v
    # Replacing absolute v by v-V would incorrectly erase this term.
    assert result.second == -v


@pytest.mark.parametrize('parity', [0, 1])
def test_parity_controls_from_actual_angular_transport_and_finite_band_moments(parity):
    theta, t = s.symbols('theta t', real=True)
    omega = s.Symbol('omega', positive=True)
    base = passive_packet(1, omega, theta, t)
    packet = base.euler_velocity
    if parity == 0:
        packet = packet.subs(theta, theta-s.pi/2)
    observed = s.integrate(s.sin(theta)*packet, (theta, 0, 2*s.pi))/s.pi
    expected = s.cos(omega*t) if parity == 0 else s.sin(omega*t)
    assert s.trigsimp(observed-expected) == 0
    # Finite-width positive band fixtures; exact integrals, not point bands.
    bands = [(s.Rational(1, 10), s.Rational(1, 8)),
             (s.Rational(1, 5), s.Rational(1, 4)),
             (s.Rational(3, 10), s.Rational(3, 8))]
    moments = s.Matrix([[s.integrate(omega**(2*j+parity), (omega, lo, hi))/(hi-lo)
                         for lo, hi in bands] for j in range(3)])
    coeff = [s.Rational(2, 3), -1, s.Rational(1, 7)]
    result = passive_output_weights(moments, coeff, parity=parity)
    actual = sum(w*s.integrate(expected, (omega, lo, hi))/(hi-lo)
                 for w, (lo, hi) in zip(result.weights, bands))
    for j, c in enumerate(coeff):
        power = 2*j+parity
        assert s.simplify(s.limit(s.diff(actual, t, power), t, 0)
                          -s.factorial(power)*c) == 0
    assert result.residuals == s.zeros(3, 1)
    assert s.simplify(actual.subs(t, -t)-(-1)**parity*actual) == 0


def test_degenerate_output_band_matrix_is_exposed():
    with pytest.raises(ValueError, match='nonsingular'):
        passive_output_weights([[1, 1], [2, 2]], [0, 1], parity=0)
    with pytest.raises(ValueError, match='parity'):
        passive_output_weights([[1]], [1], parity=2)
