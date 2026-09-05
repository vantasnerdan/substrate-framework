import pytest
import sympy as s

from substrate_framework.euler_time_control import finite_band_polynomial


def test_actual_exponential_signal_converges_with_time_derivatives_away_from_zero_frequency():
    h = s.Symbol('h', positive=True)
    t = s.Symbol('t', real=True)
    carrier = s.Rational(3, 2)
    result = finite_band_polynomial([1, -2, 3], carrier, h)
    signal = sum(w*s.exp(s.I*frequency*t)
                 for w, frequency in zip(result.weights, result.frequencies))
    target = s.exp(s.I*carrier*t)*(1-2*t+3*t*t)
    expansion = s.simplify(s.series(signal, h, 0, 2).removeO())
    for derivative in range(3):
        assert s.simplify(s.diff(expansion, t, derivative).subs(h, 0)
                          -s.diff(target, t, derivative)) == 0
    assert s.simplify(s.diff(expansion, h)) != 0
    # This fixed interval is separated from zero; no low-frequency limit occurs.
    assert all(1 < frequency.subs(h, s.Rational(1, 20)) < 2
               for frequency in result.frequencies)
    assert s.limit(result.coefficient_norm_bound, h, 0, dir='+') == s.oo


def test_small_unobserved_residual_can_give_order_one_normalized_output_error():
    eta = s.Symbol('eta', positive=True)
    omega, t = s.symbols('omega t', real=True)
    generator = s.I*omega*s.eye(3)+s.Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    state = s.Matrix([eta, eta, 1])
    observation = s.Matrix([[1, 0, 0]])
    residual = (generator-s.I*omega*s.eye(3))*state
    gain = (observation*state)[0]
    assert state.subs(eta, 0) != s.zeros(3, 1)
    assert residual.subs(eta, 0) == s.zeros(3, 1)
    actual = (observation*(t*generator).exp()*state)[0]/gain
    assert s.simplify(actual-s.exp(s.I*omega*t)) == t*s.exp(s.I*omega*t)
    # The residual/gain criterion exposes the error before claiming a clock.
    assert s.simplify(residual/gain) == s.Matrix([1, 0, 0])


def test_frequency_step_and_polynomial_domain():
    with pytest.raises(ValueError, match='nonzero'):
        finite_band_polynomial([1, 2], 1, 0)
    with pytest.raises(ValueError, match='coefficient'):
        finite_band_polynomial([], 1, 1)
    with pytest.raises(ValueError, match='real'):
        finite_band_polynomial([1], s.I, 1)
