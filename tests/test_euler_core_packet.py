"""Independent integral and differential checks for packet coefficient APIs."""

import pytest
import sympy as s

from substrate_framework.euler_core_packet import (
    common_circle_angular_rule,
    common_circle_moment_weights,
    gaussian_carrier_filter,
    laguerre_packet_angle,
    packet_material_moment_rows,
)


@pytest.mark.parametrize("n", [0, 1, 2, 8])
def test_physical_angle_from_independent_laplace_derivatives(n):
    w, t = s.symbols("w t", positive=True)
    # Derive the x**2 integral by differentiating the known x-weighted
    # Laguerre Laplace transform, independently of polynomial integration.
    transform = (n+1)*(t-1)**n/t**(n+2)
    raw = -s.diff(transform, t).subs(t, w/2)
    angle = laguerre_packet_angle(n, w)
    assert s.factor(angle-raw/raw.subs(w, 1)) == 0
    assert angle.subs(w, 1) == 1
    assert s.simplify(s.diff(angle, w).subs(w, 1)+2*n+2+s.Rational(1, n+1)) == 0
    assert s.simplify(angle-raw/(2*raw.subs(w, 1))) != 0


def test_finite_packet_gaussian_filter_and_full_action_measure():
    length, marker = s.symbols("L ell", positive=True)
    p, p0, q = s.symbols("p p0 q", real=True)
    result = gaussian_carrier_filter(length, marker, p, p0)
    cost = length**2*(q-p)**2+marker**2*(q-p0)**2
    completed = (q-result.center)**2/result.variance+result.envelope_cost*(p-p0)**2
    assert s.factor(cost-completed) == 0
    assert s.diff(result.center, p) == result.carrier_factor
    assert s.simplify(result.amplitude.subs(p, p0)-length/s.sqrt(length**2+marker**2)) == 0
    direct_norm = s.integrate(length**2*s.exp(-length**2*q*q), (q, -s.oo, s.oo))
    assert result.plancherel_weight == direct_norm
    assert s.simplify(s.diff(result.amplitude, p, 2).subs(p, p0)
                      +result.envelope_cost*result.amplitude.subs(p, p0)) == 0
    # Omitting the real envelope leaves a wrong nonzero second jet.
    assert s.diff(result.amplitude, p, 2).subs(p, p0) != 0


def test_radial_rows_derive_pressure_and_dilation_and_expose_rank_loss():
    x = s.Symbol("x", real=True)
    rows = packet_material_moment_rows(2, 1, x, reference_order=1)
    assert rows[:2] == (1, x)
    laguerre = s.assoc_laguerre(2, 1, x)
    expected = s.exp(-x/2)*(laguerre-2*s.diff(laguerre, x))
    for j in range(2):
        value, dilation = rows[2+2*j:4+2*j]
        assert s.simplify(value-x**j*expected) == 0
        assert s.simplify(dilation+value/2-s.Rational(3, 2)*x*s.diff(value, x)) == 0
    matrix = s.Matrix([[s.diff(row, x, j).subs(x, 0) for j in range(len(rows))]
                       for row in rows])
    assert matrix.det() != 0
    degenerate = packet_material_moment_rows(0, 1, x, reference_order=1)
    wrong = s.Matrix([[s.diff(row, x, j).subs(x, 0) for j in range(len(degenerate))]
                      for row in degenerate])
    assert wrong.det() == 0


@pytest.mark.parametrize("size", [1, 2, 3])
def test_positive_angular_rule_against_exact_arcsine_integrals(size):
    rule = common_circle_angular_rule(size)
    theta = s.Symbol("theta", real=True)
    for degree in range(2*size):
        moment = sum(weight*node**degree for node, weight in zip(rule.nodes, rule.weights))
        direct = s.integrate(s.cos(theta)**(2*degree), (theta, 0, s.pi))/s.pi
        assert s.simplify(moment-direct) == 0
    assert all(weight.is_positive for weight in rule.weights)
    assert common_circle_moment_weights(rule.nodes).weights == rule.weights


def test_actual_moment_solution_exposes_nonpositive_off_rule_weights():
    nodes = (s.Rational(1, 10), s.Rational(1, 5))
    rule = common_circle_moment_weights(nodes)
    assert sum(rule.weights) == 1
    assert sum(n*w for n, w in zip(nodes, rule.weights)) == s.Rational(1, 2)
    assert rule.weights[0] < 0  # No positivity is manufactured by clipping.


def test_invalid_orders_domains_and_singular_nodes():
    x = s.Symbol("x")
    for bad in (-1, s.Rational(1, 2), s.nan, True):
        with pytest.raises(ValueError, match="integer"):
            laguerre_packet_angle(bad, 1)
    for bad in (0, s.nan, s.oo, s.zoo):
        with pytest.raises(ValueError, match="finite|nonzero"):
            laguerre_packet_angle(2, bad)
    for bad in (0, -1, s.I, s.nan, s.oo):
        with pytest.raises(ValueError):
            gaussian_carrier_filter(bad, 1, 2, 2)
    with pytest.raises(ValueError, match="real"):
        gaussian_carrier_filter(1, 1, s.I, 1)
    with pytest.raises(ValueError, match="symbol"):
        packet_material_moment_rows(2, 1, x+1)
    with pytest.raises(ValueError, match="positive"):
        common_circle_angular_rule(0)
    for points in ((), (0,), (1,), (s.I,), (s.Rational(1, 2),)*2):
        with pytest.raises(ValueError):
            common_circle_moment_weights(points)
