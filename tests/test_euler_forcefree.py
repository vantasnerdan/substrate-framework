"""Independent field residuals and domain probes for the Bernoulli lift."""

import pytest
import sympy as s

from substrate_framework.euler_forcefree import planar_bernoulli_lift


def test_periodic_cell_full_cartesian_euler_and_forcefree_residuals():
    x, y, z = s.symbols("x y z", real=True)
    rho = s.Symbol("rho", positive=True)
    v = [s.sin(x)*s.cos(y), -s.cos(x)*s.sin(y)]
    p = rho*(s.cos(2*x)+s.cos(2*y))/4
    lift = planar_bernoulli_lift(v, p, rho, (x, y), 1)
    u = lift.velocity
    curl = s.Matrix([s.diff(u[2], y), -s.diff(u[2], x),
                     s.diff(u[1], x)-s.diff(u[0], y)])
    assert s.simplify(curl-lift.curl_factor*u) == s.zeros(3, 1)
    assert s.simplify(rho*u.jacobian((x, y, z))*u
                      +s.Matrix([s.diff(p, c) for c in (x, y, z)])) == s.zeros(3, 1)
    assert s.trigsimp(lift.axial_speed**2-1-2*s.sin(x)**2*s.sin(y)**2) == 0
    assert s.simplify(p/rho+u.dot(u)/2) == 1
    assert s.simplify(s.diff(lift.curl_factor, x)) != 0


def test_shaped_column_agrees_with_independently_integrated_axial_ode():
    x, y = s.symbols("x y", real=True)
    a, om, axial = s.symbols("a Omega U0", positive=True)
    r2 = x*x+y*y
    v = [-om*a*a*y/(a*a+r2), om*a*a*x/(a*a+r2)]
    p = om**2*a*a*r2/(2*(a*a+r2))
    result = planar_bernoulli_lift(v, p, 1, (x, y), axial**2/2)
    expected = axial**2-om**2*a*a+om**2*a**6/(a*a+r2)**2
    assert s.simplify(result.axial_speed**2-expected) == 0
    assert s.simplify(result.curl_factor*result.axial_speed
                      -2*om*a**4/(a*a+r2)**2) == 0
    # A caller still owes U0 > Omega*a for a positive whole-plane lift.
    assert result.domain_condition != s.true


def test_pressure_gauge_changes_no_field_when_level_is_shifted_consistently():
    x, y = s.symbols("x y", real=True)
    rho, c = s.symbols("rho C", positive=True)
    gauge = s.Symbol("gauge", real=True)
    first = planar_bernoulli_lift([y, 0], 0, rho, (x, y), c)
    shifted = planar_bernoulli_lift([y, 0], rho*gauge, rho, (x, y), c+gauge)
    assert first.velocity == shifted.velocity
    assert first.curl_factor == shifted.curl_factor
    assert s.simplify(first.planar_bernoulli+gauge-shifted.planar_bernoulli) == 0


def test_full_planar_subsystem_and_passive_axial_velocity_are_exact():
    x, y, z, t = s.symbols("x y z t", real=True)
    vx, vy, w, p = [s.Function(name)(t, x, y) for name in ("v_x", "v_y", "w", "p")]
    v, u = s.Matrix([vx, vy]), s.Matrix([vx, vy, w])
    euler = u.diff(t)+u.jacobian((x, y, z))*u+s.Matrix([s.diff(p, c) for c in (x, y, z)])
    planar = v.diff(t)+v.jacobian((x, y))*v+s.Matrix([s.diff(p, c) for c in (x, y)])
    assert euler[:2, 0] == planar
    assert s.expand(euler[2]-s.diff(w, t)-vx*s.diff(w, x)-vy*s.diff(w, y)) == 0
    # Exact local axial-energy conservation under the actual passive equation.
    density = w*w/2
    continuity = s.diff(density, t)+s.diff(vx*density, x)+s.diff(vy*density, y)
    assert s.expand(continuity-w*euler[2]
                    -density*(s.diff(vx, x)+s.diff(vy, y))) == 0


def test_harmonic_flow_is_included_without_inventing_vorticity():
    x, y = s.symbols("x y", real=True)
    lift = planar_bernoulli_lift([1, 0], 0, 1, (x, y), 1)
    assert lift.velocity == s.Matrix([1, 0, 1])
    assert lift.curl_factor == 0
    assert lift.domain_condition == s.true


def test_wrong_equations_and_invalid_domain_inputs_are_rejected():
    x, y, z = s.symbols("x y z", real=True)
    with pytest.raises(ValueError, match="residual"):
        planar_bernoulli_lift([-y, x], -(x*x+y*y)/2, 1, (x, y), 1)
    with pytest.raises(ValueError, match="residual"):
        planar_bernoulli_lift([x, y], 0, 1, (x, y), 1)
    for rho in (0, -1, s.oo, s.nan):
        with pytest.raises(ValueError, match="density"):
            planar_bernoulli_lift([0, 0], 0, rho, (x, y), 1)
    for c in (0, -1, s.oo, s.I, s.nan, s.zoo):
        with pytest.raises(ValueError, match="Bernoulli"):
            planar_bernoulli_lift([0, 0], 0, 1, (x, y), c)
    with pytest.raises(ValueError, match="spatial constants"):
        planar_bernoulli_lift([0, 0], 0, 1, (x, y), x)
    for coords in ((x, x), (x, y, z), (x, x+y)):
        with pytest.raises(ValueError, match="coordinates|coordinate"):
            planar_bernoulli_lift([0, 0], 0, 1, coords, 1)
    with pytest.raises(ValueError, match="planar column"):
        planar_bernoulli_lift([0, 0, 0], 0, 1, (x, y), 1)
