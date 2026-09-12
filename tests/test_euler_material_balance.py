"""Exposing exact-calculus checks for P253/0001, not numerical PDE tests."""

import pytest
import sympy as sp

from substrate_framework.euler_material_balance import (
    material_tag_balance,
    radial_swirl_exterior_pressure,
)


def test_material_rates_from_independently_differentiated_trajectories():
    t = sp.symbols("t", real=True)
    masses = [1, 2, 3]
    points = [sp.Matrix([1+t+t*t, t*t, 0]),
              sp.Matrix([-1+t, 2-t+t*t, t]),
              sp.Matrix([t*t, -1+2*t, 1-t*t])]
    velocities = [x.diff(t) for x in points]
    rho = sp.Rational(7, 3)
    gradients = [-rho*x.diff(t, 2) for x in points]
    at_zero = lambda rows: [x.subs(t, 0) for x in rows]
    result = material_tag_balance(masses, at_zero(points), at_zero(velocities),
                                  at_zero(gradients), density=rho)
    mass = sum(masses)
    center = sum((m*x for m, x in zip(masses, points)), sp.zeros(3, 1))/mass
    relative_x = [x-center for x in points]
    relative_v = [x.diff(t) for x in relative_x]
    second = sum((m*r*r.T for m, r in zip(masses, relative_x)), sp.zeros(3))
    first = sum((m*c*r.T for m, c, r in zip(masses, relative_v, relative_x)), sp.zeros(3))
    kinetic = sum((m*c*c.T for m, c in zip(masses, relative_v)), sp.zeros(3))
    spin = sum((m*r.cross(c) for m, r, c in zip(masses, relative_x, relative_v)), sp.zeros(3, 1))
    energy = sum(m*v.dot(v)/2 for m, v in zip(masses, velocities))
    pairs = [(result.centroid_acceleration, center.diff(t, 2)),
             (result.first_momentum_rate, first.diff(t)),
             (result.second_mass_acceleration, second.diff(t, 2)),
             (result.spin_rate, spin.diff(t)),
             (result.kinetic_covariance_rate, kinetic.diff(t))]
    for actual, expected in pairs:
        assert sp.simplify(actual-expected.subs(t, 0)) == sp.zeros(*actual.shape)
    assert result.internal_energy_rate == sp.diff(sp.trace(kinetic)/2, t).subs(t, 0)
    assert result.total_energy_rate == energy.diff(t).subs(t, 0)
    assert result.force != sp.zeros(3, 1)
    assert result.spin_rate != sp.zeros(3, 1)


def test_uniform_pressure_acceleration_changes_translation_only():
    result = material_tag_balance([1, 2], [[-1, 0, 0], [2, 1, 0]],
                                  [[0, 0, 0], [0, 0, 0]],
                                  [[3, -6, 9], [3, -6, 9]], density=3)
    assert result.centroid_acceleration == sp.Matrix([-1, 2, -3])
    assert result.second_mass_acceleration == sp.zeros(3)
    assert result.spin_rate == sp.zeros(3, 1)
    assert result.internal_energy_rate == 0


def test_exterior_pressure_from_cartesian_source_and_radial_green_integrals():
    x, y, z, s, mu = sp.symbols("x y z s mu", real=True)
    radius2 = x*x+y*y+z*z
    f = (1-radius2)**3
    u = sp.Matrix([-y*f, x*f, 0])
    derivative = u.jacobian([x, y, z])
    assert sp.trace(derivative) == 0
    source = sp.expand(sp.trace(derivative*derivative))
    angular = sp.expand(source.subs({y: 0, x: s*sp.sqrt(1-mu*mu), z: s*mu}))
    p2 = (3*mu*mu-1)/2
    source0 = sp.integrate(angular, (mu, -1, 1))/2
    source2 = 5*sp.integrate(angular*p2, (mu, -1, 1))/2
    assert sp.simplify(angular-source0-source2*p2) == 0
    monopole = sp.integrate(s*s*source0, (s, 0, 1))
    quadrupole = sp.integrate(s**4*source2, (s, 0, 1))/5
    moment = sp.integrate(s**4*(1-s*s)**6, (s, 0, 1))
    assert monopole == 0
    assert quadrupole == -2*moment/3
    R = sp.symbols("R", positive=True)
    actual = radial_swirl_exterior_pressure([0, 0, R], [0, 0, 0], [0, 0, 1],
                                           moment, density=7, support_radius=1)
    assert sp.simplify(actual-7*quadrupole/R**3) == 0


def test_exact_ambient_centroid_and_shape_acceleration_and_no_pressure_monopole():
    x, y, z = sp.symbols("x y z", real=True)
    d, rho, J = sp.symbols("d rho J", positive=True)
    pressure = radial_swirl_exterior_pressure([x, y, z], [0, 0, d], [0, 0, 1],
                                             J, density=rho, support_radius=1)
    gradient = sp.Matrix([pressure.diff(q) for q in [x, y, z]])
    hessian = sp.hessian(pressure, [x, y, z])
    center = {x: 0, y: 0, z: 0}
    assert sp.simplify(-gradient.subs(center)/rho) == sp.Matrix([0, 0, 2*J/d**4])
    assert sp.simplify(hessian.subs(center)) == rho*J/d**5*sp.diag(4, 4, -8)
    assert sp.simplify(sp.trace(hessian)) == 0
    # Reversing circulation changes u's sign but leaves this quadratic
    # pressure unchanged; it cannot be read as a signed Coulomb charge.
    reflected = radial_swirl_exterior_pressure([x, y, z], [0, 0, d], [0, 0, -1],
                                               J, density=rho, support_radius=1)
    assert sp.simplify(pressure-reflected) == 0


def test_material_balance_input_contracts():
    with pytest.raises(ValueError, match="density"):
        material_tag_balance([1], [[0, 0, 0]], [[0, 0, 0]], [[0, 0, 0]], density=0)
    with pytest.raises(ValueError, match="equal length"):
        material_tag_balance([1], [[0, 0, 0]], [[0, 0, 0]], [], density=1)
    with pytest.raises(ValueError, match="outside"):
        radial_swirl_exterior_pressure([0, 0, 0], [0, 0, 0], [0, 0, 1],
                                       1, density=1, support_radius=1)
    with pytest.raises(ValueError, match="unit"):
        radial_swirl_exterior_pressure([0, 0, 2], [0, 0, 0], [0, 0, 2],
                                       1, density=1, support_radius=1)
