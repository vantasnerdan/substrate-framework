"""Exact independent checks of the conditional smooth Euler displacement API."""

import pytest
import sympy as sp

from substrate_framework.euler_displacement import (
    euler_displacement_perturbation,
    euler_jacobi_density,
    material_derivative,
    stationarize_planar,
)


def test_general_transport_to_linearized_euler():
    x, y, t = sp.symbols("x y t", real=True)
    coords = (x, y)
    u = sp.Matrix([sp.Function(f"u{i}")(x, y) for i in range(2)])
    xi = sp.Matrix([sp.Function(f"q{i}")(t, x, y) for i in range(2)])
    v = euler_displacement_perturbation(xi, u, coords, t)
    lhs = material_derivative(v, u, coords, t) + u.jacobian(coords)*v
    acc = u.jacobian(coords)*u
    rhs = material_derivative(material_derivative(xi, u, coords, t), u, coords, t)
    rhs -= acc.jacobian(coords)*xi
    assert sp.simplify(lhs-rhs) == sp.zeros(2, 1)


def test_jacobi_variation_taylor_green_retains_gyro_and_pressure():
    x, y, t = sp.symbols("x y t", real=True)
    rho = sp.Symbol("rho", positive=True)
    coords = (x, y)
    u = sp.Matrix([sp.sin(x)*sp.cos(y), -sp.cos(x)*sp.sin(y)])
    p = rho*(sp.cos(2*x)+sp.cos(2*y))/4
    assert sp.simplify(rho*u.jacobian(coords)*u + sp.Matrix([sp.diff(p, c) for c in coords])) == sp.zeros(2, 1)
    q = [sp.Function(f"q{i}")(t, x, y) for i in range(2)]
    density = euler_jacobi_density(q, u, p, rho, coords, t)
    for i in range(2):
        el = -sp.diff(density, q[i])
        for c in (t, x, y):
            el += sp.diff(sp.diff(density, sp.diff(q[i], c)), c)
        target = rho*material_derivative(material_derivative(q, u, coords, t), u, coords, t)
        target += sp.hessian(p, coords)*sp.Matrix(q)
        assert sp.simplify(el-target[i]) == 0
    static_xi = sp.Matrix([1, 0])
    assert sp.simplify(euler_jacobi_density(static_xi, u, p, rho, coords, t)-rho*sp.cos(2*x)/2) == 0
    # A wrong pressure sign reverses a nonzero acceleration coefficient.
    assert sp.simplify((2*sp.hessian(p, coords)*static_xi)[0]) != 0


def test_planar_stationarization_general_residual_identity():
    x, y = sp.symbols("x y", real=True)
    om, rho = sp.symbols("Omega rho", positive=True)
    psi, p = sp.Function("psi")(x, y), sp.Function("p")(x, y)
    coords = (x, y)
    v = sp.Matrix([-sp.diff(psi, y), sp.diff(psi, x)])
    j = sp.Matrix([[0, -1], [1, 0]])
    w, pw = stationarize_planar(psi, p, om, rho, coords)
    rotating = v.jacobian(coords)*(v-om*j*sp.Matrix(coords)) + om*j*v
    rotating += sp.Matrix([sp.diff(p, c) for c in coords])/rho
    stationary = w.jacobian(coords)*w + sp.Matrix([sp.diff(pw, c) for c in coords])/rho
    assert sp.simplify(stationary-rotating) == sp.zeros(2, 1)
    assert sp.simplify(sp.diff(w[1]-v[1], x)-sp.diff(w[0]-v[0], y)+2*om) == 0


def test_rigid_rotation_becomes_zero_flow_constant_pressure():
    x, y = sp.symbols("x y", real=True)
    om, rho = sp.symbols("Omega rho", positive=True)
    psi = om*(x*x+y*y)/2
    w, p = stationarize_planar(psi, rho*om**2*(x*x+y*y)/2, om, rho, (x, y))
    assert w == sp.zeros(2, 1)
    assert sp.simplify(p) == 0


def test_exact_api_rejects_wrong_dimension_and_density():
    x, y, z, t = sp.symbols("x y z t")
    with pytest.raises(ValueError):
        stationarize_planar(x*y, 0, 1, 1, (x, y, z))
    with pytest.raises(ValueError):
        euler_jacobi_density([x, y], [0, 0], 0, -1, (x, y), t)
    with pytest.raises(ValueError):
        material_derivative([x, y], [0], (x, y), t)
