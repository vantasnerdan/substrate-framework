import pytest
import sympy as sp

from substrate_framework.euler_scale_causality import (
    axisymmetric_swirl_pressure_quadrupole,
    euler_similarity_weights,
    pressure_quadrupole,
)


def test_fixed_density_euler_similarity_weights_and_topology():
    a, b = sp.symbols("a b", positive=True)
    weights = euler_similarity_weights(a, b)
    assert weights.energy == a**2 / b**3
    assert weights.vorticity_impulse == a / b**3
    assert weights.angular_momentum == weights.kks_action == a / b**4
    assert weights.helicity == a**2 / b**2
    assert weights.circulation == a / b
    assert weights.tag_inertia == b**-5
    assert weights.topological_charge == 1


def test_continuous_similarity_rescales_action_inside_one_topology_class():
    target, action = sp.symbols("target action", positive=True)
    weights = euler_similarity_weights(target / action, sp.Integer(1))
    assert sp.simplify(weights.kks_action * action - target) == 0
    assert weights.topological_charge == 1


def test_axisymmetric_compact_swirl_has_nonzero_pressure_tail():
    m, rho, s, z = sp.symbols("m rho s z", positive=True)
    pressure = axisymmetric_swirl_pressure_quadrupole(m, s, z, density=rho)
    radius2 = s**2 + z**2
    expected = rho * m * (radius2 - 3 * z**2) / (
        4 * sp.pi * radius2 ** sp.Rational(5, 2)
    )
    assert sp.simplify(pressure - expected) == 0
    on_axis = sp.simplify(pressure.subs(s, 0))
    assert on_axis == -rho * m / (2 * sp.pi * z**3)
    assert sp.diff(on_axis, z) == 3 * rho * m / (2 * sp.pi * z**4)


def test_general_pressure_quadrupole_uses_euler_sign_and_cross_contraction():
    x, y, z, rho = sp.symbols("x y z rho", positive=True)
    moment = sp.Matrix([[2, 1, 0], [1, 3, 0], [0, 0, 5]])
    point = sp.Matrix([x, y, z])
    radius = sp.sqrt(x**2 + y**2 + z**2)
    green = 1 / (4 * sp.pi * radius)
    hessian = sp.hessian(green, (x, y, z))
    expected = rho * sum(
        moment[i, j] * hessian[i, j] for i in range(3) for j in range(3)
    )
    assert sp.simplify(pressure_quadrupole(moment, point, density=rho) - expected) == 0


@pytest.mark.parametrize("a,b", [(0, 1), (1, 0), (-1, 1)])
def test_similarity_rejects_nonpositive_scales(a, b):
    with pytest.raises(ValueError, match="positive"):
        euler_similarity_weights(a, b)
