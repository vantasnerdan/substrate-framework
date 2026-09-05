import pytest
import sympy as sp

from substrate_framework.euler_compact_ring import compact_ring_fields, isotropic_tag_normalization
from substrate_framework.euler_joint import prepared_joint_symbol


def test_solid_rotation_residual_and_cylindrical_factors():
    r, z, om = sp.symbols("r z om", nonzero=True)
    # psi=0, I=om*r**2 gives u_phi=om*r; p=om**2*r**2/2.
    out = compact_ring_fields(0, om*r**2, om**2*r**2/2, r, z)
    assert out.velocity == sp.ImmutableMatrix([0, om*r, 0])
    assert out.residual == sp.zeros(3, 1)
    wrong = compact_ring_fields(0, om*r**2, 0, r, z)
    assert wrong.residual[0] == -om**2*r


def test_poloidal_case_exposes_azimuthal_metric_term():
    r, z = sp.symbols("r z", positive=True)
    out = compact_ring_fields(r*z, r**2, 0, r, z)
    # ur=-1, up=r; both ur*d_r(up) and ur*up/r contribute.
    assert out.velocity == sp.ImmutableMatrix([-1, r, z/r])
    assert out.residual[1] == -2


def test_domain_and_positive_literal_normalization():
    r, z = sp.symbols("r z", positive=True)
    with pytest.raises(ValueError):
        compact_ring_fields(0, 1, 0, 0, z)
    j, a, ell = isotropic_tag_normalization(6, 4, 2, 3, density=2,
                                            mean_tag_fraction=sp.Rational(1, 2))
    assert (j, a, ell) == (sp.Rational(8, 9), sp.Rational(2, 3), 2)
    # Changing E cannot change inertia j; it changes only kinetic a.
    assert isotropic_tag_normalization(60, 4, 2, 3, density=2,
                                       mean_tag_fraction=sp.Rational(1, 2))[0] == j
    with pytest.raises(ValueError):
        isotropic_tag_normalization(1, 1, 0, 3)
    with pytest.raises(ValueError):
        isotropic_tag_normalization(1, 1, 1, 3, density=-1)
    with pytest.raises(ValueError):
        compact_ring_fields(0, 1, 0, r, r)


def test_inertia_from_an_independent_ball_integral():
    radius, polar, azimuth = sp.symbols("radius polar azimuth", real=True)
    radial_second = sp.integrate(radius**4, (radius, 0, 1))
    sphere_area = sp.integrate(sp.sin(polar), (polar, 0, sp.pi))*2*sp.pi
    # Unit ball's inertia about z comes directly from x^2+y^2.
    axial_inertia = radial_second*sp.integrate(sp.sin(polar)**3, (polar, 0, sp.pi))
    axial_inertia *= sp.integrate(1, (azimuth, 0, 2*sp.pi))
    tag_moment = radial_second*sphere_area
    tag_volume = sp.Rational(4, 3)*sp.pi
    fraction, density, cell = sp.Rational(1, 3), 5, 64
    j, _, _ = isotropic_tag_normalization(
        7, tag_moment, tag_volume, cell, density, fraction
    )
    assert sp.simplify(j-density*fraction*axial_inertia/cell) == 0


def test_measured_inputs_feed_the_same_canonical_joint_forms():
    k = sp.symbols("k", real=True)
    density, frequency = 2, 3
    j, a, length_squared = isotropic_tag_normalization(
        6, 4, 2, 3, density, sp.Rational(1, 2)
    )
    transverse = j*frequency**2/(4*density)+frequency**2*length_squared
    longitudinal = frequency**2*length_squared
    symbol = prepared_joint_symbol(k, density, j, frequency, a, transverse, longitudinal)
    assert symbol.mass == sp.diag(2, 2, sp.Rational(8, 9), sp.Rational(8, 9), sp.Rational(8, 9))
    for value in symbol.defect:
        assert all(sp.diff(value, k, order).subs(k, 0) == 0 for order in range(3))
    assert any(sp.diff(value, k, 3) != 0 for value in symbol.defect)
    zero, identity = sp.zeros(5), sp.eye(5)
    phase = sp.BlockMatrix([[zero, symbol.mass], [-symbol.mass, zero]]).as_explicit()
    generator = sp.BlockMatrix([
        [zero, identity], [-symbol.mass.inv()*symbol.stiffness, zero]
    ]).as_explicit()
    energy = sp.diag(symbol.stiffness, symbol.mass)
    assert sp.simplify(phase*generator+energy) == sp.zeros(10)
