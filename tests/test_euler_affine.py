"""Physical curl and deformation checks, independent of coefficient matching."""

import pytest
import sympy as sp

from substrate_framework.euler_affine import affine_vorticity_energy


def test_actual_beltrami_tube_energy_and_affine_lattice():
    a, b, rho, scale = sp.symbols("a b rho scale", positive=True)
    modes = {(1, 0, 0): [0, -sp.I*a/2, a/2],
             (-1, 0, 0): [0, sp.I*a/2, a/2],
             (0, 1, 0): [sp.I*b/2, 0, b/2],
             (0, -1, 0): [-sp.I*b/2, 0, b/2]}
    assert affine_vorticity_energy(sp.eye(3), modes, rho) == rho*(a*a+b*b)/2
    f = sp.diag(scale, 1/scale, 1)
    actual = affine_vorticity_energy(f, modes, rho)
    expected = rho*(a*a*(1+scale**2)+b*b*(1+scale**-2))/4
    assert sp.simplify(actual-expected) == 0


def test_nondiagonal_affine_energy_by_direct_curl_inverse():
    shear, rho = sp.symbols("shear rho", real=True, positive=True)
    f = sp.Matrix([[1, shear, 0], [0, 1, 0], [0, 0, 1]])
    wave, omega = sp.Matrix([1, 2, 3]), sp.Matrix([2, -1, 0])
    kf, wf = f.inv().T*wave, f*omega
    velocity = sp.I*kf.cross(wf)/kf.dot(kf)
    direct = rho*sp.conjugate(velocity).dot(velocity)/2
    assert sp.simplify(affine_vorticity_energy(f, {tuple(wave): omega}, rho)-direct) == 0
    assert sp.simplify(sp.I*kf.cross(velocity)-wf) == sp.zeros(3, 1)


def test_rigid_rotation_and_empty_field():
    angle = sp.Symbol("angle", real=True)
    q = sp.Matrix([[sp.cos(angle), -sp.sin(angle), 0],
                   [sp.sin(angle), sp.cos(angle), 0], [0, 0, 1]])
    modes = {(1, 0, 0): [0, 1, sp.I]}
    assert affine_vorticity_energy(q, modes, 1) == affine_vorticity_energy(sp.eye(3), modes, 1)
    assert affine_vorticity_energy(q, {}, 1) == 0


def test_wrong_physical_inputs_are_rejected():
    with pytest.raises(ValueError, match="determinant"):
        affine_vorticity_energy(2*sp.eye(3), {}, 1)
    with pytest.raises(ValueError, match="solenoidal"):
        affine_vorticity_energy(sp.eye(3), {(1, 0, 0): [1, 0, 0]}, 1)
    with pytest.raises(ValueError, match="zero vorticity"):
        affine_vorticity_energy(sp.eye(3), {(0, 0, 0): [1, 0, 0]}, 1)
    with pytest.raises(ValueError, match="positive"):
        affine_vorticity_energy(sp.eye(3), {}, -1)
