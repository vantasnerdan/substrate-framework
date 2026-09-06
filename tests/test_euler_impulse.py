"""Independent Newton-kernel and vortex-loop normalization checks."""

import pytest
import sympy as sp

from substrate_framework.euler_impulse import (
    impulse_dipole_cross_energy,
    impulse_dipole_velocity,
)


def test_velocity_is_curl_of_first_newton_moment_and_matches_circular_loop():
    x, y, z = sp.symbols("x y z", real=True)
    q = sp.Matrix([x, y, z])
    I = sp.Matrix([2, -3, 5])
    potential = I.cross(q)/(4*sp.pi*(q.dot(q))**sp.Rational(3, 2))
    curl = sp.Matrix([potential[2].diff(y)-potential[1].diff(z),
                      potential[0].diff(z)-potential[2].diff(x),
                      potential[1].diff(x)-potential[0].diff(y)])
    assert sp.simplify(curl-impulse_dipole_velocity(I, q)) == sp.zeros(3, 1)
    R, Gamma, d, epsilon = sp.symbols("R Gamma d epsilon", positive=True)
    # Direct circular line integral gives I_z=pi*Gamma*R^2 and
    # u_z=Gamma*R^2/[2*(R^2+d^2)^(3/2)]. The filament is a normalization
    # fixture, not a smooth-carrier existence theorem.
    exact_axis = Gamma*R**2/(2*(R**2+(d/epsilon)**2)**sp.Rational(3, 2))
    leading = impulse_dipole_velocity([0, 0, sp.pi*Gamma*R**2], [0, 0, d])[2]
    assert sp.limit(exact_axis/epsilon**3, epsilon, 0, dir="+") == leading


def test_cross_energy_from_second_newton_derivative_not_magnetic_potential_sign():
    x, y, z = sp.symbols("x y z", real=True)
    q = sp.Matrix([x, y, z])
    Ia, Ib = sp.Matrix([2, -3, 5]), sp.Matrix([-7, 11, 13])
    green = 1/(4*sp.pi*sp.sqrt(q.dot(q)))
    hessian = sp.hessian(green, [x, y, z])
    # M[j,i]=integral y_j*omega_i = epsilon[j,i,k] I[k].
    Ma = sp.Matrix(3, 3, lambda j, i: sum(sp.LeviCivita(j, i, k)*Ia[k] for k in range(3)))
    Mb = sp.Matrix(3, 3, lambda j, i: sum(sp.LeviCivita(j, i, k)*Ib[k] for k in range(3)))
    # The mixed y*z term in G(d+z-y) is -G_jk*y_j*z_k.
    source_energy = -sum(hessian[j, k]*Ma[j, i]*Mb[k, i]
                         for i in range(3) for j in range(3) for k in range(3))
    assert sp.simplify(impulse_dipole_cross_energy(Ia, Ib, q, density=3)-3*source_energy) == 0
    assert impulse_dipole_cross_energy([0, 0, 1], [0, 0, 1], [0, 0, 2], density=3) > 0
    assert impulse_dipole_cross_energy([0, 0, 1], [0, 0, 1], [2, 0, 0], density=3) < 0


def test_zero_impulse_and_invalid_domains():
    assert impulse_dipole_velocity([0, 0, 0], [1, 0, 0]) == sp.zeros(3, 1)
    with pytest.raises(ValueError, match="nonzero"):
        impulse_dipole_velocity([1, 0, 0], [0, 0, 0])
    with pytest.raises(ValueError, match="density"):
        impulse_dipole_cross_energy([1, 0, 0], [0, 1, 0], [1, 0, 0], density=-1)
