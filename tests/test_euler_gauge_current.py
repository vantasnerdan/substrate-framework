import pytest
import sympy as sp

from substrate_framework.euler_gauge_current import (
    coulomb_force,
    coulomb_pair_energy,
    maxwell_speed,
    transported_charge_current,
)


def test_transported_tag_current_keeps_signed_density_and_velocity():
    rho_q, current = transported_charge_current(3, -sp.Symbol("chi", positive=True), [1, 2, 0])
    assert rho_q == -3 * sp.Symbol("chi", positive=True)
    assert current == sp.Matrix([-3, -6, 0]) * sp.Symbol("chi", positive=True)


def test_coulomb_energy_and_force_have_electric_sign():
    assert coulomb_pair_energy(2, 3, 5, [0, 0, 2]) == sp.Rational(15, 16) / sp.pi
    assert coulomb_pair_energy(2, 3, -5, [0, 0, 2]) == -sp.Rational(15, 16) / sp.pi
    assert coulomb_force(2, 3, 5, [0, 0, 2]) == sp.Matrix([0, 0, sp.Rational(15, 32) / sp.pi])


def test_maxwell_speed_and_domains():
    assert maxwell_speed(2, 8) == sp.Rational(1, 4)
    with pytest.raises(ValueError):
        maxwell_speed(0, 1)
    with pytest.raises(ValueError):
        coulomb_pair_energy(1, 1, 1, [0, 0, 0])
    with pytest.raises(ValueError):
        transported_charge_current(1, 1, [1, 2])
