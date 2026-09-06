import pytest
import sympy as sp

from substrate_framework.euler_quantum_bridge import (
    exchange_characters,
    hopf_fr_phase,
    scaled_euler_angular_momentum,
    spin_orbit_quantization,
    uniform_euler_frequency,
)


def test_spin_sphere_integral_class_and_half_integer_representation():
    hbar = sp.symbols("hbar", positive=True)
    data = spin_orbit_quantization(hbar / 2, hbar)
    assert data.cern_number == 1
    assert data.hilbert_dimension == 2
    assert data.spin == sp.Rational(1, 2)
    assert data.rotation_phase == -1


def test_euler_scaling_continuously_changes_action_with_fixed_profile_topology():
    j, scale = sp.symbols("j scale", positive=True)
    assert scaled_euler_angular_momentum(j, scale, 1) == scale * j
    assert scaled_euler_angular_momentum(j, 1, scale) == j / scale**4


def test_exchange_topology_permits_both_characters_and_does_not_select_one():
    assert exchange_characters() == (1, -1)
    assert hopf_fr_phase(1, deck_character=-1) == -1
    assert hopf_fr_phase(2, deck_character=-1) == 1
    assert hopf_fr_phase(3, deck_character=1) == 1


def test_bare_uniform_euler_has_one_convective_frequency_for_both_polarizations():
    k1, k2, k3, u1, u2, u3 = sp.symbols("k1 k2 k3 u1 u2 u3", real=True)
    assert uniform_euler_frequency((k1, k2, k3), (u1, u2, u3)) == (
        k1 * u1 + k2 * u2 + k3 * u3
    )


def test_quantization_and_character_domains_are_exposing():
    with pytest.raises(ValueError, match="positive"):
        spin_orbit_quantization(0, 1)
    with pytest.raises(ValueError, match="integer"):
        spin_orbit_quantization(sp.Rational(1, 3), 1)
    with pytest.raises(ValueError, match="deck_character"):
        hopf_fr_phase(1, deck_character=sp.I)
