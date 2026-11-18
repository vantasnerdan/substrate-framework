from math import pi, sqrt

import pytest
from scipy.constants import elementary_charge, epsilon_0, hbar, m_e

from substrate_framework.cluster_screening import (
    local_tf_shift_eV, spectator_shift_eV, thomas_fermi_contact_eV,
)


def test_closer_positive_spectator_is_not_screening_and_reference_cancels():
    assert spectator_shift_eV(1.8, 2.75) < 0
    assert spectator_shift_eV(2.75, 2.75) == 0
    assert local_tf_shift_eV(.29468554754150655, .29468554754150655) == 0


@pytest.mark.parametrize('density_A3', [.016, .29468554754150655, 1.0])
def test_tf_contact_matches_independent_dos_response_in_si(density_A3):
    n = density_A3*1e30
    k_F = (3*pi*pi*n)**(1/3)
    E_F = hbar*hbar*k_F*k_F/(2*m_e)
    dos_per_J_m3 = 3*n/(2*E_F)
    kappa = sqrt(elementary_charge**2*dos_per_J_m3/epsilon_0)
    independent_eV = elementary_charge*kappa/(4*pi*epsilon_0)
    assert thomas_fermi_contact_eV(density_A3) == pytest.approx(independent_eV, rel=1e-9)
