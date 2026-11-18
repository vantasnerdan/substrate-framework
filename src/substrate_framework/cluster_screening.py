"""Conditional D-D environment shifts, not a first-principles defect dielectric.

Nuclear spectator and local-density TF are alternative closures, not additive
measurements. Both return changes relative to a matched reference environment.
"""
from math import isfinite, pi, sqrt
from scipy.constants import elementary_charge, epsilon_0, physical_constants

COULOMB_EV_A = elementary_charge / (4 * pi * epsilon_0) * 1e10
BOHR_RADIUS_A = physical_constants['Bohr radius'][0] * 1e10


def spectator_shift_eV(distance_A: float, reference_distance_A: float) -> float:
    """Frozen bare-spectator leading contact shift relative to the reference.

    Initial target/spectator energy C/R becomes2C/R at nuclear coalescence;
    closer positive spectator charge therefore raises, not lowers, the barrier.
    Nuclear turning-point distances must be small compared with both distances.
    """
    if (not isfinite(distance_A) or not isfinite(reference_distance_A) or
            min(distance_A, reference_distance_A) <= 0):
        raise ValueError('Require finite positive spectator distances')
    return COULOMB_EV_A * (1/reference_distance_A-1/distance_A)


def thomas_fermi_contact_eV(valence_density_e_A3: float) -> float:
    """Uniform unpolarized free-electron TF contact value, bare electron mass.

    Applying this to a local PAW-valence midpoint is an explicit proxy; localized
    valence density is not automatically a conduction-electron susceptibility.
    """
    if not isfinite(valence_density_e_A3) or valence_density_e_A3 < 0:
        raise ValueError('Require finite nonnegative valence electron density')
    k_F = (3*pi*pi*valence_density_e_A3)**(1/3)
    k_TF = sqrt(4*k_F/(pi*BOHR_RADIUS_A))
    return COULOMB_EV_A*k_TF


def local_tf_shift_eV(density_e_A3: float, reference_density_e_A3: float) -> float:
    """Local-density TF proxy change; avoids adding reference screening twice."""
    return thomas_fermi_contact_eV(density_e_A3)-thomas_fermi_contact_eV(reference_density_e_A3)
