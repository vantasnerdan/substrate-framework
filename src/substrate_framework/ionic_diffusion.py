"""Conditional moving-charge response of an infinite linear drift-diffusion medium.

Not an atomic-scale dielectric model. See canonical issue237 and its source-bound
construction; carrier charge magnitude is one elementary charge.
"""
from __future__ import annotations

from math import asinh, isfinite, pi, sqrt
from scipy.constants import elementary_charge, epsilon_0


def static_ionic_shift_eV(kappa_per_m: float, relative_permittivity: float = 1.0) -> float:
    """Additional ionic constant-potential reduction within this linear model."""
    if (not isfinite(kappa_per_m) or kappa_per_m < 0 or
            not isfinite(relative_permittivity) or relative_permittivity <= 0):
        raise ValueError('Require finite kappa>=0 and background permittivity>0')
    return elementary_charge * kappa_per_m / (4 * pi * epsilon_0 * relative_permittivity)


def diffusive_ionic_shift_eV(kappa_per_m: float, diffusivity_m2_s: float,
                              speed_m_s: float, relative_permittivity: float = 1.0) -> float:
    """Retarded shift: U_static asinh(s)/s, s=v/(2 D kappa).

    D is the collective diffusion coefficient in the linear response equation,
    not necessarily the tracer diffusivity. The zero-speed limit is equilibrium.
    """
    static = static_ionic_shift_eV(kappa_per_m, relative_permittivity)
    if (not isfinite(diffusivity_m2_s) or diffusivity_m2_s <= 0 or
            not isfinite(speed_m_s) or speed_m_s < 0):
        raise ValueError('Require finite collective D>0 and speed>=0')
    if kappa_per_m == 0 or speed_m_s == 0:
        return static
    ratio = speed_m_s / (2 * diffusivity_m2_s * kappa_per_m)
    return static * asinh(ratio) / ratio


def darken_charge_parameters(site_density_m3: float, occupancy: float,
                             site_capacity: float, tracer_diffusivity_m2_s: float,
                             thermal_energy_eV: float, mu_slope_eV: float,
                             relative_permittivity: float = 1.0) -> tuple[float, float]:
    """Return (kappa [m^-1], collective D [m^2/s]) from a stable Darken model.

    J=-nM D_tr x(1-x/xmax) grad(mu)/kT, mu in eV per particle.
    Frozen trapped ions are not part of the mobile occupancy or susceptibility.
    """
    values = (site_density_m3, occupancy, site_capacity, tracer_diffusivity_m2_s,
              thermal_energy_eV, mu_slope_eV, relative_permittivity)
    if (not all(map(isfinite, values)) or min(site_density_m3, site_capacity,
            tracer_diffusivity_m2_s, thermal_energy_eV, mu_slope_eV,
            relative_permittivity) <= 0 or not 0 < occupancy < site_capacity):
        raise ValueError('Require a stable finite mobile lattice-gas state')
    kappa = sqrt(elementary_charge * site_density_m3 /
                 (epsilon_0 * relative_permittivity * mu_slope_eV))
    collective = (tracer_diffusivity_m2_s * occupancy *
                  (1 - occupancy / site_capacity) * mu_slope_eV / thermal_energy_eV)
    return kappa, collective
