"""Actual QE Cartesian ionic inertia in Rydberg atomic-time units.

QE7.5 PW/src/dynamics_module.f90:867-875 uses amass(ityp)*amu_ry.
Modules/constants.f90:44-61 defines amu_ry=amu_si/(2*electronmass_si).
The supplied Native scalar, not a guessed species/isotope mass, is authoritative.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .inertial_fermion_dynamics import ionic_masses_ry
from .material_modes import AMU_RY
from .qe_response_bundle import QEResponseBundle


def ionic_mass_diagonal(bundle: QEResponseBundle) -> NDArray[np.float64]:
    """Return atom-major Cartesian M for M*u''=-gradient, tau measured in hbar/Ry."""
    background = bundle.background
    metadata = background.metadata
    nat, nsp = metadata["nat"], metadata["nsp"]
    if type(nat) is not int or type(nsp) is not int or nat < 1 or nsp < 1:
        raise ValueError("positive Native atom/species dimensions required")
    species = background.array("amass_species", (nsp,), "float64")
    types = background.array("ityp", (nat,), "int32")
    if background.descriptors["ityp"].index_base != 1:
        raise ValueError("actual Native one-based species mapping is invalid")
    conversion = metadata["amu_ry"]
    if isinstance(conversion, bool) or not np.isfinite(conversion) or conversion <= 0:
        raise ValueError("positive finite Native amu_ry is required")
    if not np.isclose(conversion, AMU_RY, rtol=8 * np.finfo(float).eps, atol=0):
        raise ValueError("actual Native Rydberg inertia convention differs from the source-bound QE constants")
    return ionic_masses_ry(species, types, conversion)
