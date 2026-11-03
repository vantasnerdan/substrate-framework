"""Same-saved-density Ni/H ionic/Hartree/internal-PBE background, in Ry.

Sources: QE7.5 setlocal, v_of_rho, gradcorr, set_vrs and setup_dgc.
The backend owns Native XC thresholds and unit conversion; this layer owns
physical Fourier derivatives, negative divergence and distinct Coulomb G0.
No saved occupations, temperature, state preparation or eigenvalue fit enters.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping
import xml.etree.ElementTree as ET

import numpy as np

from .qe_material_grid import PhysicalRhoFFT
from .qe_nc_radial import build_radial_tables
from .qe_nc_upf import read_nc_upf
from .qe_saved_native import load_density, read_bound_bytes
from .qe_saved_reference import load_saved_reference
from .spin_pbe_kernel import apply_spin_pbe_arrays

COEFFICIENT_NAMES = ("dmuxc", "dvxc_rr", "dvxc_sr", "dvxc_ss", "dvxc_s")


@dataclass(frozen=True)
class MaterialInputs:
    reference: object
    density: object
    tables: Mapping
    sources: tuple


@dataclass(frozen=True)
class MaterialReference:
    inputs: MaterialInputs
    fft: PhysicalRhoFFT
    rho_total_magnetization: np.ndarray
    rho_updown: np.ndarray
    grho: np.ndarray
    ionic_ry: np.ndarray
    hartree_ry: np.ndarray
    local_xc_ry: np.ndarray
    gradient_flux_ry_bohr: np.ndarray
    xc_ry: np.ndarray
    vrs: np.ndarray
    coefficients: Mapping
    sources: tuple
    imaginary_maxima: Mapping

    def apply_kernel(self, delta_n, *, q_bohr_inv=(0., 0., 0.)):
        """Full longitudinal spin K_q on complex physical up/down envelopes."""
        return apply_spin_pbe_arrays(
            self.fft.at_q(q_bohr_inv), delta_n, grho=self.grho,
            coefficients=self.coefficients, tpiba=2*np.pi/self.inputs.reference.alat_bohr)


def load_material_inputs(xml_binding, density_binding, upf_bindings, *, qmax_bohr_inv=None):
    """Bind exact decoding buffers; build radial tables once per species.

    upf_bindings maps XML species names to path/bytes/SHA256 bindings.
    Default radial range covers density |G| and the declared PW cutoff.
    """
    reference = load_saved_reference(xml_binding)
    xml, _ = read_bound_bytes(xml_binding)
    functional = ET.fromstring(xml).findtext("output/dft/functional")
    if functional is None or functional.strip().upper() != "PBE":
        raise ValueError("saved reference must explicitly declare PBE")
    density = load_density(density_binding)
    if (reference.gamma_only or not reference.lsda or reference.noncolin or reference.spinorbit
            or reference.fft_grid != reference.fft_smooth or density.nspin != 2
            or density.gamma_only or density.ngm != reference.ngm):
        raise ValueError("requires full-complex collinear PBE and equal dense/smooth physical grids")
    if not np.allclose(density.reciprocal_bohr_inv, reference.reciprocal_bohr_inv, atol=2e-13, rtol=0):
        raise ValueError("saved density/XML physical reciprocal bases disagree")
    volume = abs(float(np.linalg.det(reference.cell_bohr)))
    momenta = density.reciprocal_bohr_inv @ density.miller
    required = max(float(np.linalg.norm(momenta, axis=0).max()), np.sqrt(2*reference.ecutwfc_ha))
    qmax = required if qmax_bohr_inv is None else float(qmax_bohr_inv)
    if qmax < required:
        raise ValueError("radial table range does not cover saved density/PW cutoff")
    tables = {}
    sources = [dict(reference.source), dict(density.source)]
    for species in reference.species:
        binding = upf_bindings[species.name]
        upf = read_nc_upf(binding["path"], expected_bytes=binding["bytes"], expected_sha256=binding["sha256"])
        if upf.element != species.name:
            raise ValueError("UPF/XML species identity mismatch")
        tables[species.name] = build_radial_tables(upf, volume, qmax)
        sources.append(dict(path=upf.path, bytes=str(len(upf.original_bytes)), sha256=upf.sha256))
    return MaterialInputs(reference, density, MappingProxyType(tables), tuple(sources))


def reconstruct_reference(inputs: MaterialInputs, backend) -> MaterialReference:
    """Reconstruct raw fields; retain imaginary roundoff diagnostics, never fit G0."""
    reference, density = inputs.reference, inputs.density
    volume = abs(float(np.linalg.det(reference.cell_bohr)))
    fft = PhysicalRhoFFT.create(reference.fft_grid, density.miller, density.reciprocal_bohr_inv, volume)
    imaginary = {}

    def real_field(name, value):
        imaginary[name] = float(np.max(np.abs(value.imag)))
        if imaginary[name] > 1e-10 * (1 + float(np.max(np.abs(value.real)))):
            raise ValueError(f"{name} is not a real saved-background field: {imaginary[name]}")
        return np.asfortranarray(value.real)

    rho = real_field("rho_total_magnetization", fft.from_g(density.components))
    updown = np.column_stack(((rho[:, 0] + rho[:, 1])/2, (rho[:, 0] - rho[:, 1])/2))
    grho = real_field("grho", fft.gradient(updown))
    g = fft.wavevectors
    moduli = np.linalg.norm(g, axis=0)
    ionic_g = np.zeros(density.ngm, dtype=np.complex128)
    for species in reference.species:
        # Radial interpolation once per species, not once per atom.
        local = inputs.tables[species.name].vloc(moduli)
        structure = np.zeros(density.ngm, dtype=np.complex128)
        for atom, name in enumerate(reference.atom_names):
            if name == species.name:
                structure += np.exp(-1j * (g.T @ reference.positions_bohr[:, atom]))
        ionic_g += local * structure
    ionic = real_field("ionic_ry", fft.from_g(ionic_g[:, None]))[:, 0]
    hartree_g = np.zeros(density.ngm, dtype=np.complex128)
    nonzero = np.any(density.miller != 0, axis=0)
    hartree_g[nonzero] = 8*np.pi*density.components[nonzero, 0]/moduli[nonzero]**2
    hartree = real_field("hartree_ry", fft.from_g(hartree_g[:, None]))[:, 0]
    evaluated = backend.evaluate(rho, grho)
    divergence = real_field("xc_flux_divergence", fft.divergence(evaluated.gradient_flux_ry_bohr))
    xc = evaluated.local_xc_ry - divergence
    coefficients = MappingProxyType({name: getattr(evaluated, name) for name in COEFFICIENT_NAMES})
    return MaterialReference(inputs, fft, rho, updown, grho, ionic, hartree,
                             evaluated.local_xc_ry, evaluated.gradient_flux_ry_bohr,
                             xc, ionic[:, None] + hartree[:, None] + xc,
                             coefficients, inputs.sources + (dict(backend.source),)
                             + tuple(dict(source) for source in backend.dependencies),
                             MappingProxyType(imaginary))
