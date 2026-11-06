"""Saved-reference frozen NC curvature extension, not screened phonon modes.

C(q)=C_A22(0)+C_II^QE(q)-C_II^QE(0). The complete imported electronic
remainder (including residual offsite blocks) is kept unchanged. It belongs
to the original numerical MV reference, not A27's prepared populations.
"""
from __future__ import annotations

from dataclasses import dataclass
import io
import json
import re

import numpy as np

from .material_force_lattice import validated_grid_shape, q_grid_to_force_lattice
from .qe_ionic_curvature import QEIonicCurvature, qe_alpha
from .qe_nc_upf import read_nc_upf
from .qe_saved_native import read_bound_bytes, load_density
from .qe_saved_reference import load_saved_reference


class FrozenLattice:
    def __init__(self, saved_gamma, ionic):
        self.saved_gamma = np.array(saved_gamma, dtype=complex, copy=True)
        self.ionic = ionic
        if self.saved_gamma.shape != (3 * ionic.nat, 3 * ionic.nat) or not np.isfinite(self.saved_gamma).all():
            raise ValueError("complete finite atom-major Gamma matrix required")
        self.ionic_gamma = ionic.matrix(np.zeros(3))
        self.electronic_remainder = self.saved_gamma - self.ionic_gamma

    def matrix(self, q_cart_bohr_inv):
        # Add the difference, not a repartitioned/repaired electronic matrix.
        return self.saved_gamma + (self.ionic.matrix(q_cart_bohr_inv) - self.ionic_gamma)

    def remainder_report(self):
        tensor = self.electronic_remainder.reshape(self.ionic.nat, 3, self.ionic.nat, 3)
        onsite = np.stack([tensor[a, :, a, :] for a in range(self.ionic.nat)])
        offsite = tensor.copy()
        for a in range(self.ionic.nat):
            offsite[a, :, a, :] = 0
        return dict(onsite_blocks_ry_per_bohr2=onsite,
                    onsite_frobenius=float(np.linalg.norm(onsite)),
                    offsite_frobenius=float(np.linalg.norm(offsite)),
                    offsite_max_abs=float(np.max(np.abs(offsite))),
                    hermiticity_frobenius=float(np.linalg.norm(self.electronic_remainder - self.electronic_remainder.conj().T)),
                    imaginary_frobenius=float(np.linalg.norm(self.electronic_remainder.imag)))

    def complete_q_grid(self, grid_shape):
        """Return every slot [n1,n2,n3,nat,3,nat,3], with no atomic tau phase.

        Representatives are centered with QE ties-away folding; finite G
        truncation need not be exactly periodic at the mesh boundary. No
        opposite-q averaging hides that residual. No screened-mode claim.
        """
        shape = validated_grid_shape(grid_shape)
        reciprocal = 2 * np.pi * np.linalg.inv(self.ionic.cell).T
        result = np.empty(shape + (self.ionic.nat, 3, self.ionic.nat, 3), complex)
        for index in np.ndindex(shape):
            q = np.asarray(index) / shape
            q -= np.floor(q + 0.5)
            result[index] = self.matrix(reciprocal @ q).reshape(result.shape[3:])
        return result

    def force_lattice(self, grid_shape):
        return q_grid_to_force_lattice(self.complete_q_grid(grid_shape))


@dataclass(frozen=True)
class MaterialFrozenLattice:
    lattice: FrozenLattice
    reference: object
    atom_masses_amu: np.ndarray
    source_manifest: dict
    producer_alpha: float
    computed_alpha: float


def load_material_frozen_lattice(source_manifest):
    """Consume A22's saved arrays and admitted XML/UPF/density, never WFCs.

    source_manifest is the explicit A28 ionic-q source registration path.
    All material values come from authenticated original/admitted bytes.
    This reader does not reimport or project A22's matrix.
    """
    with open(source_manifest, encoding='utf-8') as stream:
        manifest = json.load(stream)
    sources = manifest['inputs']
    receipt_bytes, _ = read_bound_bytes(sources['a22_receipt'], max_bytes=32 << 20)
    receipt = json.loads(receipt_bytes)
    if receipt['status'] != 'complete' or receipt['archive'] != sources['a22_archive']:
        raise ValueError("saved A22 output binding differs")
    component = receipt['material_component']
    identity = component['source_identity']
    if component['mass_weighted'] or component['induced_density_screening_included'] or component['units'] != 'Ry/bohr2':
        raise ValueError("unmassweighted unscreened imported curvature required")
    ref = load_saved_reference(sources['xml'])
    if (tuple(identity['atom_symbols']) != ref.atom_names
            or identity['ground_state_xml']['sha256'] != ref.source['sha256']
            or not np.array_equal(np.asarray(identity['cell']['vectors_rows']).T, ref.cell_bohr)
            or not np.array_equal(np.asarray(identity['cell']['positions_atom_rows']).T, ref.positions_bohr)):
        raise ValueError("A22 and original XML geometry/reference differ")
    valences = {}
    for species in ref.species:
        binding = sources['upfs'][species.name]
        pp = read_nc_upf(binding['path'], expected_bytes=binding['bytes'], expected_sha256=binding['sha256'])
        if pp.element != species.name:
            raise ValueError("NC pseudovalence species mismatch")
        valences[species.name] = pp.z_valence
    charges = np.asarray([valences[s] for s in ref.atom_names])
    density = load_density(sources['density'])
    if density.gamma_only or density.ngm != ref.ngm or not np.allclose(density.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, atol=2e-14, rtol=2e-14):
        raise ValueError("original full signed density G support differs")
    stdout, _ = read_bound_bytes(sources['ph_stdout'])
    matches = re.findall(rb'Alpha used in Ewald sum\s*=\s*([0-9.]+)', stdout)
    if len(matches) != 1:
        raise ValueError("one original PH alpha report required")
    printed_alpha = float(matches[0])
    alpha = qe_alpha(total_valence=charges.sum(), density_cutoff_ry=2 * ref.ecutrho_ha)
    if f'{alpha:.4f}' != f'{printed_alpha:.4f}':
        raise ValueError("source cutoff/valences do not reproduce original PH alpha")
    raw, _ = read_bound_bytes(sources['a22_archive'])
    with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
        gamma = archive['cartesian_matrix']
        masses = archive['PH_atom_mass_amu']
    if not np.array_equal(masses, np.asarray(identity['PH_atom_mass_amu'])):
        raise ValueError("saved PH masses differ from imported provenance")
    ionic = QEIonicCurvature(cell_bohr=ref.cell_bohr, positions_bohr=ref.positions_bohr,
                            valence_charges=charges, g_cart_bohr_inv=density.reciprocal_bohr_inv @ density.miller,
                            alat_bohr=ref.alat_bohr, alpha_bohr_inv2=alpha)
    return MaterialFrozenLattice(FrozenLattice(gamma, ionic), ref, masses, manifest, printed_alpha, alpha)
