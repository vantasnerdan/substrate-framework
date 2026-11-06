"""One finite periodic ion-ion scalar and its ionic LINEAR expansion term.

QE7.5 PW/src/ewald.f90:94-174, force_ew.f90:98-179 and Modules/rgen.f90:47-79
supply the ordinary 3D valence-charge conventions, self/background constants
and real-shell refolding. This is NOT either Native finite-cutoff path:
ewald chooses alpha from 2.8 down with 1e-7 energy bound and rmax=4/sqrt(alpha);
force_ew starts at 1.0 with 1e-6 bound and rmax=5/sqrt(alpha). Here alpha and
both physical cutoffs are explicit, fixed numerical inputs to ONE scalar.

U = e2/2 sum_ijR' Zi Zj erfc(sqrt(alpha)|Ri-Rj-R|)/|Ri-Rj-R|
  + e2*2pi/Omega sum_G!=0 exp(-G^2/(4alpha)) |sum_i Zi exp(iG.Ri)|^2/G^2
  - e2*sqrt(alpha/pi) sum_i Zi^2 - e2*pi*(sum_i Zi)^2/(2alpha*Omega).
e2=2 (Ry); uniform neutralizing background, zero G=0 potential/tinfoil 3D
boundary, no surface dipole, 2D cutoff, ESM or MT. The finite hard real cutoff
makes U piecewise smooth; F=-dU/dR is returned only off its crossing surfaces.
No convergence/Native equality follows merely from evaluating this contract.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from math import pi, sqrt
from types import MappingProxyType

import numpy as np

from .qe_ewald_sums import real_sum, reciprocal_sum


def _frozen_real(value, shape, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) or raw.shape != shape or not np.isfinite(raw).all():
        raise ValueError(name + ' must be a finite real array of shape ' + str(shape))
    result = np.array(raw, dtype=float, copy=True)
    result.flags.writeable = False
    return result


def _authority(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('explicit ' + name + ' required')


@dataclass(frozen=True)
class EwaldSumContract:
    """Numerical sum definition, not a physical screening parameter.

    alpha_bohr_inv2 is QE alpha (NOT its square root). Both spherical sums
    include their cutoff boundary; differentiation holds cell/alpha/cutoffs
    fixed. Real cutoff crossings are rejected at floating-point resolution.
    """
    alpha_bohr_inv2: float
    real_cutoff_bohr: float
    reciprocal_cutoff_bohr_inv: float

    def __post_init__(self):
        for name in ('alpha_bohr_inv2', 'real_cutoff_bohr', 'reciprocal_cutoff_bohr_inv'):
            value = getattr(self, name)
            if isinstance(value, (bool, complex)) or not np.isscalar(value) or not np.isfinite(value) or value <= 0:
                raise ValueError('positive finite numerical ' + name + ' required')


@dataclass(frozen=True)
class IonicLinearTerm:
    """Source-bound result at Rbar; use ONLY ionic_force in PreparedReference.

    Coordinates and forces are [nat,3] Cartesian atom-major, cell columns in
    bohr. ionic_force/gradient flatten as atom0-x,y,z, atom1-x,y,z,... .
    energy_ry is a reference diagnostic, not an extra displaced potential.
    """
    energy_ry: float
    force_ry_bohr: np.ndarray
    real_force_ry_bohr: np.ndarray
    reciprocal_force_ry_bohr: np.ndarray
    cell_bohr: np.ndarray
    positions_bohr: np.ndarray
    valence_charges: np.ndarray
    coordinate_identity: str
    source_authority: str
    sources: tuple
    contract: EwaldSumContract
    energy_terms_ry: object
    sum_metadata: object

    @property
    def ionic_force(self):
        """Read-only full atom-major Cartesian F_II [3*nat], Ry/bohr."""
        return self.force_ry_bohr.reshape(-1)

    @property
    def gradient_ry_bohr(self):
        """ell_II=dU/dR, atom-major [3*nat]; no fitted offset."""
        result = -self.ionic_force
        result.flags.writeable = False
        return result

    @property
    def ionic_authority(self):
        """Serializable binding for PreparedReference.ionic_authority."""
        return json.dumps(dict(
            scalar='qe_ionic_ewald ordinary 3D finite scalar; F_II=-dU/dR',
            coordinate_identity=self.coordinate_identity,
            positions_bohr=self.positions_bohr.tolist(), cell_bohr=self.cell_bohr.tolist(),
            valence_charges=self.valence_charges.tolist(), source_authority=self.source_authority,
            sources=[dict(source) for source in self.sources],
            contract=vars(self.contract), metadata=dict(self.sum_metadata)), sort_keys=True)



def periodic_ionic_ewald(cell_bohr, positions_bohr, valence_charges, *, contract,
                         coordinate_identity, source_authority, sources):
    """Evaluate a finite scalar/gradient from explicit physical inputs.

    positions_bohr is [nat,3], unlike SavedReference's [3,nat]. sources is an
    explicit tuple of original input byte bindings; synthetic callers supply
    () and identify their defining literals in source_authority. No decoding,
    alpha optimization, force centering or material-state choice occurs here.
    """
    if not isinstance(contract, EwaldSumContract):
        raise TypeError('explicit EwaldSumContract required')
    _authority(coordinate_identity, 'Cartesian atom-order identity')
    _authority(source_authority, 'geometry/UPF valence authority')
    charges = np.asarray(valence_charges)
    if charges.ndim != 1 or len(charges) < 1:
        raise ValueError('one valence charge per ion required')
    charges = _frozen_real(charges, charges.shape, 'UPF valence charges')
    if np.any(charges <= 0):
        raise ValueError('positive UPF valence charges required, not nuclear numbers')
    cell = _frozen_real(cell_bohr, (3, 3), 'cell columns in bohr')
    positions = _frozen_real(positions_bohr, (len(charges), 3), 'Cartesian positions in bohr')
    volume = abs(float(np.linalg.det(cell)))
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError('nonsingular finite 3D periodic cell required')
    alpha = contract.alpha_bohr_inv2
    real, fr, nr, gap = real_sum(cell, positions, charges, alpha, contract.real_cutoff_bohr)
    reciprocal, fg, ng = reciprocal_sum(cell, positions, charges, alpha, contract.reciprocal_cutoff_bohr_inv)
    terms = dict(real=real, reciprocal=reciprocal,
                 self=-2 * sqrt(alpha / pi) * float(charges @ charges),
                 background=-pi * float(charges.sum())**2 / (alpha * volume))
    metadata = dict(boundary='ordinary 3D periodic, neutralizing background, tinfoil',
                    energy_units='Ry', force_units='Ry/bohr', e2_ry_bohr=2.,
                    cell_vectors='columns', coordinate_order='atom-major Cartesian xyz',
                    real_ordered_pair_terms=nr, reciprocal_signed_vectors=ng,
                    real_cutoff_gap_bohr=gap, derivative='fixed cell/alpha/cutoffs; no shell crossing',
                    zero_exclusion='exact i=j,R=0 only; coincident distinct ions rejected',
                    qe_rgen_difference='no 1e-10 alat-squared near-collision omission or mxr=50 cap',
                    native_match='not asserted; Native energy and force have distinct finite sums',
                    convergence='unmeasured; explicit cutoffs are not an error certificate',
                    processing='none: no net-force subtraction, symmetrization or fitted support')
    return IonicLinearTerm(
        float(sum(terms.values())), _frozen_real(fr + fg, positions.shape, 'ionic force'),
        _frozen_real(fr, positions.shape, 'real force'), _frozen_real(fg, positions.shape, 'reciprocal force'),
        cell, positions, charges, coordinate_identity, source_authority,
        tuple(MappingProxyType(dict(source)) for source in sources), contract,
        MappingProxyType(terms), MappingProxyType(metadata))


def ionic_linear_term(material_inputs, *, contract, coordinate_identity):
    """Consume admitted MaterialInputs, using actual NCUPF.z_valence only.

    Reuses decoded SavedReference/UPFs; does not read WFCs, density, XML forces,
    nuclear charges or saved ewald labels. Original atom labels and source
    bindings are carried. Caller owns the admitted input-object provenance.
    """
    reference = material_inputs.reference
    charges, sources = [], [dict(reference.source)]
    for species in reference.species:
        upf = material_inputs.tables[species.name].upf
        if upf.element != species.name:
            raise ValueError('UPF/XML species mismatch')
        sources.append(dict(path=upf.path, bytes=str(len(upf.original_bytes)), sha256=upf.sha256))
    for name in reference.atom_names:
        charges.append(material_inputs.tables[name].upf.z_valence)
    if len(charges) != reference.nat or len(reference.atom_indices) != reference.nat:
        raise ValueError('saved atom order/count mismatch')
    atom_identity = tuple(zip(reference.atom_indices, reference.atom_names))
    return periodic_ionic_ewald(
        reference.cell_bohr, reference.positions_bohr.T, charges, contract=contract,
        coordinate_identity=coordinate_identity,
        source_authority='SavedReference geometry and NCUPF VALENCE charges; atom labels=' + repr(atom_identity),
        sources=tuple(sources))
