"""Source-bound saved-label population policy, separate from raw-D instantiation.

This conditional research family uses Q=Psi(Psi†Psi)^(-1/2), diagonal saved Ry
energy labels in Q, free M/common multiplier, independent spin/full-k sectors,
and an exactly empty initial full-PW complement. Future actions retain every
PW row. A canonical unitary extension is NOT little-group invariance. Neither
saved MV occupations nor degauss enters. No H0 stationarity is asserted.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import fsum, isfinite, prod
from typing import Mapping

import numpy as np

from .physical_factor_preparation import (
    FactorFunctional, FactorPreparation, coefficients_from_orthonormal_populations,
)
from .qe_saved_reference import SavedReference
from .zero_temperature_preparation import (
    FillingPolicy, ModeLevel, ZeroTemperatureFilling, fill_zero_temperature,
)


@dataclass(frozen=True)
class SourceBinding:
    path: str
    bytes: int
    sha256: str

    def __post_init__(self):
        if (not self.path.strip() or type(self.bytes) is not int or self.bytes <= 0
                or len(self.sha256) != 64
                or any(c not in '0123456789abcdef' for c in self.sha256)):
            raise ValueError('original path/count/SHA256 binding required')

    @classmethod
    def from_mapping(cls, source):
        return cls(str(source['path']), int(source['bytes']), str(source['sha256']))


@dataclass(frozen=True)
class FullGridSector:
    """One distinct CAR sector, not a weighted copy of a single Gamma mode.

    The caller supplies the admitted canonical extension and frame identities;
    this module does not import or implement Bloch unfolding. Band order is the
    unchanged source order; arbitrary additional band gauges need transported
    operators/populations and are not accepted by this diagonal-label adapter.
    """
    spin: int
    mesh_index: tuple[int, int, int]
    source_ik: int
    source_binding: SourceBinding
    frame_identity: str
    pw_basis_identity: str
    full_pw_dimension: int
    extension_identity: str

    @property
    def key(self):
        return (self.spin, self.mesh_index)


@dataclass(frozen=True)
class MaterialPreparationPolicy:
    """All authorities are caller-supplied original source/decision identities.

    extension_authority adopts XML-order unitary canonical transport (including
    saved inversion), no antiunitary/spin swap, identity on representatives and
    otherwise first valid saved operation. It does not certify path invariance.
    coverage_authority identifies complete source/header/row/slot admission, not
    an assertion that a selected-sector output is the neutral material total.
    """
    decision_authority: str
    neutral_number_authority: str
    symmetric_frame_authority: str
    exact_degeneracy_authority: str
    extension_authority: str
    independence_authority: str
    empty_complement_authority: str
    coverage_authority: str

    def __post_init__(self):
        if any(not isinstance(v, str) or not v.strip() for v in vars(self).values()):
            raise ValueError('all conditional material policy authorities required')


@dataclass(frozen=True)
class MaterialPopulationSpecification:
    """Population metadata only: no Gram, raw D, material rho or equilibrium."""
    reference_binding: SourceBinding
    charge_bindings: tuple[tuple[str, SourceBinding], ...]
    valence_charges: tuple[tuple[str, float], ...]
    mesh: tuple[int, int, int]
    sectors: tuple[FullGridSector, ...]
    band_slices: tuple[tuple[int, int], ...]
    filling: ZeroTemperatureFilling
    policy: MaterialPreparationPolicy


def _neutral_charge(reference, upfs):
    charges, bindings = {}, []
    for species in reference.species:
        upf = upfs[species.name]
        if upf.element != species.name or not isfinite(upf.z_valence) or upf.z_valence <= 0:
            raise ValueError('actual NC UPF valence/species identity required')
        charges[species.name] = float(upf.z_valence)
        bindings.append((species.name, SourceBinding(upf.path, len(upf.original_bytes), upf.sha256)))
    neutral = fsum(charges[name] for name in reference.atom_names)
    # Uses decoded actual valence and saved nelec, not nominal nuclear charges
    # or a Gamma trace. Their difference is the saved net-charge consistency.
    if reference.nelec != neutral:
        raise ValueError('saved nelec must equal atom-resolved neutral UPF valence sum')
    return neutral, tuple(charges.items()), tuple(bindings)


def _complete_sectors(reference, sectors):
    if not reference.lsda or reference.noncolin or reference.spinorbit or len(reference.nbnd_by_spin) != 2:
        raise ValueError('two explicit collinear spin sectors required; no implicit spin two')
    mesh = tuple(reference.monkhorst_pack.get('nk' + str(i), 0) for i in (1, 2, 3))
    if min(mesh) < 1 or any(reference.monkhorst_pack.get('k' + str(i), -1) != 0 for i in (1, 2, 3)):
        raise ValueError('explicit complete unshifted Monkhorst-Pack mesh required')
    expected = {(spin, index) for spin in (1, 2) for index in product(*(range(n) for n in mesh))}
    if len(sectors) != len(expected) or {s.key for s in sectors} != expected:
        raise ValueError('every full-grid slot exactly once per explicit spin required')
    for sector in sectors:
        reference.select(sector.source_ik, sector.spin)
        if (not isinstance(sector.source_binding, SourceBinding)
                or type(sector.full_pw_dimension) is not int
                or sector.full_pw_dimension < reference.nbnd_by_spin[sector.spin - 1]
                or any(not value.strip() for value in (sector.frame_identity,
                       sector.pw_basis_identity, sector.extension_identity))):
            raise ValueError('complete original source/frame/full-row extension identity required')
    for spin in (1, 2):
        if {s.source_ik for s in sectors if s.spin == spin} != set(range(1, reference.nks + 1)):
            raise ValueError('all original spin-local source records must be represented')
    return mesh


def specify_material_preparation(reference: SavedReference, upfs: Mapping, sectors, *,
                                 policy: MaterialPreparationPolicy) -> MaterialPopulationSpecification:
    """Fill neutral mean N using complete full-grid metadata, without a raw D.

    upfs maps existing XML species names to decoded bound NCUPF instances.
    nelec is checked against the actual atom-resolved UPF valence sum. One full
    k mode has measure 1/Nk per primitive cell, once per explicit spin; XML IBZ
    weights are never used as probabilities or as fluctuation weights. Caller
    supplies complete canonical slot metadata; no actual decoding happens here.
    """
    if not isinstance(policy, MaterialPreparationPolicy):
        raise TypeError('explicit conditional material policy required')
    sectors = tuple(sectors)
    mesh = _complete_sectors(reference, sectors)
    neutral, charges, charge_bindings = _neutral_charge(reference, upfs)
    levels, slices = [], []
    for sector in sectors:
        _, energies, _ = reference.select(sector.source_ik, sector.spin)
        start = len(levels)
        if len(energies) != reference.nbnd_by_spin[sector.spin - 1]:
            raise ValueError('source band count differs from saved spin metadata')
        for band, energy in enumerate(energies):
            levels.append(ModeLevel(f'{sector.key}:band={band + 1}', sector.spin,
                                    float(energy), 1.0 / prod(mesh)))
        slices.append((start, len(levels)))
    filling_policy = FillingPolicy(
        policy.symmetric_frame_authority, policy.neutral_number_authority,
        'explicit independent full-k slots; 1/Nk per primitive cell; ' + policy.coverage_authority,
        policy.exact_degeneracy_authority, 'Ry', 'electrons per primitive cell')
    filling = fill_zero_temperature(levels, neutral, policy=filling_policy)
    return MaterialPopulationSpecification(SourceBinding.from_mapping(reference.source),
        charge_bindings, charges, mesh, sectors, tuple(slices), filling, policy)


def instantiate_material_sector(specification: MaterialPopulationSpecification, sector_key,
                                functional: FactorFunctional, raw_frame, *,
                                source_binding: SourceBinding, reference_binding: SourceBinding,
                                coverage_authority: str, frame_authority: str) -> FactorPreparation:
    """Instantiate ONE admitted sector lazily; never label it the neutral total.

    Raw frame is the unchanged all-row/all-column canonical saved-label frame.
    Its recomputed Gram must match the functional. Byte verification and Bloch
    row/header admission are producer responsibilities bound by these required
    identities. No PW-square density, truncation or band projection is formed.
    """
    if (reference_binding != specification.reference_binding
            or coverage_authority != specification.policy.coverage_authority
            or not frame_authority.strip()):
        raise ValueError('matching original XML, complete coverage and raw frame authority required')
    try:
        index = next(i for i, sector in enumerate(specification.sectors) if sector.key == sector_key)
    except StopIteration:
        raise ValueError('sector not in complete population specification') from None
    sector = specification.sectors[index]
    start, stop = specification.band_slices[index]
    psi = np.asarray(raw_frame)
    if (source_binding != sector.source_binding or functional.factor_identity != sector.frame_identity
            or functional.pw_basis_identity != sector.pw_basis_identity
            or functional.full_pw_dimension != sector.full_pw_dimension
            or psi.shape != (sector.full_pw_dimension, stop - start)
            or not np.isfinite(psi).all()):
        raise ValueError('source/frame/full-PW rows/band identities differ from specification')
    gram = psi.conj().T @ psi
    allowance = 128 * np.finfo(float).eps * max(psi.shape) * (1 + np.max(np.abs(gram)))
    if functional.gram.shape != gram.shape or np.max(np.abs(functional.gram - gram)) > allowance:
        raise ValueError('functional Gram is not that of the admitted raw frame')
    mapping = specification.policy.symmetric_frame_authority + '; ' + frame_authority
    coefficients = coefficients_from_orthonormal_populations(
        functional.gram, specification.filling.populations[start:stop], mapping_authority=mapping)
    authority = (specification.policy.decision_authority + '; sector=' + str(sector.key)
                 + '; independent sectors: ' + specification.policy.independence_authority
                 + '; initially empty complement: ' + specification.policy.empty_complement_authority
                 + '; complete coverage: ' + coverage_authority + '; ' + mapping
                 + '; conditional mean-N preparation, not exact-N or H0 equilibrium; sector is not total')
    return FactorPreparation(functional, coefficients, authority)
