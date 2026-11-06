"""Exact stationary-primary intake, without loaded-only production authorities.

Callers supply bound case manifest, native cache receipt, and primary admission.
No case name is mapped to an implicit path. In particular gamma-continuation-002
is NOT substituted for loaded-primary. Original XML geometry, species, cutoff,
spin, raw occupations and native quadrature remain separately identifiable.
"""
from dataclasses import dataclass
import json
from pathlib import Path
from types import MappingProxyType

import numpy as np

from .qe_saved_native import read_bound_bytes, load_wavefunctions
from .qe_saved_reference import load_saved_reference
from .qe_nc_upf import read_nc_upf
from .qe_bloch_frame import admit_source, _digest
from .qe_real_bloch_grid import (
    REAL_POLICY, admit_real_hamiltonian, build_real_bloch_grid, real_basis,
)
from .material_preparation import (
    SourceBinding, FullGridSector, MaterialPreparationPolicy, specify_material_preparation,
)


def clean_binding(value):
    return {field: str(value[field]) for field in ('path', 'bytes', 'sha256')}


def bound_json(value):
    raw, _ = read_bound_bytes(value, max_bytes=32 << 20)
    return json.loads(raw)


def _same_binding(left, right):
    return clean_binding(left) == clean_binding(right)


@dataclass(frozen=True)
class MaterialControlInputs:
    case: str
    dataset: str
    manifest_binding: object
    native_receipt_binding: object
    admission_binding: object
    final_attempt: str
    geometry_binding: object
    attempt_geometry_binding: object
    xml_binding: object
    density_binding: object
    upf_bindings: object
    wfc_bindings: object
    reference: object
    upfs: object
    scope: object
    grid: object
    model_settings: object

    def load_material_inputs(self, *, qmax_bohr_inv=None):
        """Original generic H0/K input API; no loaded wrappers or new SCF."""
        from .qe_material_reference import load_material_inputs
        return load_material_inputs(self.xml_binding, self.density_binding, self.upf_bindings,
                                    qmax_bohr_inv=qmax_bohr_inv)

    def provenance(self):
        return dict(case=self.case, dataset=self.dataset, final_attempt=self.final_attempt,
            case_manifest=dict(self.manifest_binding), native_receipt=dict(self.native_receipt_binding),
            primary_admission=dict(self.admission_binding), geometry=dict(self.geometry_binding),
            attempt_geometry=dict(self.attempt_geometry_binding), xml=dict(self.xml_binding),
            density=dict(self.density_binding), upfs={k: dict(v) for k, v in self.upf_bindings.items()},
            wavefunctions=[dict(ispin=s, source_ik=k, binding=dict(b))
                           for (s, k), b in sorted(self.wfc_bindings.items())],
            model_settings=dict(self.model_settings), real_scope_identity=self.scope.identity,
            symmetry_policy=REAL_POLICY,
            geometry_authority='original selected final XML atom order/cell/positions; separately bound accepted final geometry',
            continuation_equivalence='not asserted; no Gamma-continuation bytes substituted')


def load_control_inputs(case_manifest_binding, native_receipt_binding, admission_binding):
    """Authenticate accepted native bindings; WFC byte decoding stays lazy.

    The admission authenticates the supplied manifest and receipt; its selected
    final attempt chooses files from completed_cache_before, never glob order.
    Binary readers verify the very buffers decoded on subsequent consumption.
    """
    manifest = bound_json(case_manifest_binding)
    receipt = bound_json(native_receipt_binding)
    admission = bound_json(admission_binding)
    case = manifest['id']
    entry = admission['materials'][case]
    if (entry['status'] != 'admit'
            or not _same_binding(entry['native_receipt'], native_receipt_binding)
            or not _same_binding(entry['baseline_spec'], case_manifest_binding)
            or receipt['status'] != 'actual-completed-material-force-native-source-science-and-cache-verified'
            or receipt['whole_completed_cache_exactly_preserved'] is not True
            or receipt['whole_publication_exactly_preserved'] is not True):
        raise ValueError('explicit accepted primary manifest/native receipt bindings required')
    published = [r for r in receipt['original_frozen_publication']
                 if r['path'] == case_manifest_binding['path']]
    if len(published) != 1 or not _same_binding(published[0], case_manifest_binding):
        raise ValueError('case manifest is not the original frozen publication')
    attempt = Path(entry['final_scf']['attempt_directory'])
    case_dir = Path(manifest['directory'])
    if (str(attempt) != receipt['completed_attempt'] or case_dir not in attempt.parents
            or case_dir.parent.name != 'primary'
            or entry['final_scf']['returncode'] != 0 or not entry['final_scf']['job_done']):
        raise ValueError('different selected stationary-primary attempt; continuation is not interchangeable')
    rows = receipt['completed_cache_before']
    cache = {row['path']: clean_binding(row) for row in rows}
    if len(cache) != len(rows):
        raise ValueError('duplicate native cache path')
    save = attempt / 'scratch/scf.save'

    def original(name):
        path = str(save / name)
        if path not in cache:
            raise ValueError('required original saved file absent from native bindings: ' + name)
        return MappingProxyType(cache[path])

    xml_binding = original('data-file-schema.xml')
    expected_xml = entry['final_scf']['artifacts']['data-file-schema.xml']
    if (xml_binding['sha256'] != expected_xml['sha256']
            or int(xml_binding['bytes']) != int(expected_xml['bytes'])):
        raise ValueError('selected saved XML differs from admitted final attempt XML')
    reference = load_saved_reference(xml_binding)
    settings = manifest['scf_settings']
    mesh = [reference.monkhorst_pack['nk' + str(i)] for i in (1, 2, 3)]
    if (2 * reference.ecutwfc_ha != settings['ecutwfc_ry']
            or 2 * reference.ecutrho_ha != settings['ecutrho_ry']
            or reference.degauss_ry != settings['degauss_ry']
            or reference.smearing != settings['smearing'] or settings['nspin'] != 2
            or len(reference.nbnd_by_spin) != 2
            or any(n != settings['nbnd'] for n in reference.nbnd_by_spin)
            or mesh != settings['kpts'] or any(reference.monkhorst_pack['k' + str(i)] for i in (1, 2, 3))):
        raise ValueError('saved cutoff/spin/band/k/smearing model differs from exact case manifest')
    geometry = clean_binding(entry['final_artifacts']['final.extxyz'])
    attempt_geometry = clean_binding(entry['final_scf']['artifacts']['final.extxyz'])
    # Preserve BOTH accepted geometry artifacts; never use the registered seed.
    for b in (geometry, attempt_geometry):
        read_bound_bytes(b, max_bytes=1 << 20)
    if not _same_binding(cache[str(attempt / 'final.extxyz')], attempt_geometry):
        raise ValueError('attempt final geometry is not native-cache bound')
    upf_bindings, upfs = {}, {}
    for species in reference.species:
        if manifest['execution']['pseudopotentials'][species.name] != species.pseudo_file:
            raise ValueError('manifest/XML UPF name mismatch')
        b = original(species.pseudo_file)
        original_path = str(Path(manifest['execution']['pseudo_dir']) / species.pseudo_file)
        if manifest['input_files_sha256'][original_path] != b['sha256']:
            raise ValueError('copied UPF differs from registered model')
        upf_bindings[species.name] = b
        upfs[species.name] = read_nc_upf(b['path'], expected_bytes=b['bytes'], expected_sha256=b['sha256'])
    wfc = {(spin, point.ik): original(f'wfc{label}{point.ik}.dat')
           for spin, label in ((1, 'up'), (2, 'dw')) for point in reference.kpoints}
    saved_wfc = {Path(p).name for p in cache if Path(p).parent == save and Path(p).name.startswith('wfc')}
    if saved_wfc != {Path(b['path']).name for b in wfc.values()}:
        raise ValueError('native WFC set differs from complete explicit original spin/source set')
    scope = admit_real_hamiltonian(reference, upfs)
    grid = build_real_bloch_grid(reference, wfc, scope=scope, symmetry_policy=REAL_POLICY)
    return MaterialControlInputs(case, 'stationary-primary',
        MappingProxyType(clean_binding(case_manifest_binding)), MappingProxyType(clean_binding(native_receipt_binding)),
        MappingProxyType(clean_binding(admission_binding)), str(attempt), MappingProxyType(geometry),
        MappingProxyType(attempt_geometry), xml_binding, original('charge-density.dat'),
        MappingProxyType(upf_bindings), MappingProxyType(wfc), reference, MappingProxyType(upfs),
        scope, grid, MappingProxyType(dict(settings)))


def iter_control_source_admissions(control):
    """One original WFC at a time; every row/header/cutoff and target basis.

    No source Gram is computed here. Consumers must not call this full-grid
    state/observable execution. All original raw columns are retained.
    """
    grid = control.grid
    for key, binding in sorted(control.wfc_bindings.items()):
        source = load_wavefunctions(binding)
        slots = tuple(s for s in grid.slots if (s.ispin, s.source_ik) == key)
        audit = admit_source(control.reference, slots[0], source)
        coefficient_digest = _digest(source.coefficients)
        descriptors = []
        for slot in slots:
            basis = real_basis(grid, slot, source, audit)
            descriptors.append(dict(spin=slot.ispin, mesh_index=list(slot.mesh_index),
                source_ik=slot.source_ik, source_binding=dict(slot.source_binding),
                basis_identity=basis.identity, frame_identity=basis.identity + ':raw-saved-columns',
                rows=source.igwx, columns=source.nbnd, operation=slot.operation.index,
                orbital_conjugation=slot.orbital_conjugation,
                reciprocal_shift=slot.reciprocal_shift.tolist(), k_fractional=slot.k_fractional.tolist(),
                source_k_fractional=slot.source_k_fractional.tolist(), measure=slot.weight,
                original_xml_weight=control.reference.kpoints[source.ik - 1].weight,
                source_rows_sha256=audit['source_rows_sha256'], target_rows_sha256=_digest(basis.miller),
                source_coefficients_sha256=coefficient_digest))
        yield source, audit, tuple(descriptors)
        del source


@dataclass(frozen=True)
class RealMaterialPreparationPolicy(MaterialPreparationPolicy):
    """Explicit extension of the caller-supplied policy, NOT unitary authority.

    The neutral filler is symmetry-agnostic: FullGridSector metadata supplies
    the extension. This subclass overrides the old policy docstring's unitary
    scope with A28's bound SAME-spin orbital conjugation and raw Gram/D rule.
    """


def control_population_specification(control, descriptors, *, coverage_authority, decision_authority):
    """Existing source-backed neutral filling; MV never used as probabilities."""
    by_key = {(int(d['spin']), tuple(d['mesh_index'])): d for d in descriptors}
    if len(by_key) != len(descriptors) or set(by_key) != {(s.ispin, s.mesh_index) for s in control.grid.slots}:
        raise ValueError('complete independently admitted slot descriptor set required')
    sectors = []
    for slot in control.grid.slots:
        d = by_key[slot.ispin, slot.mesh_index]
        if (not _same_binding(d['source_binding'], slot.source_binding)
                or d['source_ik'] != slot.source_ik or d['rows'] != slot.expected_rows
                or d['columns'] != control.reference.nbnd_by_spin[slot.ispin - 1]
                or d['orbital_conjugation'] != slot.orbital_conjugation
                or d['measure'] != slot.weight or d['operation'] != slot.operation.index
                or d['reciprocal_shift'] != slot.reciprocal_shift.tolist()):
            raise ValueError('slot descriptors do not match original source/extension')
        sectors.append(FullGridSector(slot.ispin, slot.mesh_index, slot.source_ik,
            SourceBinding.from_mapping(slot.source_binding), d['frame_identity'], d['basis_identity'],
            slot.expected_rows, REAL_POLICY + ':' + control.scope.identity))
    policy = RealMaterialPreparationPolicy(decision_authority,
        'actual atom-resolved bound NC UPF valence sum equals saved nelec',
        'Q=Psi(Gram)^(-1/2); original saved Ry labels diagonal in Q; conjugate Gram AND D on orbital conjugation',
        'exact equal saved energy labels only; symmetric fractional filling at exact Fermi degeneracy',
        REAL_POLICY + '; same explicit spin; no physical spin-flipping time reversal; scope=' + control.scope.identity,
        'independent quasifree explicit spin/full-k CAR sectors; free magnetization/common multiplier',
        'exactly empty initial full-PW orthogonal complement; later actions retain all PW rows',
        coverage_authority)
    return specify_material_preparation(control.reference, control.upfs, sectors, policy=policy)


def control_population_arrays(control, specification):
    """Case-independent complete ledger (ragged bands use explicit offsets)."""
    slots = control.grid.slots
    if tuple(s.key for s in specification.sectors) != tuple((s.ispin, s.mesh_index) for s in slots):
        raise ValueError('population ordering differs from complete admitted grid')
    return dict(populations=np.asarray(specification.filling.populations),
        energies_ry=np.asarray([level.energy for level in specification.filling.levels]),
        band_slices=np.asarray(specification.band_slices, dtype=np.int64),
        spin=np.asarray([s.ispin for s in slots]), mesh_index=np.asarray([s.mesh_index for s in slots]),
        source_ik=np.asarray([s.source_ik for s in slots]),
        source_sha256=np.asarray([s.source_binding['sha256'] for s in slots]),
        measures=np.asarray([s.weight for s in slots]),
        orbital_conjugation=np.asarray([s.orbital_conjugation for s in slots]),
        reciprocal_shift=np.asarray([s.reciprocal_shift for s in slots]),
        basis_identity=np.asarray([s.pw_basis_identity for s in specification.sectors]),
        frame_identity=np.asarray([s.frame_identity for s in specification.sectors]))
