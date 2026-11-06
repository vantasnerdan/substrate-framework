"""Restore original material geometry/full state and compute ONLY missing third jets.

A29's persisted original basis is reused; active12 overlaps cannot represent
full-column state D. The necessary full-column moment overlaps are contracted
in blocks against the established physical placed-projector factors. No H/F/J,
SCF, Gaussian preparation, lower curvature or native WFC replay is invoked.
"""
import json
from pathlib import Path
import sys
from time import perf_counter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'),str(HERE.parent/'0041'),str(HERE.parent/'0029'),str(HERE.parent/'0028')]
import matched_reference_io as io
from matched_joint_sources import registration, validate_aggregate, source_catalog
from reduced_prepared_projection import restore_basis
from reduced_projection_storage import ImmutableItems
from reduced_projection_cache import ProjectionCache, atomic_json
from substrate_framework.material_control_reference import bound_json, bound_arrays
from substrate_framework.material_physical_operators import physical_coupling
from substrate_framework.streamed_material_state import DiskMean, authenticate, file_binding
from substrate_framework.partitioned_q_storage import release_mapped_pages
from cubic_electronic_donor import full_state_third
from cubic_position_moments import position_factors

BOUNDARY = dict(path='/var/tmp/p254-a28-matched-loaded-rebind-001/original/boundary.json',
                bytes=277032,sha256='6fb559b067146ce2a44e5ecdf29a0f6045477dec4af769714830049af55752e2')
MATCHING = dict(
    definition='E_ext=E_original+(E_donor-Taylor_0_1_2 E_donor); only third jet retained here',
    lower_order='EXACT original force and complete C(q), including old d2ionq; no donor Hessian substitution',
    ionic='ordinary 3D finite Ewald, alpha1.0 bohr^-2, real12 bohr, reciprocal12 bohr^-1; spherical physical BvK K support',
    electronic='fixed original full128 explicit spin/k state; NC local plus all8 nonlocal position-phase product legs; all raw columns/complex D',
    normalization='u_Ra=N^-1/2 sum_q exp(+2pi i q.R)u_qa, N64; Phi carries 1/sqrt64; no masses or atomic q phase',
    omitted='density-dependent quadratic electron-ion operator, selfconsistent/passive dynamics, rates, screening and field-cubic terms')


def restore_original(manifest):
    _,data = io.registration()
    control,_,ledger,rows,context = io.load_case('loaded-ni31h6',data,background=False)
    _,joint,_ = registration()
    aggregate_binding = joint['cases']['loaded-ni31h6']
    aggregate = bound_json(aggregate_binding)
    rows = tuple(sorted(rows,key=lambda row:(row[1].ispin,row[1].mesh_index)))
    validate_aggregate(aggregate,context,rows)
    states,_ = source_catalog(aggregate,context,rows,ledger)
    original = DiskMean(authenticate(BOUNDARY),verify=False)
    spec = manifest['specification']
    if (spec['original_boundary'] != BOUNDARY or original.time != 0
            or np.any(original.u) or np.any(original.p) or len(original.slots) != 128):
        raise ValueError('exact full original time0 boundary required')
    inputs = control.load_material_inputs()
    return control,inputs,rows,states,original,aggregate_binding,context


def _sector_lineage(index,slot,record,item,source_state):
    metadata = item['metadata']
    if (metadata['frame'] != record.frame or metadata['raw_d'] != record.coefficient
            or metadata['basis_identity'] != record.basis_identity
            or metadata['frame_identity'] != record.gauge_identity
            or metadata['source']['sha256'] != record.source_identity
            or record.weight != slot.weight or record.component != slot.ispin-1
            or metadata['source_state'] != source_state):
        raise ValueError('original cached basis/raw state/source/explicit measure differs')
    return dict(index=index,identity=record.identity,weight=record.weight,
        component=record.component,frame=record.frame,coefficient=record.coefficient,
        preparation_archive=item['archive'],source_state=source_state,
        source_archive=metadata['source_archive'],basis_identity=record.basis_identity,
        source=dict(slot.source_binding))


def _new_sector(index,slot,record,lineage,reader,control,inputs,original,electronic_atoms=None):
    geometry_names = ('miller','k_cart_bohr_inv','reciprocal_bohr_inv','row_map','source_miller')
    key = 'preparation/%03d'%index
    arrays = reader.get(key,geometry_names)
    basis = restore_basis(arrays,reader.items[key]['metadata'],slot,control.grid)
    source = bound_json(lineage['source_state'])
    if source['archive'] != lineage['source_archive']:
        raise ValueError('original full raw D source archive differs')
    source_arrays = bound_arrays(source['archive'])
    authenticate(record.frame); authenticate(record.coefficient)
    psi,coefficient = original.arrays(record)
    try:
        expected_d = source_arrays['density_coefficients']
        if slot.orbital_conjugation:
            expected_d = expected_d.conj()
        if not np.array_equal(coefficient,expected_d):
            raise ValueError('original full raw D source transport differs')
        factors,compact = position_factors(inputs,basis,electronic_atoms)
        coupling = physical_coupling(inputs,basis,factors,block_size=8)
        coupling.cubic_radials = compact
        if electronic_atoms is not None:
            coupling.cubic_atom_indices = tuple(electronic_atoms)
        local,nonlocal_term = full_state_third(psi,coefficient,coupling)
        return dict(local=local,nonlocal_term=nonlocal_term)
    finally:
        release_mapped_pages(psi); release_mapped_pages(coefficient)


def admit_equivalent_algorithm(cache_path,specification,admission_path):
    """Explicitly admit only a proved contraction association, never physics/state."""
    admission_binding = file_binding(admission_path)
    admission = bound_json(admission_binding)
    previous_binding = file_binding(cache_path)
    previous = json.loads(Path(cache_path).read_bytes())
    old_spec = previous['specification']
    expected = dict(old_spec,source_bindings=specification['source_bindings'])
    if 'electronic_atom_indices' in specification:
        if admission.get('electronic_atom_indices') != specification['electronic_atom_indices']:
            raise ValueError('no exact explicit external-atom coverage admission')
        expected['electronic_atom_indices'] = specification['electronic_atom_indices']
    if (admission['schema'] != 'p254-a43-contraction-association-admission-v1'
            or admission['status'] != 'passed'
            or admission['original_partial_cache'] != previous_binding
            or admission['original_source_bindings'] != old_spec['source_bindings']
            or admission['current_source_bindings'] != specification['source_bindings']
            or expected != specification or previous['status'] != 'partial'
            or not previous['items'] or 'assembled' in previous['items']):
        raise ValueError('no exact source-bound association proof for partial-cache reuse')
    authenticate(admission['oracle'])
    for binding in admission['original_source_snapshots']:
        authenticate(binding)
    # The original manifest bytes are retained by the proof; existing NPZ items
    # and their original timings/attribution are never modified.
    authenticate(admission['original_partial_cache_snapshot'])
    previous['source_migrations'] = previous.get('source_migrations',[]) + [dict(admission=admission_binding,
        original_manifest=previous_binding,original_source_bindings=old_spec['source_bindings'],
        current_source_bindings=specification['source_bindings'],
        reused_items={key:item['archive'] for key,item in previous['items'].items()})]
    previous['specification'] = specification
    atomic_json(cache_path,previous)


def prepare_full_third(manifest,manifest_binding,output,*,resume=False,association_admission=None,electronic_atoms=None):
    """Commit immutable unweighted third jets; sum authenticated saved128 in order."""
    started = perf_counter(); output = Path(output).resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError('scientific arrays must remain outside repository')
    output.mkdir(parents=True,exist_ok=resume)
    control,inputs,rows,states,original,lower,context = restore_original(manifest)
    reader = ImmutableItems(manifest['items'])
    records = {record.identity:record for record in original.slots}
    expected = {str(slot.ispin)+':'+','.join(map(str,slot.mesh_index)) for _,slot,_ in rows}
    if set(records) != expected:
        raise ValueError('original full128 identities differ')
    electronic_atoms = list(range(inputs.reference.nat)) if electronic_atoms is None else list(electronic_atoms)
    if (not electronic_atoms or len(set(electronic_atoms)) != len(electronic_atoms)
            or min(electronic_atoms)<0 or max(electronic_atoms)>=inputs.reference.nat):
        raise ValueError('distinct explicit original electronic atoms required')
    lineage = []
    for index,(_,slot,_) in enumerate(rows):
        record = records[str(slot.ispin)+':'+','.join(map(str,slot.mesh_index))]
        lineage.append(_sector_lineage(index,slot,record,manifest['items']['preparation/%03d'%index],
                                      states[(slot.ispin,slot.source_ik)]))
    sources = [file_binding(HERE/name) for name in
               ('cubic_material_source.py','cubic_electronic_donor.py','cubic_position_moments.py')]
    specification = dict(kind='full-original-state-new-onsite-thirds-only',coefficient_manifest=manifest_binding,
        original_boundary=BOUNDARY,original_lower_order_authority=lower,original_context=context,
        material_sources=list(inputs.sources),matching=MATCHING,full_sector_lineage=lineage,
        source_bindings=sources,electronic_atom_indices=electronic_atoms)
    cache_path = output/'sector-thirds'/'manifest.json'
    if association_admission is not None:
        if not resume:
            raise ValueError('association admission applies only to exact partial resume')
        admit_equivalent_algorithm(cache_path,specification,association_admission)
    cache = ProjectionCache(output/'sector-thirds',specification)
    local = np.zeros((len(electronic_atoms),3,3,3),complex)
    nonlocal_term = np.zeros_like(local); newly_computed = reused = 0
    sector_archives = []
    for index,(_,slot,_) in enumerate(rows):
        unit_started = perf_counter(); key = 'sector/%03d'%index
        if cache.has(key):
            reused += 1
        else:
            record = records[lineage[index]['identity']]
            arrays = _new_sector(index,slot,record,lineage[index],reader,control,inputs,original,electronic_atoms)
            if any(value.shape != local.shape or not np.isfinite(value).all() for value in arrays.values()):
                raise ValueError('finite complete declared-atom original-sector thirds required')
            cache.put(key,arrays,metadata=dict(lineage[index],electronic_atom_indices=electronic_atoms),
                      timings=dict(seconds=perf_counter()-unit_started))
            newly_computed += 1
            del arrays
        metadata = cache.manifest['items'][key]['metadata']
        attribution = {name:value for name,value in metadata.items() if name!='electronic_atom_indices'}
        if attribution != lineage[index]:
            raise ValueError('cached third sector has different actual original source attribution')
        saved_atoms = metadata.get('electronic_atom_indices',list(range(inputs.reference.nat)))
        if any(atom not in saved_atoms for atom in electronic_atoms):
            raise ValueError('cached third sector lacks a requested original atom')
        saved = cache.get(key)
        if (set(saved) != {'local','nonlocal_term'}
                or any(value.shape != (len(saved_atoms),3,3,3) or not np.isfinite(value).all()
                       for value in saved.values())):
            raise ValueError('invalid authenticated saved sector-third payload')
        selected = [saved_atoms.index(atom) for atom in electronic_atoms]
        local += lineage[index]['weight']*saved['local'][selected]
        nonlocal_term += lineage[index]['weight']*saved['nonlocal_term'][selected]
        sector_archives.append(cache.manifest['items'][key]['archive'])
        print(json.dumps(dict(completed_sector=index,identity=lineage[index]['identity'],
            weight=lineage[index]['weight'],archive=sector_archives[-1],newly_computed=newly_computed,
            reused=reused,sector_elapsed_seconds=perf_counter()-unit_started,
            original_sector_compute_seconds=cache.manifest['timings'][key]['seconds'],
            elapsed_seconds=perf_counter()-started)),flush=True)
        del saved
    assembly_metadata = dict(ordered_sector_archives=sector_archives,
        electronic_atom_indices=electronic_atoms,weighting='original explicit weights applied once')
    if not cache.has('assembled'):
        cache.put('assembled',dict(local=local,nonlocal_term=nonlocal_term,
            electronic_atom_indices=np.asarray(electronic_atoms,int),
            cell=inputs.reference.cell_bohr,positions=inputs.reference.positions_bohr.T,
            charges=np.asarray([inputs.tables[n].upf.z_valence for n in inputs.reference.atom_names])),
            metadata=assembly_metadata)
    if cache.manifest['items']['assembled']['metadata'] != assembly_metadata:
        raise ValueError('assembled third binding differs from saved original128-sector sum')
    if cache.manifest['status'] != 'complete':
        cache.manifest['status'] = 'complete'; cache.publish()
    source_path = output/'source.json'
    record = dict(schema='p254-a43-full-state-third-v2',status='complete',
        coefficient_manifest=manifest_binding,original_boundary=BOUNDARY,
        original_lower_order_authority=lower,original_context=context,
        material_sources=list(inputs.sources),atom_names=list(inputs.reference.atom_names),
        matching=MATCHING,full_sector_lineage=lineage,electronic_atom_indices=electronic_atoms,
        electronic_atom_coverage='Only explicitly declared atoms; unavailable external readouts fail closed',
        arrays=cache.manifest['items']['assembled']['archive'],sector_cache=file_binding(cache.path),
        source_bindings=sources,original_settings_changed=False,
        recomputed_old_coefficients=False,computed_higher_jet_only=True)
    if cache.manifest.get('source_migrations'):
        record['source_migrations'] = cache.manifest['source_migrations']
    if source_path.exists():
        if not resume or json.loads(source_path.read_bytes()) != record:
            raise ValueError('completed source output differs from exact saved third-sector cache')
    else:
        atomic_json(source_path,record)
    return file_binding(source_path)
