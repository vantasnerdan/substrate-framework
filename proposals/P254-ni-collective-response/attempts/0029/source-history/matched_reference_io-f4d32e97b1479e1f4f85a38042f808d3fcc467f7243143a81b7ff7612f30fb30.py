"""Bound, resumable matched-control receipts; numerical arrays stay outside Git."""
import hashlib
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0025'), str(HERE.parent/'0024')]


def binding(path):
    path = Path(path).resolve()
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(path=str(path), bytes=str(path.stat().st_size), sha256=digest)


def save_json(path, value):
    path = Path(path)
    with path.open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')
    return binding(path)


def runtime(array_root=None):
    affinity = sorted(os.sched_getaffinity(0))
    if len(affinity) != 1 or affinity[0] not in range(6, 11):
        raise ValueError('allocate exactly one CPU in 6..10')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise ValueError(name+'=1 required before startup')
    if array_root is not None:
        path = Path(array_root).resolve()
        if path == ROOT or ROOT in path.parents:
            raise ValueError('numerical arrays must remain outside Git')
        path.mkdir(parents=True, exist_ok=True)


def registration():
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.qe_saved_native import read_bound_bytes
    record = binding(HERE/'matched-reference-registration.json')
    data = bound_json(record)
    for source in data['implementation_sources']:
        read_bound_bytes(source, max_bytes=32 << 20)
    return record, data


def load_case(case, data, *, background=True):
    from substrate_framework.material_control_reference import bound_json, bound_arrays, restore_control_background
    from substrate_framework.material_control_inputs import load_control_inputs
    from substrate_framework.material_matched_reference import ledger_slots
    selected = data['cases'][case]
    intake = bound_json(selected['reference'])
    if intake['case'] != case or intake['status'] != 'completed-background-restoration-and-bounded-H0-K-intake':
        raise ValueError('completed original H0/K intake required')
    prep_binding = intake['preparation']
    prep = bound_json(prep_binding)
    p = prep['provenance']
    control = load_control_inputs(p['case_manifest'], p['native_receipt'], p['primary_admission'])
    if control.case != case or control.provenance() != p:
        raise ValueError('original control provenance changed')
    ledger = bound_arrays(prep['population_ledger'])
    slots = ledger_slots(control.grid, ledger)
    bg = restore_control_background(prep_binding, intake['background_receipt'])[0] if background else None
    context = dict(case=case, model=data['model'], reference=selected['reference'],
                   preparation=prep_binding, population_ledger=prep['population_ledger'],
                   background_receipt=intake['background_receipt'])
    return control, bg, ledger, slots, context


def save_arrays(path, arrays, max_bytes):
    import numpy as np
    path = Path(path)
    if path.exists():
        raise ValueError('unreceipted archive already exists: '+str(path))
    if sum(np.asarray(a).nbytes for a in arrays.values()) > max_bytes:
        raise ValueError('explicit array byte cap exceeded')
    temporary = path.with_suffix('.partial.npz')
    with temporary.open('wb') as stream:
        np.savez(stream, **arrays)
    if temporary.stat().st_size > max_bytes:
        temporary.unlink()
        raise ValueError('archive byte cap exceeded')
    temporary.rename(path)
    return binding(path)


def receipt_arrays(record):
    from substrate_framework.material_control_reference import bound_json, bound_arrays
    value = bound_json(record)
    if value['status'] != 'complete':
        raise ValueError('incomplete coefficient receipt')
    return value, bound_arrays(value['archive'])


def loaded_field_cache(control, state, source, registration_data, *, expected_populations):
    """Reuse only exactly identical loaded identity-target optical/mixed means.

    Historical unitary and current real-policy basis identities have different
    namespaces. Prove equality through original XML/UPF/WFC content, actual
    identity target, all source rows and exact Gram/D/populations instead.
    Cached A25 spatial terms get their stored physical difference exactly once.
    """
    if control.case != 'loaded-ni31h6':
        return None, None
    import numpy as np
    from substrate_framework.material_control_reference import bound_json
    import run_full_field_moments_material as old
    old.preload_parsers()
    independent = registration_data['loaded_independent_evidence']
    _, manifest_binding, manifest, saved, original, grid = old.load_inputs(
        independent['conditional']['path'])
    def content(record):
        return int(record['bytes']), record['sha256']
    if (content(original['xml']) != content(control.xml_binding)
            or content(original['density']) != content(control.density_binding)
            or {k:content(v) for k,v in original['upfs'].items()} !=
               {k:content(v) for k,v in control.upf_bindings.items()}):
        raise ValueError('loaded cache geometry/UPF content differs')
    key = source.ispin, source.ik
    records = [r for r in manifest['sources'] if (r['spin'],r['source_ik']) == key]
    if len(records) != 1:
        raise ValueError('missing or duplicate loaded source state')
    record = records[0]
    if record['reuse_selected_gamma']:
        if key != (1,1):
            raise ValueError('unexpected cached Gamma state')
        arrays = saved.selected_gamma
        authority = manifest['preparation_terminal']
    else:
        arrays, state_receipt = old.receipt_arrays(record['state'])
        authority = record['state']
    if not np.array_equal(arrays['populations'],expected_populations):
        raise ValueError('original loaded source populations differ from the control ledger')
    if state is None:
        from substrate_framework.material_reference_state import restore_source
        state = restore_source({k:arrays[k] for k in
            ('gram','density_coefficients','populations','density_eigenvalues')},
            source_identity=source.source['sha256'])
    elif any(not np.array_equal(arrays[name],getattr(state,name))
             for name in ('gram','density_coefficients','populations')):
        raise ValueError('loaded original state differs; no silent cache reuse or recomputation')
    representatives = [s for s in grid.slots if (s.ispin,s.source_ik)==key
                       and s.operation.index==1 and not np.any(s.reciprocal_shift)]
    if len(representatives)!=1 or content(representatives[0].source_binding)!=content(source.source):
        raise ValueError('loaded identity representative/source differs')
    representative = representatives[0]
    if (not np.array_equal(representative.operation.rotation,np.eye(3))
            or np.any(representative.operation.fractional_translation)):
        raise ValueError('cached representative is not the original raw identity frame')
    joint = old.unpack(bound_json(independent['joint']))
    optical_aggregate = old.unpack(bound_json(joint['optical_aggregate']))
    if joint['status']!='complete' or optical_aggregate['status']!='complete':
        raise ValueError('completed loaded optical/mixed evidence required')
    selected = []
    for collection, field in ((optical_aggregate['source_checkpoints'],'moments'),
                               (joint['mixed_source_checkpoints'],'mixed')):
        matches = []
        for checkpoint in collection:
            value = old.unpack(bound_json(checkpoint))
            if (value['spin'],value['source_ik'])==key:
                matches.append((checkpoint,value))
        if len(matches)!=1:
            raise ValueError('loaded field cache source coverage differs')
        checkpoint,value = matches[0]
        if (value['status']!='complete' or value['mode']!='source'
                or value['manifest']!=manifest_binding or value['state_authority']!=authority
                or content(value['source_binding'])!=content(source.source)
                or value['original_source_weight_applied'] is not False):
            raise ValueError('loaded cache source/state/weight contract differs')
        selected.append((checkpoint,value[field]))
    import full_reference_inputs as spatial_reader
    density_receipt, spatial = spatial_reader.checkpoint(record['density_gradient'])
    curvature_receipt, curvature = spatial_reader.checkpoint(record['curvature'])
    dm = density_receipt['metadata']
    if not record['reuse_selected_gamma'] and any(
            state_receipt['metadata'].get(k)!=dm[k] for k in spatial_reader.SOURCE_FIELDS):
        raise ValueError('loaded source state and moments have different authorities')
    cm = (curvature_receipt['metadata']['source_metadata']
          if record['curvature_kind']=='expansion' else curvature_receipt['metadata'])
    if record['curvature_kind'] not in ('source','expansion'):
        raise ValueError('unknown loaded curvature codec')
    if (any(cm.get(k)!=dm[k] for k in spatial_reader.SOURCE_FIELDS)
            or (dm['spin'],dm['source_ik'])!=key or not dm['all_raw_columns']
            or content(dm['source'])!=content(source.source)
            or dm['preparation_terminal']!=manifest['preparation_terminal']
            or dm['coverage_identity']!=saved.coverage_identity
            or dm['mesh_index']!=list(representative.mesh_index)
            or not np.array_equal(curvature['conditional_density_bohr3'],spatial['density_bohr3'])):
        raise ValueError('loaded density/absolute nonlocal curvature authority differs')
    if record['curvature_kind']=='expansion' and (
            not np.array_equal(curvature['source_gram'],state.gram) or
            not np.array_equal(curvature['source_density_coefficients'],state.density_coefficients)):
        raise ValueError('loaded expansion contains different original Gram/D')
    spatial['conditional_nonlocal_ry_bohr2'] = curvature['conditional_nonlocal_ry_bohr2']
    reuse = dict(optical=selected[0][0], mixed=selected[1][0], state=authority,
                 density_gradient=record['density_gradient'], curvature=record['curvature'],
                 original_basis_identity=dm['original_basis_identity'],
                 original_frame_identity=dm['original_frame_identity'],
                 identity='exact XML/UPF/WFC content and Gram/D/populations; original identity target only',
                 cached_A25_spatial_conversion='stored physical-minus-A25 gradient and nonlocal onsite added ONCE; no correction to fresh physical factors')
    return state, (representative, (selected[0][1],selected[1][1],reuse,spatial))


def cache_for_slot(cache, slot):
    if cache is None:
        return None
    import numpy as np
    representative, coefficients = cache
    if (slot.mesh_index == representative.mesh_index and not slot.orbital_conjugation
            and slot.operation.index == 1 and not np.any(slot.reciprocal_shift)
            and np.array_equal(slot.operation.rotation,np.eye(3))
            and np.array_equal(slot.k_cart_bohr_inv,representative.k_cart_bohr_inv)
            and not np.any(slot.operation.fractional_translation)):
        return coefficients
    return None
