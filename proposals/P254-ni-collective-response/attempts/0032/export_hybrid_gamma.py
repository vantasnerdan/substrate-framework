"""Run one NEW moving source0/Gamma hybrid from immutable paid saved inputs.

CLI: --bath-result --original-result --output (a fresh directory). No physical
Hamiltonian/force action, DFT, old baseline solve, state refill or projection.
"""
import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
A29 = HERE.parent/'0029'
A30 = HERE.parent/'0030'
sys.path[:0] = [str(ROOT/'src'), str(A29), str(A30), str(HERE.parent/'0028')]

import numpy as np
import scipy
from reduced_coefficients import ReducedCoefficients, QCoefficients
from reduced_material_readout import original_geometry
from reduced_nested_q_views import NestedQBanks
from reduced_projection_algebra import q_feature_geometry
from reduced_projection_storage import ImmutableItems
from saved_active_baseline_inputs import _execution_bindings
from substrate_framework.material_clamped_array import ClampedArray
from substrate_framework.streamed_material_state import authenticate, file_binding
from hybrid_gamma_evolution import evolve_hybrid_gamma, RECOVERY, save_npz_exclusive


SCOPE = ('ONE Gamma moving Hamiltonian/Gaussian model: source0 replaced once by '
         'the paid354-dimensional finite bath; other127 original12-band sectors '
         'remain dynamic; original111 ionic and3 field coordinates/preparation '
         'unchanged; source0 density/field active embedding plus two paid ionic '
         'force-mode corrections. Not allq/allsource/full-PW physical calibration, '
         'spatial coherence correction, physical loss, damping, lifetime, '
         'larger-basis qualification or campaign completion. Original16 runs unchanged.')


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _same(left, right):
    return (Path(left['path']).resolve() == Path(right['path']).resolve()
            and int(left['bytes']) == int(right['bytes'])
            and left['sha256'] == right['sha256'])


def _digest(value):
    return sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _json(binding):
    with authenticate(binding).open() as stream:
        return json.load(stream)


def _arrays(binding):
    with np.load(authenticate(binding), allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def _encoded(value):
    """Preserve complex parts; represent undefined numbers explicitly in JSON."""
    if isinstance(value, dict):
        return {key: _encoded(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_encoded(item) for item in value]
    if isinstance(value, np.ndarray):
        return _encoded(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=_encoded(float(value.real)), imag=_encoded(float(value.imag)))
    if isinstance(value, np.generic):
        return _encoded(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return 'NaN' if np.isnan(value) else ('+Infinity' if value > 0 else '-Infinity')
    return value


def publish_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(_encoded(value), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    return file_binding(path)


@dataclass(frozen=True)
class SavedHybridInputs:
    coefficients: ReducedCoefficients
    initial_mechanical: np.ndarray
    bath_bank: dict
    mode_directions: np.ndarray
    coordinate: int
    source: dict
    work: dict


def restore_saved_hybrid_inputs(bath_result_path, original_result_path, *, output):
    """Restore the A30 fixed/Gamma cache path, WITHOUT integrating its baseline.

    The fresh derivative Gamma manifest isolates the old mapper's working files.
    Its arrays still authenticate against the original immutable archives.
    """
    started = perf_counter()
    original_binding = file_binding(original_result_path)
    original = _json(original_binding)
    _require(original.get('schema') == 'P254-activeband-material-run-v1'
             and original.get('status') == 'complete', 'complete original A29 result required')
    config = original['configuration']
    _require(_digest(config) == original['identity'], 'original configuration identity differs')
    spec = config['coefficient_specification']
    _require(spec['case'] == 'loaded-ni31h6' and config['bands'] == 12
             and config['features'] == 27 and not config['baseline_only']
             and config['atom'] == 1 and config['axis'] == 'x'
             and config['times'] == [0., 1., 10., 30.], 'original loaded12/27 Ni1/x0/1/10/30 required')
    _require([item['time_inverse_ry'] for item in original['completed']] == config['times'],
             'original complete requested-time coverage required')
    snapshot = original['original_coefficient_manifest_snapshot']
    input_manifest = original['original_input_manifest']
    _require(snapshot['sha256'] == input_manifest['sha256']
             and int(snapshot['bytes']) == int(input_manifest['bytes']), 'input snapshot differs')
    authenticate(input_manifest)
    manifest = _json(snapshot)
    _require(manifest['schema'] == 'p254-physical-activeband-cache-v1'
             and manifest['status'] == 'complete' and manifest['specification'] == spec,
             'original complete coefficient manifest/specification required')
    _require(spec['bands'] == 12 and spec['features'] == 27 and spec['mesh'] == [4,4,4]
             and spec['case_identity'] == _digest(spec['context']), 'original maximum model differs')
    authorities, execution = _execution_bindings(original, spec, manifest)
    boundary = _json(spec['original_boundary'])
    _require(boundary['time_ry_inverse'] == 0 and len(boundary['slots']) == 128,
             'original complete128 t0 boundary required')
    baseline0 = original['completed'][0]['baseline_state']
    initial = _arrays(baseline0)
    _require(set(initial) == {'rho', 'x'}, 'literal saved baseline0 rho/x required')
    keys = (['slot/%03d' % slot for slot in range(128)] + ['mechanical']
            + ['edge/00/%03d' % slot for slot in range(128)] + ['kernel/00'])
    _require(all(key in manifest['items'] for key in keys)
             and 'passive/12' in manifest['items'], 'complete load-bearing Gamma/passive inputs required')
    items = {key: manifest['items'][key] for key in keys}
    archives = {key: item['archive'] for key, item in items.items()}

    bath_binding = file_binding(bath_result_path)
    contract_binding = file_binding(HERE/'implementation-contract-001.json')
    contract = _json(contract_binding)
    _require(_same(bath_binding, contract['existing_paid_data']['bath_result']),
             'exact preregistered paid A31 bath result required')
    bath = _json(bath_binding)
    _require(bath.get('schema') == 'P254-A31-actual-frozen-bath-comparator-v1'
             and bath.get('status') == 'complete', 'complete physical finite-bath result required')
    source = bath['source']
    action = _json(bath['original_action_packet'])
    mode1 = _json(bath['new_mode_packet'])
    filter_packet = _json(bath['filter_packet'])
    _require(action['status'] == 'complete' and mode1['status'] == 'complete'
             and _same(filter_packet['source']['original_result'], original_binding)
             and _same(filter_packet['source']['baseline0'], baseline0)
             and _same(filter_packet['source']['original_boundary'], spec['original_boundary'])
             and _same(action['filter_packet'], bath['filter_packet'])
             and _same(mode1['original_action_packet'], bath['original_action_packet']),
             'paid bath/action/mode/filter/original-result lineage differs')
    metadata = items['slot/000']['metadata']
    record = boundary['slots'][0]
    _require(metadata['slot'] == 0 and metadata['identity'] == record['identity'] == '1:0,0,0'
             and source['original_record'] == 0 and source['coefficient_item'] == 'slot/000'
             and source['spin'] == 1 and source['mesh_index'] == [0,0,0]
             and spec['indexing']['slot_spin'][0] == 1 and spec['indexing']['slot_k'][0] == 0,
             'exact source0/Gamma spin/slot identity differs')
    for left, right in ((source['coefficient_archive'], archives['slot/000']),
                        (source['coefficient_manifest'], input_manifest),
                        (source['original_boundary'], spec['original_boundary']),
                        (metadata['frame'], record['frame']),
                        (metadata['raw_d'], record['coefficient']),
                        (source['original_raw_d'], record['coefficient']),
                        (action['frame'], record['frame']),
                        (source['coordinates'], action['coordinates']),
                        (source['source_state'], metadata['source_state']),
                        (source['source_archive'], metadata['source_archive']),
                        (source['source'], metadata['source'])):
        _require(_same(left, right), 'source0 original saved binding differs')
        authenticate(left)
    for name in ('selected', 'orbital_conjugation', 'operation', 'reciprocal_shift', 'source_admission'):
        _require(source[name] == metadata[name], 'source0 saved lineage differs: '+name)
    _require(source['coefficient_case_identity'] == spec['case_identity']
             and metadata['basis_identity'] == record['basis_identity']
             and metadata['frame_identity'] == record['gauge_identity'], 'source0 model/basis/gauge differs')
    _require(source['original_metadata'] == {key:value for key,value in record.items()
             if key not in ('frame', 'coefficient')}, 'original source0 boundary metadata differs')
    state = _json(metadata['source_state'])
    _require(_same(state['archive'], metadata['source_archive'])
             and _same(spec['context']['source_states'][str((1,1))], metadata['source_state']),
             'original physical source-state authority differs')
    authority = _digest(dict(initial=source['initial_result'], preparation=source['preparation'],
        boundary=source['original_boundary'], record=source['original_metadata'],
        frame=record['frame'], raw_d=record['coefficient'], state=source['source_state'],
        archive=state['archive']))
    _require(authority == bath['original_state_authority'] == action['state_authority']
             and record['basis_identity'] == bath['original_pw_basis_identity'] == action['pw_basis_identity']
             and action['factor_identity'] == bath['original_factor_identity'],
             'physical original state authority differs')
    for name in ('initial_result', 'preparation'):
        authenticate(source[name])
    bank = _arrays(bath['basis'])
    coordinates = _arrays(source['coordinates'])
    state_arrays = _arrays(state['archive'])
    raw_d = np.load(authenticate(record['coefficient']), allow_pickle=False)
    _require(not metadata['orbital_conjugation'], 'this settled source0 bank is direct unconjugated')
    _require(np.array_equal(bank['density_coefficients'], raw_d)
             and np.array_equal(raw_d, state_arrays['density_coefficients'])
             and np.array_equal(bank['gram'], state_arrays['gram'])
             and np.array_equal(coordinates['gram'], bank['gram'])
             and np.array_equal(coordinates['density_coefficients'], raw_d)
             and np.array_equal(coordinates['selected'], metadata['selected']),
             'unchanged original S/D/selection required')
    projected = bank['projected_frame']
    q = bank['projected_active_coordinates']
    embedding_defect = float(np.max(np.abs(q-projected@coordinates['active_coordinates'])))
    density_defect = float(np.max(np.abs(bank['projected_density']-projected@raw_d@projected.conj().T)))
    channels = _arrays(action['channel_archive'])
    directions = channels['ionic_basis'][:, :2]
    saved_mode1 = np.asarray(mode1['direction_real'])+1j*np.asarray(mode1['direction_imag'])
    _require(directions.shape == (111,2) and np.array_equal(directions[:,1], saved_mode1),
             'exact saved two-mode complex directions required')
    for name in ('hamiltonian', 'mode0', 'mode1'):
        entry = bath['basis_metadata']['operators'][name]
        _require(entry['pw_basis_identity'] == record['basis_identity']
                 and entry['source_authority'] == authority, 'bank operator state/basis differs')
    _require(mode1['actions']['source_authority'] == authority, 'mode1 authority differs')

    work_root = Path(output).resolve()
    protected = [original_binding, snapshot, input_manifest, baseline0, config['geometry_archive'],
                 bath_binding, bath['basis'], bath['original_action_packet'], bath['new_mode_packet'],
                 bath['filter_packet'], action['channel_archive']] + list(archives.values())
    _require(not any(work_root.is_relative_to(Path(binding['path']).resolve().parent)
                     for binding in protected), 'fresh output must be outside immutable input namespaces')
    directory = work_root/'saved-hybrid-inputs'
    directory.mkdir(exist_ok=False)
    view = dict(schema=manifest['schema'], status='gamma-complete', specification=spec,
                items={key:dict(archive=item['archive']) for key,item in items.items()})
    view_path = directory/'gamma-input-manifest.json'
    view_binding = publish_json(view_path, view)
    reader = ImmutableItems(items)
    fixed = [reader.get('slot/%03d' % slot, ('h0','rho0','field','field_reverse')) for slot in range(128)]
    mechanical_arrays = reader.get('mechanical')
    banks = NestedQBanks(view_path, view, 12, 27, QCoefficients, q_feature_geometry)

    def gamma_only(index):
        if index != 0:
            raise ValueError('oneGamma saved-input view forbids nonGamma banks')
        return banks(0)

    indexing = spec['indexing']
    coefficients = ReducedCoefficients(mesh=tuple(spec['mesh']),
        k_indices=np.asarray(indexing['k_indices']), q_indices=np.asarray(indexing['q_indices']),
        slot_spin=np.asarray(indexing['slot_spin']), slot_k=np.asarray(indexing['slot_k']),
        target=np.asarray(indexing['target']), minus=np.asarray(indexing['minus']),
        h0=np.asarray([item['h0'] for item in fixed]), rho0=np.asarray([item['rho0'] for item in fixed]),
        field=np.asarray([item['field'] for item in fixed]),
        field_reverse=np.asarray([item['field_reverse'] for item in fixed]),
        mechanical=ClampedArray(**mechanical_arrays), q_provider=gamma_only,
        field_polarizations=mechanical_arrays['polarizations'])
    _require(np.array_equal(initial['rho'], coefficients.rho0)
             and initial['x'].shape == (coefficients.mechanical_sizes[0],)
             and not np.any(initial['x']) and np.isfinite(initial['x']).all(),
             'literal original unexcited mechanical/electronic state differs')
    active_density_defect = float(np.max(np.abs(q.conj().T@bank['projected_density']@q-coefficients.rho0[0])))
    with np.load(authenticate(config['geometry_archive']), allow_pickle=False) as archive:
        geometry = original_geometry({key:archive[key] for key in
            ('atom_indices','atom_names','cell_bohr','valence_charges')}, coefficients.mesh)
    chosen = np.flatnonzero(geometry['ids'] == 1)
    _require(len(chosen) == 1 and geometry['names'][chosen[0]] == 'Ni' and geometry['d'] == 111,
             'original Ni1/x coordinate differs')
    coordinate = 3*int(chosen[0])
    _require(config.get('method', 'DOP853') == 'DOP853', 'original DOP853 required')
    rtol, atol = config['rtol'], config['atol']
    _require(np.isfinite(rtol) and rtol > 0 and np.isfinite(atol) and atol > 0,
             'finite positive original tolerances required')
    cache_binding = file_binding(banks._maximum.root/'q00/completed.json')
    cache = _json(cache_binding)
    _require(cache['source_archives'] == [archives['edge/00/%03d' % slot] for slot in range(128)]
             + [archives['kernel/00']], 'Gamma cache source bindings differ')
    bindings = dict(schema='P254-A32-saved-hybrid-inputs-v1', original_result=original_binding,
        original_result_identity=original['identity'], original_coefficient_manifest_snapshot=snapshot,
        original_input_manifest=input_manifest, baseline0=baseline0,
        original_boundary=spec['original_boundary'], coefficient_source_archives=archives,
        original_initial_archive=manifest['items']['initial']['archive'],
        original_passive_archive=manifest['items']['passive/12']['archive'],
        initial_passive_usage='identities only, neither reopened to redo preparation or projection',
        bath_result=bath_binding, basis=bath['basis'], mode_packet=bath['new_mode_packet'],
        action_packet=bath['original_action_packet'], direction_archive=action['channel_archive'],
        filter_packet=bath['filter_packet'], implementation_contract=contract_binding,
        original_source0=source, state_authority=authority, source_slot=0, q_index=0,
        coordinate=coordinate, atom=1, axis='x', method='DOP853', rtol=rtol, atol=atol,
        method_authority='authenticated A29 producer/integrate_baseline DOP853 default carried by A30 execution binding',
        transport_defects=dict(active_embedding_max_abs=embedding_defect,
            full_density_max_abs=density_defect, active_density_max_abs=active_density_defect),
        original_implementation_authorities=authorities, restored_implementation=execution,
        scope=SCOPE, initial_mechanical_convention='literal saved A29 canonical x unchanged')
    files = [file_binding(path) for path in sorted(directory.rglob('*')) if path.is_file()]
    work = dict(loading_seconds=perf_counter()-started, gamma_view=view_binding,
        gamma_cache_receipt=cache_binding, files=files,
        total_disk_bytes=sum(int(binding['bytes']) for binding in files),
        old_baseline_integrations=0, physical_hamiltonian_actions=0, physical_force_actions=0,
        dft_actions=0, non_gamma_banks_loaded=0, original_namespaces_written=False)
    return SavedHybridInputs(coefficients, initial['x'], bank, directions, coordinate, bindings, work)


def _implementation_bindings():
    paths = {Path(__file__).resolve(), HERE/'hybrid_gamma_evolution.py', HERE/'hybrid_gamma_model.py'}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return [file_binding(path) for path in sorted(paths)]


def run_hybrid_gamma(bath_result, original_result, output):
    """Fresh exclusive actual run. Failure leaves all prior checkpoints intact."""
    directory = Path(output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    started, cpu_started = perf_counter(), process_time()
    log_path = directory/'progress.jsonl'
    with log_path.open('x', buffering=1) as log:
        def progress(stage, **details):
            record = dict(stage=stage, utc=datetime.now(timezone.utc).isoformat(),
                          elapsed_seconds=perf_counter()-started, **details)
            text = json.dumps(_encoded(record), sort_keys=True, allow_nan=False)
            log.write(text+'\n')
            log.flush()
            os.fsync(log.fileno())
            print(text, flush=True)

        try:
            progress('restoring-saved-inputs', recovery=RECOVERY, scope=SCOPE)
            inputs = restore_saved_hybrid_inputs(bath_result, original_result, output=directory)
            # Delayed import keeps input restoration independent of the model;
            # implementation work never needs to import an unsettled model.
            from hybrid_gamma_model import HybridGammaModel
            model = HybridGammaModel(inputs.coefficients, inputs.initial_mechanical,
                                     inputs.bath_bank, inputs.mode_directions, source_slot=0)
            coordinates = save_npz_exclusive(directory/'interaction-coordinates.npz',
                eigenvalues=model.energies, right_eigenvectors=model.v,
                inverse_right_eigenvectors=model.vi)
            implementation = _implementation_bindings()
            source_receipt = publish_json(directory/'input-bindings.json', dict(
                source=inputs.source, implementation=implementation, loading=inputs.work,
                model_metadata=model.metadata, model_diagnostics=model.diagnostics,
                interaction_coordinates=coordinates,
                software=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
                recovery=RECOVERY, scope=SCOPE))
            progress('inputs-restored', input_bindings=source_receipt,
                     state_size=int(model.size), old_baseline_integrations=0)
            evolution = evolve_hybrid_gamma(model, inputs.coordinate, output=directory,
                rtol=inputs.source['rtol'], atol=inputs.source['atol'], progress=progress,
                publish_json=lambda path, value: publish_json(path, dict(
                    value, input_bindings=source_receipt, scope=SCOPE)))
            result = dict(schema='P254-A32-moving-hybrid-gamma-v1', status='complete',
                source=inputs.source, implementation=implementation, input_bindings=source_receipt,
                model_metadata=model.metadata, model_diagnostics=model.diagnostics,
                artifacts=evolution.artifacts, diagnostics=evolution.diagnostics,
                times_inverse_ry=[0.,1.,10.,30.], atom=1, axis='x', source_slot=0, q_index=0,
                ordered_covariance=evolution.ordered,
                normalized_gamma_coherence=evolution.coherence,
                coherence_abs=np.abs(evolution.coherence), coherence_phase=np.angle(evolution.coherence),
                work=dict(loading=inputs.work, elapsed_seconds=perf_counter()-started,
                    cpu_seconds=process_time()-cpu_started,
                    peak_process_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024,
                    new_hybrid_forward_mean_integrations=1, backward_mean_dual_passes=1,
                    dense_history_saved=False, old_baseline_integrations=0,
                    resource_scope='actual CPU/memory/runtime limits are recorded by the caller startup and Native outcome; no automatic restart'),
                normalization='Gamma-only ordered covariance and normalized coherence; no spatial Fourier/source/spin factor',
                recovery=RECOVERY, scope=SCOPE)
            result_binding = publish_json(directory/'result.json', result)
            progress('complete', result=result_binding)
            return result
        except BaseException as error:
            # Never delete a completed baseline because the Gaussian pass failed.
            progress('failed', exception=type(error).__name__, message=str(error),
                     baseline_complete=(directory/'baseline-complete.json').is_file(), recovery=RECOVERY)
            publish_json(directory/'failure.json', dict(status='failed',
                exception=type(error).__name__, message=str(error),
                baseline_complete=(directory/'baseline-complete.json').is_file(),
                elapsed_seconds=perf_counter()-started, recovery=RECOVERY, scope=SCOPE))
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bath-result', required=True)
    parser.add_argument('--original-result', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    run_hybrid_gamma(arguments.bath_result, arguments.original_result, arguments.output)


if __name__ == '__main__':
    main()
