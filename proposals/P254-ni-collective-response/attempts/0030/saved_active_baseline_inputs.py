"""Necessary NEW Gamma baseline from the immutable complete A29 loaded12/27 run.

Only saved load-bearing Gamma coefficients and the saved t=0 unexcited state
are opened. The caller owns CPU/memory/time admission. Dense baseline output
exists in memory only; it is a new integration, NOT recovered adaptive history.
No finite-dose evolution, projection, preparation, Native/DFT or filter import.
"""
from dataclasses import dataclass
from hashlib import sha256
import inspect
import json
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
A29 = HERE.parent / '0029'
sys.path[:0] = [str(ROOT / 'src'), str(A29), str(HERE.parent / '0028')]

import numpy as np
from reduced_coefficients import ReducedCoefficients, QCoefficients
from reduced_electronic_dynamics import ActiveBandDynamics, BaselineState, BaselineSolution
from reduced_electronic_correlations import ActiveBandCorrelations
from reduced_material_readout import original_geometry
from reduced_nested_q_views import NestedQBanks
from reduced_projection_algebra import q_feature_geometry
from reduced_projection_storage import ImmutableItems
from substrate_framework.material_clamped_array import ClampedArray
from substrate_framework.streamed_material_state import authenticate, file_binding


@dataclass(frozen=True)
class SavedActiveBaselineInputs:
    coefficients: ReducedCoefficients
    baseline: BaselineSolution
    correlations: ActiveBandCorrelations
    coordinate: int
    source: dict
    work: dict


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _digest(value):
    return sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _read_json(binding):
    # Do not use scientific encode/decoded: Native integers remain Python ints.
    with authenticate(binding).open() as stream:
        return json.load(stream)


def _implementation_authority(binding, snapshots):
    """Authenticate exact historical bytes, never relabel current source as old."""
    try:
        authenticate(binding)
        return dict(original_binding=binding, authenticated_binding=binding,
                    method='original bound source bytes')
    except (FileNotFoundError, ValueError):
        key = (binding['path'], int(binding['bytes']), binding['sha256'])
        if key in snapshots:
            saved = snapshots[key]
            _require(int(saved['bytes']) == int(binding['bytes'])
                     and saved['sha256'] == binding['sha256'],
                     'explicit historical source snapshot differs')
            authenticate(saved)
            return dict(original_binding=binding, authenticated_binding=saved,
                        method='exact original_binding/archived_snapshot pair; not new execution provenance')
        original = Path(binding['path'])
        checksum = binding['sha256']
        _require(original.suffix == '.py' and len(checksum) == 64
                 and all(c in '0123456789abcdef' for c in checksum),
                 'changed source requires exact canonical Python source history')
        saved = dict(path=str(original.parent / 'source-history' /
                              (original.stem + '-' + checksum + original.suffix)),
                     bytes=binding['bytes'], sha256=checksum)
        authenticate(saved)
        return dict(original_binding=binding, authenticated_binding=saved,
                    method='retrospective exact source-history authentication; not original execution binding')


def _execution_bindings(result, specification, manifest):
    # The equations/storage used NOW must be the original admitted A29 bytes.
    # Historical recovery below authenticates provenance only; it never silently
    # executes a newer dynamics/correlation implementation in place of the old.
    snapshots = {}
    for entry in manifest.get('fullq', {}).get('original_producer_history', []):
        original, saved = entry['original_binding'], entry['archived_snapshot']
        key = (original['path'], int(original['bytes']), original['sha256'])
        snapshots[key] = saved
    authorities = {}
    for binding in (specification['implementation']
                    + specification['physical_implementation']
                    + result['consumer_sources']):
        key = (binding['path'], int(binding['bytes']), binding['sha256'])
        if key not in authorities:
            authorities[key] = _implementation_authority(binding, snapshots)
    by_path = {binding['path']: binding for binding in
               specification['implementation'] + result['consumer_sources']}
    names = ('run_reduced_material_dynamics.py', 'reduced_coefficients.py',
             'reduced_electronic_dynamics.py', 'reduced_electronic_correlations.py',
             'reduced_correlation_transpose.py', 'reduced_ode_initial_step.py',
             'reduced_material_readout.py', 'reduced_nested_q_views.py',
             'reduced_projection_storage.py', 'reduced_projection_algebra.py')
    executed = []
    for name in names:
        path = A29 / name
        expected = by_path[str(path)]
        authenticate(expected)
        executed.append(expected)
        # Imported modules cannot come from a shadow implementation on sys.path.
        module = sys.modules.get(path.stem)
        if module is not None:
            _require(Path(module.__file__).resolve() == path.resolve(),
                     'active numerical module is not the authenticated A29 source: ' + name)
    executed.extend(file_binding(path) for path in (
        Path(__file__), ROOT / 'src/substrate_framework/material_clamped_array.py',
        ROOT / 'src/substrate_framework/streamed_material_state.py'))
    return list(authorities.values()), executed


def prepare_saved_active_baseline(original_result_path, *, work_directory):
    """Return coefficients/baseline/correlations/coordinate/source/work.

    Makes exactly one NEW unexcited Gamma integration [0,30], using the actual
    result's tolerances, authenticated producer's declared DOP853 default, and
    literal baseline-000 state. A Gamma-only mapped working view is isolated
    below work_directory/saved-active-baseline-inputs; it cannot populate the
    immutable original result/projection namespaces. Use a fresh work directory
    for this admitted pilot, not this API to replay the finite-dose campaign.
    """
    started, cpu_started = perf_counter(), process_time()
    original_binding = file_binding(original_result_path)
    result = _read_json(original_binding)
    _require(result.get('schema') == 'P254-activeband-material-run-v1'
             and result.get('status') == 'complete', 'complete original A29 finite result required')
    config = result['configuration']
    _require(_digest(config) == result['identity'], 'original result/configuration identity differs')
    specification = config['coefficient_specification']
    _require(specification['case'] == 'loaded-ni31h6' and config['bands'] == 12
             and config['features'] == 27 and not config['baseline_only']
             and config['atom'] == 1 and config['axis'] == 'x'
             and config['times'] == [0., 1., 10., 30.],
             'original loaded12/27 Ni1/x complete0/1/10/30 finite result required')
    records = result['completed']
    _require([record['time_inverse_ry'] for record in records] == config['times'],
             'complete original requested-time coverage required')
    snapshot_binding = result['original_coefficient_manifest_snapshot']
    input_binding = result['original_input_manifest']
    _require(snapshot_binding['sha256'] == input_binding['sha256']
             and int(snapshot_binding['bytes']) == int(input_binding['bytes']),
             'original input and immutable snapshot differ')
    authenticate(input_binding)
    manifest = _read_json(snapshot_binding)
    _require(manifest['schema'] == 'p254-physical-activeband-cache-v1'
             and manifest['status'] == 'complete'
             and manifest['specification'] == specification,
             'original complete coefficient snapshot/specification required')
    _require(specification['bands'] == 12 and specification['features'] == 27
             and specification['mesh'] == [4, 4, 4]
             and specification['case_identity'] == _digest(specification['context']),
             'original maximum loaded model identity differs')
    source_authorities, execution = _execution_bindings(result, specification, manifest)
    authenticate(specification['original_boundary'])  # Receipt only; no full frames/D are loaded.

    # Same fixed-array/storage path as load_projected_coefficients, but no
    # all-q authentication and no mapper rooted at an immutable input manifest.
    keys = (['slot/%03d' % slot for slot in range(128)] + ['mechanical']
            + ['edge/00/%03d' % slot for slot in range(128)] + ['kernel/00'])
    _require(all(key in manifest['items'] for key in keys)
             and 'passive/12' in manifest['items'], 'load-bearing Gamma inputs/passive split absent')
    items = {key: manifest['items'][key] for key in keys}
    source_archives = {key: item['archive'] for key, item in items.items()}
    initial_binding = records[0]['baseline_state']
    geometry_binding = config['geometry_archive']
    work_root = Path(work_directory).resolve()
    immutable_parents = {Path(binding['path']).resolve().parent for binding in
                         [original_binding, snapshot_binding, input_binding,
                          initial_binding, geometry_binding] + list(source_archives.values())}
    _require(not any(work_root.is_relative_to(parent) for parent in immutable_parents),
             'work directory must be outside original result/coefficient/source namespaces')
    directory = work_root / 'saved-active-baseline-inputs'
    directory.mkdir(parents=True, exist_ok=False)
    # Derivative receipt: original spec and original NPZ bindings remain exact.
    # Remove unneeded repeated item metadata, not any scientific coefficient.
    view = dict(schema=manifest['schema'], status='gamma-complete',
                specification=specification,
                items={key: dict(archive=item['archive']) for key, item in items.items()})
    view_path = directory / 'gamma-input-manifest.json'
    with view_path.open('x') as stream:
        json.dump(view, stream, sort_keys=True, allow_nan=False)
        stream.write('\n')
    view_binding = file_binding(view_path)
    reader = ImmutableItems(items)
    fixed = [reader.get('slot/%03d' % slot, ('h0', 'rho0', 'field', 'field_reverse'))
             for slot in range(128)]
    mechanical_arrays = reader.get('mechanical')
    mechanical = ClampedArray(**mechanical_arrays)
    banks = NestedQBanks(view_path, view, 12, 27, QCoefficients, q_feature_geometry)

    def gamma_only(index):
        if index != 0:
            raise ValueError('this necessary baseline/filter working view contains Gamma only')
        return banks(0)

    indexing = specification['indexing']
    coefficients = ReducedCoefficients(
        mesh=tuple(specification['mesh']), k_indices=np.asarray(indexing['k_indices']),
        q_indices=np.asarray(indexing['q_indices']), slot_spin=np.asarray(indexing['slot_spin']),
        slot_k=np.asarray(indexing['slot_k']), target=np.asarray(indexing['target']),
        minus=np.asarray(indexing['minus']), h0=np.asarray([item['h0'] for item in fixed]),
        rho0=np.asarray([item['rho0'] for item in fixed]),
        field=np.asarray([item['field'] for item in fixed]),
        field_reverse=np.asarray([item['field_reverse'] for item in fixed]),
        mechanical=mechanical, q_provider=gamma_only,
        field_polarizations=mechanical_arrays['polarizations'])
    with np.load(authenticate(geometry_binding), allow_pickle=False) as archive:
        geometry = original_geometry({key: archive[key] for key in
            ('atom_indices', 'atom_names', 'cell_bohr', 'valence_charges')}, coefficients.mesh)
    selected = np.flatnonzero(geometry['ids'] == config['atom'])
    _require(len(selected) == 1 and geometry['names'][selected[0]] == 'Ni'
             and geometry['d'] == coefficients.d, 'original Ni1 geometry/coordinate differs')
    coordinate = int(3 * int(selected[0]) + 'xyz'.index(config['axis']))
    with np.load(authenticate(initial_binding), allow_pickle=False) as archive:
        _require(set(archive.files) == {'rho', 'x'}, 'literal saved baseline0 rho/x required')
        initial = BaselineState(archive['rho'], archive['x'])
    model = ActiveBandDynamics(coefficients)
    model.pack_baseline(initial)  # Original shape contract, without initial-state replacement.
    _require(np.array_equal(initial.rho, coefficients.rho0) and not np.any(initial.x)
             and np.isfinite(initial.rho).all() and np.isfinite(initial.x).all(),
             'saved original unexcited initial state differs from coefficient initial state')
    method = inspect.signature(ActiveBandDynamics.integrate_baseline).parameters['method'].default
    _require(method == 'DOP853' and config.get('method', method) == method,
             'original producer declared DOP853 solver required')
    rtol, atol = config['rtol'], config['atol']
    _require(np.isfinite(rtol) and rtol > 0 and np.isfinite(atol) and atol > 0,
             'original finite positive solver tolerances required')
    cache_root = banks._maximum.root
    cache_receipt = cache_root / 'q00/completed.json'
    cache_binding = file_binding(cache_receipt)
    cache = _read_json(cache_binding)
    _require(cache['source_archives'] == [source_archives['edge/00/%03d' % slot]
             for slot in range(128)] + [source_archives['kernel/00']],
             'Gamma working cache source bindings differ')
    source = dict(
        schema='P254-saved-active-baseline-inputs-v1', original_result=original_binding,
        original_result_identity=result['identity'], original_coefficient_manifest_snapshot=snapshot_binding,
        original_input_manifest=input_binding, geometry_archive=geometry_binding,
        original_boundary=specification['original_boundary'], baseline0=initial_binding,
        coefficient_source_archives=source_archives,
        original_initial_coefficient_archive=manifest['items']['initial']['archive'],
        original_passive_split_archive=manifest['items']['passive/12']['archive'],
        initial_and_passive_archive_usage='identities retained from original manifest; neither archive is opened or used to redo preparation/projection; literal baseline0 supplies the initial state',
        original_consumer_sources=result['consumer_sources'],
        original_implementation_authorities=source_authorities, execution_implementation=execution,
        case='loaded-ni31h6', bands=12, features=27, q_index=0,
        atom=1, axis='x', coordinate=coordinate, interval_inverse_ry=[0., 30.],
        rtol=rtol, atol=atol, method=method,
        method_authority='authenticated original producer call uses authenticated integrate_baseline declared default; method was not a saved configuration key',
        history='NEW necessary unexcited Gamma baseline reconstruction from unchanged original inputs; NOT recovered original dense/adaptive history',
        scope='conditional approved moving active model; not a new complete material/full-PW simulation, passive noise, loss or lifetime')
    # Release the enormous all-q receipt and unused arrays BEFORE dense integration.
    del manifest, result, records, fixed, mechanical_arrays, geometry, items, reader, indexing
    loading_seconds = perf_counter() - started
    baseline_started = perf_counter()
    baseline = model.integrate_baseline((0., 30.), initial=initial,
                                        rtol=rtol, atol=atol, method=method)
    baseline_seconds = perf_counter() - baseline_started
    correlations = ActiveBandCorrelations(coefficients, baseline, rtol=rtol, atol=atol, method=method)
    work_files = [file_binding(path) for path in sorted(directory.rglob('*')) if path.is_file()]
    work = dict(
        schema='P254-new-necessary-baseline-work-v1', directory=str(directory),
        gamma_input_manifest=view_binding, gamma_cache_root=str(cache_root),
        gamma_cache_receipt=cache_binding, gamma_cache_arrays=cache['arrays'],
        files=work_files, total_disk_bytes=sum(int(binding['bytes']) for binding in work_files),
        new_baseline_integrations=1, finite_dose_integrations=0, non_gamma_banks_loaded=0,
        source_npz_archives_used=len(source_archives), finite_state_archives_opened=0,
        dense_state_history_written=False, original_namespaces_written=False,
        loading_seconds=loading_seconds, baseline_seconds=baseline_seconds,
        elapsed_seconds=perf_counter() - started, cpu_seconds=process_time() - cpu_started,
        accepted_baseline_intervals=int(len(baseline.result.t) - 1),
        baseline_rhs_evaluations=int(baseline.result.nfev),
        retained_baseline_boundary_bytes=int(baseline.result.y.nbytes),
        peak_process_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024,
        history='NEW necessary baseline reconstruction; in-memory dense output is not original adaptive-history recovery',
        resource_scope='caller must enforce admitted CPU7/BLAS1, MemoryMax4GiB/swap0/RuntimeMax600s and shared5GiB pool; Gamma mmap/cache plus in-memory adaptive baseline; no all-q/full-state/action bank',
        scientific_scope='one conditional active-model baseline only; not a new complete material/full-PW simulation')
    return SavedActiveBaselineInputs(coefficients, baseline, correlations, coordinate, source, work)
