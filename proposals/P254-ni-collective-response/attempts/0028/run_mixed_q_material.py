#!/usr/bin/env python3
"""One original selected-q raw-coherent Native mixed trajectory, never all64q.

Required authority: --case, --initialization/--initialization-sha256,
--history/--history-sha256, --execution-proof/--execution-proof-sha256, and
repeatable --frame-admission PATH SHA256 naming completed actual-QWork frame
admissions. The executed canonical local-pulse Fourier source fixes amplitude/N;
there is no additional per-run normalization-proof gate.
Select --q QX QY QZ --atom N --axis x|y|z --amplitude J --cell RX RY RZ,
--time T --subdivisions N, --cpus (unique CPUs6..10), --max-actions,
--max-seconds, --max-output-bytes (owned live logical paths), --max-memory-bytes
and --worker-memory-bytes. Readonly inputs/aliases are separately inventoried.

--mode plan performs metadata accounting only. --mode forward authors sparse
raw-zero jets and calls the genuine Native mixed RK4 producer. Available
admitted frames are reused exactly/read-only. Genuinely missing interior times
use the original qualified history/QWork exponential within the SAME caps;
owned frames retire after their final committed step/readout consumer.
There is no symmetry inference or extrapolation. --energy requests the fixed
mixed_energy_readout API at the final physical state; its complex eleven-term
mixed derivatives are NOT finite-dose energy, noise, heat or covariance labels.
The default rolls nonzero raw endpoint payloads: an obsolete step base is
retired only after its next actual endpoint and requested readouts are durable,
using qualified per-edge successor journals. All source/checkpoint/small/cost
metadata and the original sparse-zero seed survive. --retain-endpoints is an
explicit alternative whose full retained payloads are charged in the plan.

--directory is exclusive and --output must be a NEW JSON path inside it.
--resume reuses the same authenticated owner/source/input/caps/step identities,
original absolute wall deadline and durable action admissions. Choose a new
--output path on each resume; earlier partial/error receipts are immutable.
An invocation that dies without a terminal receipt retains its owner, original
deadline, action ledger, native/QWork attempts and source-bound step journals.
Root, not this program's metadata plan, owns all production authorization and
changed-path numerical/storage/refinement admission.
"""
import argparse
from contextlib import contextmanager
from copy import copy
from datetime import datetime, timezone
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import resource
import signal
import sys
from time import perf_counter, process_time, time
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from matched_joint_case import load_case
from matched_joint_state import preparation_arrays
from run_joint_q_material import bindings as existing_dependencies, admit_oracle, bind_saved_model, selected_initialization
from mixed_material_driver_support import (
    ContinuationBudget, RetainedFrameWork, canonical,
    check_loaded_sources, checked_json, exact_publish, frame_requirements, initial_small,
    initialize_raw_zero, metadata_geometry, nested_bindings, owned_directory,
    owned_tree, resource_plan, retire_previous_endpoint, source_closure, time_steps)
from mixed_material_dynamics import NativeMixedDynamics
from mixed_response_rk4 import advance_step
from mixed_response_storage import COMPONENTS, MappedMixedState
from substrate_framework.material_control_reference import bound_arrays
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.partitioned_q_storage import sync_directory
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-native-mixed-q-material-v1'


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--case', required=True)
    for label in ('initialization', 'history', 'execution-proof'):
        p.add_argument('--'+label, type=Path, required=True)
        p.add_argument('--'+label+'-sha256', required=True)
    p.add_argument('--frame-admission', nargs=2, action='append', required=True, metavar=('PATH', 'SHA256'))
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mode', choices=('plan', 'forward'), required=True)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--energy', action='store_true')
    p.add_argument('--mean-energy-ledger', type=Path)
    p.add_argument('--mean-energy-ledger-sha256')
    p.add_argument('--retain-endpoints', action='store_true')
    p.add_argument('--q', nargs=3, type=float, required=True)
    p.add_argument('--atom', type=int, required=True)
    p.add_argument('--axis', choices=('x', 'y', 'z'), required=True)
    p.add_argument('--amplitude', type=float, required=True)
    p.add_argument('--cell', nargs=3, type=int, required=True)
    p.add_argument('--time', type=float, required=True)
    p.add_argument('--subdivisions', type=int, required=True)
    p.add_argument('--cpus', nargs='+', type=int, choices=range(6, 11), required=True)
    p.add_argument('--max-actions', type=int, required=True)
    p.add_argument('--max-seconds', type=float, required=True)
    p.add_argument('--max-output-bytes', type=int, required=True)
    p.add_argument('--max-memory-bytes', type=int, required=True)
    p.add_argument('--worker-memory-bytes', type=int, required=True)
    p.add_argument('--block-size', type=int, default=8)
    return p


def caps(args):
    return {k: getattr(args, k) for k in ('max_actions', 'max_seconds', 'max_output_bytes',
                                         'max_memory_bytes', 'worker_memory_bytes')}


def require_arguments(args):
    if len(set(args.cpus)) != len(args.cpus) or not set(args.cpus) <= os.sched_getaffinity(0):
        raise ValueError('unique explicitly allocated available CPUs6-10 required')
    if (args.max_actions < 1 or not math.isfinite(args.max_seconds) or args.max_seconds <= 0
            or not 0 < args.worker_memory_bytes <= args.max_memory_bytes
            or args.max_output_bytes < 1 or args.block_size < 1 or args.subdivisions < 1):
        raise ValueError('positive explicit action/wall/live-byte/RSS/block/subdivision allocations required')
    if (args.max_memory_bytes > 18*1024**3 or args.worker_memory_bytes > 4*1024**3
            or args.max_output_bytes > 150*1024**3):
        raise ValueError('original aggregate18GiB/worker4GiB/live150GiB ceilings cannot be enlarged')
    if args.energy != (args.mean_energy_ledger is not None
                       and args.mean_energy_ledger_sha256 is not None):
        raise ValueError('--energy requires an explicit --mean-energy-ledger and digest')
    if not args.energy and (args.mean_energy_ledger is not None
                            or args.mean_energy_ledger_sha256 is not None):
        raise ValueError('mean energy ledger is only consumed with --energy')
    q = np.asarray(args.q)
    if not np.isfinite(q).all() or not np.allclose(4*q, np.rint(4*q), atol=1e-10, rtol=0):
        raise ValueError('original commensurate four-cube q required')
    if not math.isfinite(args.amplitude) or not math.isfinite(args.time) or args.time <= 0:
        raise ValueError('finite original local pulse amplitude and positive selected time required')
    output, directory = args.output.resolve(), args.directory.resolve()
    if not output.is_relative_to(directory) or output == directory:
        raise ValueError('--output must be an explicitly NEW JSON file within --directory')
    if output.exists() or output.is_symlink():
        raise FileExistsError('immutable output already exists; choose a new resume receipt: '+str(output))
    for parent in output.parents:
        if parent.is_symlink():
            raise ValueError('output cannot traverse a symlink')


@contextmanager
def watchdog(budget):
    """Bound the parent too, including initialization and native parent kernels."""
    parent = os.getpid()
    previous = signal.getsignal(signal.SIGALRM)
    def check(signum, frame):
        # Timers are not inherited by fork children. QWork separately supervises
        # its children, including forced termination and terminal cost receipts.
        if os.getpid() == parent:
            budget.check_processes(parent, [p.pid for p in mp.active_children() if p.is_alive()])
    signal.signal(signal.SIGALRM, check)
    signal.setitimer(signal.ITIMER_REAL, .2, .2)
    try:
        budget.check_processes(parent, [])
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def metadata_inputs(args):
    initialization, saved = selected_initialization(args)
    hb, history = checked_json(args.history, args.history_sha256)
    proof_binding, proof = checked_json(args.execution_proof, args.execution_proof_sha256)
    admit_oracle(proof, existing_dependencies())
    mean_energy = None
    if args.energy:
        mean_energy, _ = checked_json(args.mean_energy_ledger, args.mean_energy_ledger_sha256)
    normalization = dict(source=file_binding(HERE/'local_pulse_energy.py'),
        function='fourier_pulse; physical inverse-DFT coefficient = collective/sqrt(N)',
        convention='physical p_plus=J exp(-2 pi i q.cell)/N; p_minus=conj(p_plus)',
        derivative='independent analytic q/-q; eta=lambda*; mixed derivative has no half')
    if (history.get('schema') != 'streamed-joint-mean-v1' or history.get('complete') is not True
            or history.get('model_identity') != saved['model_identity'] or not history.get('history')):
        raise ValueError('completed actual original matched mean history required')
    first_binding = history['history'][0]['prior_boundary_binding']
    if first_binding != saved['joint_boundary']:
        raise ValueError('actual history does not start at the selected saved initialization')
    initial = json.loads(authenticate(first_binding).read_text())
    if (initial.get('time_ry_inverse') != 0 or initial.get('model_identity') != saved['model_identity']
            or initial.get('original_electronic_boundary') != saved['original_boundary']):
        raise ValueError('initial actual electronic/physical preparation authority differs')
    shapes, target_plus, target_minus = metadata_geometry(initial, args.q)
    steps = time_steps(history, args.time, args.subdivisions)
    small = initial_small(initial, q=args.q, atom=args.atom, axis=args.axis,
                          amplitude=args.amplitude, cell=args.cell)
    admissions, frame_metadata, frame_bindings = [], [], []
    for path, sha in args.frame_admission:
        binding, admission = checked_json(path, sha)
        if (admission.get('schema') != 'P254-real-frame-consumer-admission-v1'
                or admission.get('complete') is not True):
            raise ValueError('completed actual frame consumer admission required')
        cache = Path(admission['cache']).resolve()
        descriptor_binding = file_binding(cache/'cache.json')
        descriptor = json.loads(authenticate(descriptor_binding).read_text())
        complete_binding = file_binding(cache/'cache-complete.json')
        complete = json.loads(authenticate(complete_binding).read_text())
        original, expected = descriptor['original_request'], admission['expected_request']
        if (descriptor.get('owner') != admission['owner'] or descriptor.get('history') != hb
                or descriptor.get('cache') != str(cache) or admission['original_request'] != original
                or expected['input'] != original['input'] or expected['count'] != len(initial['slots'])
                or original['count'] != expected['count'] or complete.get('complete') is not True
                or complete.get('cache') != descriptor_binding):
            raise ValueError('explicit immutable original frame history/input/count/producer metadata differs')
        item_records = [json.loads(authenticate(b).read_text()) for b in complete['items']]
        # Original producer links and cache/adopted aliases are all retained;
        # none is assumed reclaimed because another namespace retired a link.
        present = [b for b in nested_bindings(
            admission['items'], admission['expected_request_file'],
            admission['original_request_binding'], admission['cache_complete'],
            admission['adoption'], descriptor, complete, item_records,
            descriptor_binding, complete_binding, binding)]
        frame_bindings.extend(present)
        frame_metadata.append(dict(proof=binding, cache_descriptor=descriptor_binding,
            cache_complete=complete_binding, original_request=original, admitted_request=expected))
        admissions.append((binding, admission))
    # Historical constructor source bindings in saved.specification are NOT
    # current physical payloads (the admitted saved-model adapter preserves
    # their original authority). Bind that specification through its JSON,
    # and inventory actual baseline payloads separately.
    baseline = [history['slots'], initial['slots'], saved['preparation'],
                saved['joint_boundary'], saved['original_boundary']]
    for record in history['history']:
        prior_binding = record['prior_boundary_binding']
        prior = json.loads(authenticate(prior_binding).read_text())
        baseline.extend((prior_binding, prior['slots'], record['arrays']))
    borrowed = nested_bindings(initialization, hb, proof_binding, normalization['source'],
                               baseline, frame_bindings, mean_energy)
    return dict(initialization=initialization, saved=saved, history_binding=hb, history=history,
        initial=initial, execution_binding=proof_binding, execution=proof, mean_energy=mean_energy,
        normalization=normalization, admissions=admissions, frame_metadata=frame_metadata,
        borrowed=borrowed, shapes=shapes, target_plus=target_plus, target_minus=target_minus,
        steps=steps, small=small, frame_requirements=frame_requirements(
            history, steps, block_size=args.block_size, admissions=admissions))


def validate_physical_preparation(case, histories, saved):
    history = histories['plus']
    if (file_binding(history.initial.path) != saved['joint_boundary']
            or history.initial.origin != saved['original_boundary']
            or case.identity != saved['model_identity'] or case.specification != saved['specification']):
        raise ValueError('saved initialization/history is not this exact original matched scalar')
    if (not np.array_equal(history.initial.u, np.zeros_like(case.reference.masses))
            or not np.array_equal(history.initial.p, np.zeros_like(case.reference.masses))
            or not np.array_equal(history.initial.coordinate, case.field_preparation['coordinate_mean'])
            or not np.array_equal(history.initial.momentum, case.field_preparation['momentum_mean'])):
        raise ValueError('actual initial canonical means differ from original product preparation')
    actual, expected = bound_arrays(saved['preparation']), preparation_arrays(case)
    if set(actual) != set(expected) or any(not np.array_equal(actual[k], v) for k, v in expected.items()):
        raise ValueError('actual saved ionic/field product preparation differs')


def execute(args, inputs, sources, budget, directory, result):
    case = bind_saved_model(load_case(args.case), inputs['saved'], inputs['execution'].get('source_compatibility'))
    original = case.make_history(inputs['history_binding'], budget, block_size=args.block_size)
    histories, geometries = {}, {}
    for component, q in (('plus', np.asarray(args.q)), ('minus', -np.asarray(args.q)), ('mixed', np.zeros(3))):
        key = tuple(q)
        if key not in geometries:
            # Mean/basis/D authority is q-independent and already authenticated.
            # bind_response replaces only the view's geometry and lazy actions.
            history = original if not geometries else copy(original)
            case.bind_response(history, q)
            geometries[key] = history
        histories[component] = geometries[key]
    validate_physical_preparation(case, histories, inputs['saved'])
    if args.atom == 37 and case.grid.reference.atom_names[36] != 'H':
        raise ValueError('original H37 physical atom authority differs')
    frames = RetainedFrameWork(inputs['admissions'], history_binding=inputs['history_binding'], budget=budget)
    shared_records = {}
    for history in histories.values():
        history.frame_work, history.frame_records = frames, shared_records
    dynamics = NativeMixedDynamics(case, histories, budget=budget, block_size=args.block_size)
    if args.energy:
        dynamics.energy_mean_readout = inputs['mean_energy']
    if (list(dynamics.target_plus) != inputs['target_plus']
            or list(dynamics.target_minus) != inputs['target_minus']):
        raise ValueError('metadata target plan differs from genuine full-PW physical geometry')
    for c in COMPONENTS:
        r = histories[c].response
        if (len(inputs['small'][c]['u']) != r.ndof
                or len(inputs['small'][c]['coordinate']) != r.nfield):
            raise ValueError('metadata physical per-cell/global-field dimensions differ from actual response')
    check_loaded_sources(sources)
    owned_frame_identity = digest(dict(owner=file_binding(directory/'owner.json'),
        dynamics=dynamics.source_identity, sources=sources, history=inputs['history_binding'],
        operation='original-qualified-missing-only-frozen-mean-exponential'))
    owned_frames = QWork(directory/'baseline-frame-work', budget, args.cpus, owned_frame_identity)
    frames.configure_owned(owned_frames, inputs['frame_requirements'], dynamics.count)
    result['owned_baseline_frames'] = dict(identity=owned_frame_identity,
        directory=str(owned_frames.directory), required=inputs['frame_requirements'],
        actual_requests=frames.execution_log,
        retention='only until last actual step/readout consumer commit; immutable costs/metadata survive')
    # Check planned inputs against ALL genuine active history stages before
    # authoring large seeds. This does not execute missing frame generation.
    for step in inputs['steps']:
        for history in histories.values():
            history.set_interval(step['start'], step['dt'])
        for node in (step['start'], step['start']+step['dt']/2, step['end']):
            if histories['plus']._saved_frames(node) is None:
                identity = histories['plus'].stage_identity(node)
                key = digest(identity)
                if inputs['frame_requirements'].get(key, {}).get('input') != canonical(identity):
                    raise ValueError('metadata missing-time plan differs from genuine supported history input')
                if key in frames.entries:
                    frames.execute('actual-baseline-frames', identity, dynamics.count, None)
    for history in histories.values():
        history.clear_interval()
    provenance = dict(model=case.identity, matched_case=args.case,
        initialization=inputs['initialization'], completed_mean=inputs['history_binding'],
        history_execution_proof=inputs['execution_binding'], physical_inputs=inputs['borrowed'],
        source_bindings=sources, normalization=inputs['normalization'],
        original_local_pulse=dict(atom=args.atom, axis=args.axis, amplitude=args.amplitude, cell=args.cell),
        electronic_initialization='new literal raw first/mixed zeros; no weighted-response import or D inverse',
        physical_initialization='first per-cell p q/-q kicks; all u/Maxwell/mixed small initially zero',
        fixed_raw_D_and_slots=dynamics.source_identity['slots'])
    initial = initialize_raw_zero(directory/'initial', shapes=inputs['shapes'], small=inputs['small'],
        q=args.q, target_plus=dynamics.target_plus, target_minus=dynamics.target_minus,
        identity=case.identity, provenance=provenance, budget=budget)
    scientific_identity = digest(dict(owner=file_binding(directory/'owner.json'),
        dynamics=dynamics.source_identity, initial=initial, dependencies=sources,
        actual_frame_reuse=frames.provenance))
    work = QWork(directory/'scientific-work', budget, args.cpus, scientific_identity)
    result.update(initial=initial, model_identity=case.identity,
        native_source_identity=dynamics.source_identity, qwork_identity=scientific_identity,
        frame_reuse=frames.provenance, steps_completed=[], final=initial)
    current, start_index = initial, 0
    owner_binding = file_binding(directory/'owner.json')
    # Resume follows compact durable endpoint commits, NOT replay from t0.
    # Old checkpoints remain readable as metadata after raw bases are retired.
    for index, step in enumerate(inputs['steps']):
        path = directory/'step-receipts'/f'{index:06d}.json'
        if not path.exists():
            if any((directory/'step-receipts'/f'{later:06d}.json').exists()
                   for later in range(index+1, len(inputs['steps']))):
                raise ValueError('noncontiguous durable endpoint commit chain')
            break
        commit_binding = file_binding(path)
        commit = json.loads(authenticate(commit_binding).read_text())
        if (commit.get('schema') != SCHEMA+'-step' or commit.get('complete') is not True
                or commit.get('owner') != owner_binding or commit.get('index') != index
                or commit.get('step') != step or commit.get('predecessor') != current
                or commit.get('qwork_identity') != scientific_identity):
            raise ValueError('authenticated resume endpoint/source/predecessor chain differs')
        endpoint = commit['final']
        metadata = json.loads(authenticate(endpoint).read_text())
        if (metadata.get('identity') != case.identity or metadata.get('time') != step['end']
                or metadata.get('shapes') != inputs['shapes'] or metadata.get('q') != args.q):
            raise ValueError('committed physical endpoint has changed raw geometry/time/model')
        expected_readouts = {'energy'} if args.energy and index == len(inputs['steps'])-1 else set()
        if set(commit.get('readouts', {})) != expected_readouts:
            raise ValueError('requested endpoint readout commit differs')
        for name, binding in commit['readouts'].items():
            if json.loads(authenticate(binding).read_text()).get('status') != 'complete':
                raise ValueError('resume endpoint readout is not complete')
            result[name+'_readout'] = binding
        previous, current = current, endpoint
        result['final'] = current
        summary = dict(index=index, interval=step['interval'], final=current, commit=commit_binding,
                       execution='reused authenticated actual endpoint; no RK4 replay')
        if previous != initial and not args.retain_endpoints:
            summary['previous_raw_retirement'] = retire_previous_endpoint(directory,
                previous=previous, successor=current, successor_receipt=commit_binding, budget=budget)
        if not args.retain_endpoints:
            summary['owned_frame_retirements'] = frames.retire_after_commit(index, commit_binding, shared_records)
        result['steps_completed'].append(summary)
        start_index = index+1
    # Authenticate only the latest actual endpoint payloads. Retired earlier
    # bases are never reopened and need no fake zero/reconstruction fallback.
    if start_index:
        with MappedMixedState(current, budget=budget) as state:
            for i in range(dynamics.count):
                for c in COMPONENTS:
                    state.leg(c, i)
                state.release_edge(i)
    for index in range(start_index, len(inputs['steps'])):
        step = inputs['steps'][index]
        budget.check_processes(os.getpid(), [])
        if json.loads(authenticate(current).read_text())['time'] != step['start']:
            raise ValueError('authenticated Native predecessor does not match the selected step start')
        for history in histories.values():
            history.set_interval(step['start'], step['dt'])
        previous = current
        current = advance_step(previous, step['dt'], dynamics=dynamics, work=work,
            directory=directory/'steps'/f'{index:06d}', identity=case.identity, budget=budget)
        with MappedMixedState(current, budget=budget) as state:
            if state.time != step['end'] or state.identity != case.identity:
                raise ValueError('Native producer did not return the actual selected physical endpoint')
        result['final'] = current
        readouts = {}
        if args.energy and index == len(inputs['steps'])-1:
            from mixed_material_energy import mixed_energy_readout
            readout = mixed_energy_readout(current, dynamics=dynamics, work=work,
                directory=directory/'energy', budget=budget)
            report = json.loads(authenticate(readout).read_text())
            result['energy_readout'] = readout
            if report.get('status') != 'complete':
                result.update(status='partial', evolution_complete=True,
                    reason='selected physical energy readout unavailable; predecessor retained',
                    energy_status=report.get('status'))
                return
            readouts['energy'] = readout
        commit_binding = exact_publish(directory/'step-receipts'/f'{index:06d}.json',
            dict(schema=SCHEMA+'-step', complete=True, owner=owner_binding, index=index,
                step=step, predecessor=previous, final=current, readouts=readouts,
                qwork_identity=scientific_identity, observed_budget=budget.snapshot(),
                completed_at=datetime.now(timezone.utc).isoformat()), budget)
        summary = dict(index=index, interval=step['interval'], final=current, commit=commit_binding,
                       execution='actual Native RK4; original per-edge QWork costs retained')
        if previous != initial and not args.retain_endpoints:
            summary['previous_raw_retirement'] = retire_previous_endpoint(directory,
                previous=previous, successor=current, successor_receipt=commit_binding, budget=budget)
        if not args.retain_endpoints:
            summary['owned_frame_retirements'] = frames.retire_after_commit(index, commit_binding, shared_records)
        result['steps_completed'].append(summary)
    result['evolution_complete'] = True
    check_loaded_sources(sources)
    result['status'] = 'complete'


def main(argv=None):
    args = parser().parse_args(argv)
    require_arguments(args)
    os.sched_setaffinity(0, {args.cpus[0]})
    started, cpu_started = perf_counter(), process_time()
    invocation_started_unix = time()
    result = dict(schema=SCHEMA, status='error', mode=args.mode, q=args.q, case=args.case,
        selected_time=args.time, subdivisions=args.subdivisions, cpus=args.cpus, blas_threads=1,
        whole_spatial_reconstruction=False, required_q_for_whole_local=64,
        numerical_accuracy_certified=False, production_authorization=False,
        scope='one selected-q raw coherent mixed derivative; not noise, finite-dose energy, heat or full E/D/C')
    budget = None
    code = 2
    with owned_directory(args.directory, args.resume) as directory:
        args.directory = directory
        args.output = args.output.resolve()
        invocation = uuid4().hex
        receipt = directory/'invocations'/invocation/'receipt.json'
        receipt.parent.mkdir(parents=True)
        try:
            inputs = metadata_inputs(args)
            source_seeds = existing_dependencies()+[file_binding(HERE/name) for name in (
                'run_mixed_q_material.py', 'mixed_material_driver_support.py', 'mixed_material_dynamics.py',
                'mixed_response_rk4.py', 'mixed_response_edge_transition.py', 'mixed_response_storage.py',
                'mixed_payload_output.py', 'local_pulse_energy.py',
                'mixed_material_edge.py', 'mixed_material_moments.py', 'mixed_coordinate_response.py',
                'mixed_frame_response.py', 'mixed_physical_field_dyads.py', 'actual_frame_cache_io.py',
                'actual_frame_cache_provenance.py', 'reuse_actual_baseline_frames.py')]
            if args.energy:
                source_seeds.append(file_binding(HERE/'mixed_material_energy.py'))
            sources = source_closure(source_seeds)
            request = canonical(dict(schema=SCHEMA+'-owner', case=args.case,
                initialization=inputs['initialization'], history=inputs['history_binding'],
                execution_proof=inputs['execution_binding'], normalization=inputs['normalization'],
                frames=inputs['frame_metadata'], q=args.q, atom=args.atom, axis=args.axis,
                frame_requirements=inputs['frame_requirements'],
                amplitude=args.amplitude, cell=args.cell, time=args.time, subdivisions=args.subdivisions,
                cpus=args.cpus, caps=caps(args), block_size=args.block_size, energy=args.energy,
                retain_endpoints=args.retain_endpoints, mean_energy_ledger=inputs['mean_energy'],
                source_bindings=sources, inputs=inputs['borrowed']))
            owner_path = directory/'owner.json'
            if owner_path.exists():
                owner = json.loads(owner_path.read_text())
                if owner.get('request') != request:
                    raise ValueError('resume owner/source/input/selection/caps differ; no fresh continuation admitted')
            else:
                if args.resume:
                    raise ValueError('missing original owner request')
                owner = dict(request=request, started_unix=invocation_started_unix,
                    started_at=datetime.now(timezone.utc).isoformat())
                # Owner metadata precedes the budget so its original absolute
                # clock is never reset, including a crash before first action.
                exact_publish(owner_path, owner, None)
            # Accepted input authority is immutable across producer retirement.
            # Existing readonly aliases consume disk, not this run's output cap.
            inputs['borrowed'] = nested_bindings(inputs['borrowed'], sources)
            existing = owned_tree(directory)
            result.update(owner=file_binding(owner_path), source_bindings=sources,
                initialization=inputs['initialization'], history=inputs['history_binding'],
                history_execution_proof=inputs['execution_binding'], normalization=inputs['normalization'],
                borrowed_inputs=inputs['borrowed'])
            plan = resource_plan(inputs['shapes'], inputs['steps'], inputs['borrowed'],
                existing=existing, workers=len(args.cpus), caps=caps(args), history=inputs['history'],
                gamma=bool(np.all(np.asarray(args.q) == 0)), retain_endpoints=args.retain_endpoints,
                frame_requirements=inputs['frame_requirements'])
            result['resource_plan'] = plan
            budget = ContinuationBudget(directory, started_at=owner['started_unix'],
                max_actions=args.max_actions, max_seconds=args.max_seconds,
                max_bytes=args.max_output_bytes, max_memory_bytes=args.max_memory_bytes,
                worker_memory_bytes=args.worker_memory_bytes,
                existing_bytes=sum(existing.values()))
            # Only owned output consumes the live-output cap. The plan reports
            # readonly aliases and their physical allocation separately.
            terminal_reservation = 2*(4*len(json.dumps(request))+len(json.dumps(plan))
                                      +8192*len(inputs['steps'])+1024**2)
            budget.reserve(max(0, args.max_actions*budget.RECORD.size-budget.ledger.stat().st_size)
                           +terminal_reservation)
            with watchdog(budget):
                if args.mode == 'plan':
                    result.update(status='planned', physics_executed=False)
                else:
                    execute(args, inputs, sources, budget, directory, result)
            code = 0 if result['status'] in ('complete', 'planned') else 2
        except BaseException as error:
            result.update(status='partial' if result.get('final') else 'error',
                          error=type(error).__name__+': '+str(error))
        finally:
            result.update(invocation=invocation, invocation_wall_seconds=perf_counter()-started,
                invocation_parent_cpu_seconds=process_time()-cpu_started,
                parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                budget=budget.snapshot() if budget is not None else None,
                terminal_at=datetime.now(timezone.utc).isoformat())
            # Terminal receipts must survive a cap/deadline failure. Their space
            # was reserved before execution; never invoke an exhausted budget.
            receipt_binding = exact_publish(receipt, result, None)
            exact_publish(args.output, dict(result, invocation_receipt=receipt_binding), None)
            sync_directory(args.output.parent)
    print(json.dumps(dict(schema=SCHEMA, status=result['status'], output=str(args.output),
                          final=result.get('final'), error=result.get('error')), sort_keys=True))
    return code


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())
