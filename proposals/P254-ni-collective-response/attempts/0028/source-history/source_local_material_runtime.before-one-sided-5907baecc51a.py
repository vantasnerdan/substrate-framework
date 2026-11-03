"""Actual saved-material construction and one-edge measurement for source-local q.

Metadata restoration authenticates receipts, not a newly evolved trajectory. Only
``build_runtime`` constructs the original physical providers. Two-sided
covariance and bounded Native measurements are separate consumers.
"""
from dataclasses import dataclass
from functools import partial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]

# This is the original physical provider's import/BLAS policy, not a second
# implementation of its model or executable entry point.
import run_joint_q_material as physical
import numpy as np

from actual_frame_cache_provenance import _complete_cache, _history_slots, _request
from mixed_material_driver_support import (
    RetainedFrameWork, checked_json, frame_requirements, nested_bindings, time_steps)
from substrate_framework.material_control_reference import bound_arrays, bound_json
from substrate_framework.partitioned_q_response import PartitionedQResponse, frame_times
from substrate_framework.source_local_q_plan import source_local_initial_plan, source_local_plan
from substrate_framework.streamed_material_state import authenticate, file_binding


SOURCE_LOCAL_MODULES = (
    'source_local_q_state.py', 'source_local_q_context.py',
    'source_local_q_response.py', 'source_local_q_covariance.py',
    'source_local_q_plan.py', 'source_local_q_edge.py', 'source_local_q_contraction.py')


@dataclass(frozen=True)
class MaterialRuntime:
    executor: PartitionedQResponse
    preparation: physical.JointProductPreparation
    probe: object
    times: np.ndarray
    second_times: np.ndarray
    provenance: dict
    resource_plan: dict


def material_sources():
    """Current physical-provider and source-local seed bindings for the owner."""
    paths = [HERE/'source_local_material_runtime.py', HERE/'run_source_local_q_oracle.py',
             HERE/'source_local_material_probe.py',
             HERE/'mixed_material_driver_support.py', HERE/'actual_frame_cache_provenance.py']
    paths += [ROOT/'src/substrate_framework'/name for name in SOURCE_LOCAL_MODULES]
    sources = physical.bindings()+[file_binding(path) for path in paths]
    return list({source['path']: source for source in sources}.values())


def _selectors(args):
    q = np.asarray(args.q, float)
    if q.shape != (3,) or not np.isfinite(q).all():
        raise ValueError('three finite original-mesh q coordinates required')
    physical.validate_selection(args.observable, args.field_axis, q)
    if (type(args.subdivisions) is not int or args.subdivisions < 1
            or type(args.block_size) is not int or args.block_size < 1
            or type(args.atom) is not int or args.atom < 1 or args.axis not in ('x', 'y', 'z')
            or args.dynamics not in ('screened-joint', 'bare-electronic')):
        raise ValueError('original atom/axis and positive subdivision/block sizes required')
    cpus = list(args.cpus)
    if (not cpus or any(type(cpu) is not int or cpu < 0 for cpu in cpus)
            or len(set(cpus)) != len(cpus)):
        raise ValueError('explicit unique allocated CPUs required')
    return dict(case=args.case, q=q.tolist(), time=float(args.time),
        second_time=float(args.second_time), subdivisions=args.subdivisions,
        atom=args.atom, axis=args.axis, observable=args.observable, field_axis=args.field_axis,
        dynamics=args.dynamics, block_size=args.block_size)


def _metadata_grid(history, final_time, subdivisions):
    if not np.isfinite(final_time) or not 0 <= final_time <= history['time_ry_inverse']:
        raise ValueError('finite selected time within completed original history required')
    if final_time == 0:
        return np.asarray([0.])
    steps = time_steps(history, final_time, subdivisions)
    return np.asarray([0.]+[step['end'] for step in steps])


def _binding_equal(left, right):
    return (left['path'] == right['path'] and left['sha256'] == right['sha256']
            and int(left['bytes']) == int(right['bytes']))


def _admissions(args, history_binding, slots):
    admissions, metadata, borrowed = [], [], []
    for path, sha in getattr(args, 'frame_admission', None) or ():
        binding, proof = checked_json(path, sha)
        if (proof.get('schema') != 'P254-real-frame-consumer-admission-v1'
                or proof.get('complete') is not True):
            raise ValueError('completed original real-frame consumer admission required')
        for source in proof['implementation_sources']:
            authenticate(source)
        descriptor, items = _complete_cache(Path(proof['cache']), proof['owner'])
        original, expected = descriptor['original_request'], proof['expected_request']
        _request(expected, expected['source_identity'])
        if (not _binding_equal(descriptor['history'], history_binding)
                or descriptor['slots'] != slots
                or proof['original_request'] != original
                or expected['input'] != original['input']
                or expected['count'] != len(slots) or original['count'] != len(slots)
                or expected['input']['block_size'] != args.block_size):
            raise ValueError('admitted cache differs from exact saved history/input/full geometry/block size')
        cache = Path(proof['cache'])
        descriptor_binding = file_binding(cache/'cache.json')
        complete_binding = file_binding(cache/'cache-complete.json')
        if not _binding_equal(proof['cache_complete'], complete_binding):
            raise ValueError('frame admission complete barrier differs from original live cache')
        admissions.append((binding, proof))
        metadata.append(dict(proof=binding, cache_descriptor=descriptor_binding,
            cache_complete=complete_binding, original_request=original, admitted_request=expected))
        borrowed.extend(nested_bindings(binding, proof['items'], proof['expected_request_file'],
            proof['original_request_binding'], proof['adoption'], descriptor_binding,
            complete_binding, items))
    return admissions, metadata, borrowed


def material_inputs(args):
    """Authenticate actual saved metadata and existing proofs without Native work.

    Cache admissions use the original read-only provenance validator. It may
    authenticate full cached payload bytes, but never reconstructs a frame.
    Original baseline PW payloads are authenticated by the original physical
    history restoration in build_runtime, not misrepresented as loaded here.
    """
    selectors = _selectors(args)
    initialization, saved = physical.selected_initialization(args)
    hb, history = checked_json(args.history, args.history_sha256)
    qualification_binding, qualification = checked_json(
        args.physical_qualification, args.physical_qualification_sha256)
    physical.admit_oracle(qualification, physical.bindings(), args.observable)
    if any(check.get('passed') is not True for check in qualification.get('checks', ())):
        raise ValueError('passed current-source physical qualification required')
    oracle_binding, oracle = checked_json(args.source_local_oracle, args.source_local_oracle_sha256)
    if (oracle.get('schema') != 'P254-source-local-q-oracle-v1' or oracle.get('status') != 'passed'
            or not oracle.get('checks')
            or any(check.get('passed') is not True for check in oracle['checks'])):
        raise ValueError('executed current-source numerical proof required')
    oracle_sources = {source['path']: source for source in oracle['sources']}
    for source in oracle['sources']:
        authenticate(source)
    for name in SOURCE_LOCAL_MODULES:
        current = file_binding(ROOT/'src/substrate_framework'/name)
        if current != oracle_sources.get(current['path']):
            raise ValueError('source-local proof does not bind the current implementation: '+name)
    if (history.get('schema') != 'streamed-joint-mean-v1' or history.get('complete') is not True
            or history.get('model_identity') != saved['model_identity'] or not history.get('history')):
        raise ValueError('completed original matched actual moving history required')
    record = history['history'][0]
    if not _binding_equal(record['prior_boundary_binding'], saved['joint_boundary']):
        raise ValueError('actual history must begin at the exact saved initialization')
    initial = bound_json(saved['joint_boundary'])
    if (initial.get('time_ry_inverse') != 0 or initial.get('model_identity') != saved['model_identity']
            or initial.get('original_electronic_boundary') != saved['original_boundary']):
        raise ValueError('saved original t0 physical/electronic preparation authority differs')
    # The original cache provenance API validates the entire prior/final chain,
    # original basis/gauge/D and actual midpoint inputs, without loading PW data.
    identity = dict(time=record['time_ry_inverse'], interval=0, stage=record,
        baseline_boundary=record['prior_boundary_binding'], model=saved['model_identity'],
        origin=saved['original_boundary'], block_size=args.block_size)
    slots = _history_slots(identity, 128, hb)
    indices = [tuple(int(i) for i in slot['identity'].split(':')[1].split(','))
               for slot in initial['slots']]
    mesh = [1+max(index[axis] for index in indices) for axis in range(3)]
    if (len(initial['slots']) != 128 or mesh != [4, 4, 4]
            or {slot['identity'].split(':')[0] for slot in initial['slots']} != {'1', '2'}):
        raise ValueError('literal full 128 spin/mesh sectors of the original material required')
    if args.atom > len(initial['u_bohr'])//3:
        raise ValueError('probe atom is outside the original saved material')
    background_binding = saved['specification']['context']['background_receipt']
    background = bound_json(background_binding)
    density_shape = background['arrays']['rho_updown']['shape']
    times = _metadata_grid(history, args.time, args.subdivisions)
    second_times = _metadata_grid(history, args.second_time, args.subdivisions)
    # This also validates complete Cartesian/spin coverage and commensurate q.
    forecast = source_local_initial_plan(saved['joint_boundary'], q=args.q, mesh=mesh,
        density_shape=density_shape, times=times, second_times=second_times,
        workers=len(args.cpus), dynamics=args.dynamics)
    if forecast['model_identity'] != saved['model_identity']:
        raise ValueError('metadata geometry is not the selected original saved model')
    admissions, frame_metadata, frame_borrowed = _admissions(args, hb, slots)
    baseline = [saved['joint_boundary'], saved['original_boundary'], saved['preparation'],
                initial['slots'], history['slots']]
    for stage in history['history']:
        prior = bound_json(stage['prior_boundary_binding'])
        baseline.extend((stage['prior_boundary_binding'], stage['arrays'], prior['slots']))
    borrowed = nested_bindings(initialization, hb, qualification_binding, oracle_binding,
        background_binding, baseline, frame_borrowed)
    return dict(initialization=initialization, saved=saved, history_binding=hb, history=history,
        initial=initial, physical_qualification_binding=qualification_binding,
        physical_qualification=qualification, source_local_oracle_binding=oracle_binding,
        source_local_oracle=oracle, admissions=admissions, frame_metadata=frame_metadata,
        borrowed=borrowed, mesh=mesh, density_shape=density_shape,
        density_shape_metadata=background_binding, times=times, second_times=second_times,
        selectors=selectors, numerical_work_executed=False,
        baseline_payload_authentication='deferred to original physical history restoration')


def metadata_plan(args, inputs):
    """Forecast only: saved receipts do not turn proposed work into execution."""
    if inputs['selectors'] != _selectors(args):
        raise ValueError('metadata selectors changed after input admission')
    plan = source_local_initial_plan(inputs['saved']['joint_boundary'], q=args.q,
        mesh=inputs['mesh'], density_shape=inputs['density_shape'], times=inputs['times'],
        second_times=inputs['second_times'], workers=len(args.cpus), dynamics=args.dynamics)
    plan.update(status='metadata-only-forecast-no-physical-construction-or-propagation',
        saved_history_metadata=inputs['history_binding'], actual_moving_history_constructed=False,
        positive_time_execution_evidence=False, density_shape_metadata=inputs['density_shape_metadata'],
        authenticated_frame_admissions=inputs['frame_metadata'],
        frame_forecast_policy='Conservative original initial-boundary plan; no positive-time boundary or cache aliases are inferred.',
        scientific_scope='Proposed selected-q two-sided source-local work, not material noise or a production allocation.')
    return plan


def build_runtime(args, inputs, *, budget, directory, source_bindings):
    """Restore the exact original model/preparation and actual moving history."""
    if inputs['selectors'] != _selectors(args):
        raise ValueError('physical selectors changed after input admission')
    for key in ('initialization', 'history_binding', 'physical_qualification_binding',
                'source_local_oracle_binding'):
        authenticate(inputs[key])
    sources = {source['path']: source for source in source_bindings}
    if len(sources) != len(source_bindings):
        raise ValueError('duplicate runtime source bindings')
    for source in material_sources():
        current = sources.get(source['path'])
        if current is None or not _binding_equal(source, current):
            raise ValueError('runtime owner does not bind the current physical/source-local providers')
    for source in source_bindings:
        budget.check()
        authenticate(source)
    saved = inputs['saved']
    context = physical.bind_saved_model(physical.load_case(args.case), saved,
        inputs['physical_qualification'].get('source_compatibility'))
    if context.identity != saved['model_identity'] or context.specification != saved['specification']:
        raise ValueError('restored physical scalar/preparation differs from selected saved model')
    history = context.make_history(inputs['history_binding'], budget, block_size=args.block_size)
    if (not _binding_equal(file_binding(history.initial.path), saved['joint_boundary'])
            or history.initial.origin != saved['original_boundary']):
        raise ValueError('moving history differs from the selected literal original preparation')
    if (not np.array_equal(history.initial.u, np.zeros_like(context.reference.masses))
            or not np.array_equal(history.initial.p, np.zeros_like(context.reference.masses))
            or not np.array_equal(history.initial.coordinate, context.field_preparation['coordinate_mean'])
            or not np.array_equal(history.initial.momentum, context.field_preparation['momentum_mean'])):
        raise ValueError('initial canonical means differ from selected matched product preparation')
    prepared, expected = bound_arrays(saved['preparation']), physical.preparation_arrays(context)
    if set(prepared) != set(expected) or any(not np.array_equal(prepared[k], v) for k, v in expected.items()):
        raise ValueError('saved ionic/field product preparation differs from the original physical model')
    response = context.bind_response(history, np.asarray(args.q, float))
    if (len(response.sectors) != 128 or len(response.edges) != 128
            or any(tuple(sector.shape) != tuple(slot['shape'])
                   for sector, slot in zip(response.sectors, inputs['initial']['slots']))):
        raise ValueError('original full-PW rows/columns or full128 response geometry changed')
    times = history.time_grid(args.time, args.subdivisions)
    second_times = history.time_grid(args.second_time, args.subdivisions)
    if (not np.array_equal(times, inputs['times'])
            or not np.array_equal(second_times, inputs['second_times'])):
        raise ValueError('actual restored time grids differ from admitted saved history metadata')
    if not 1 <= args.atom <= context.grid.reference.nat:
        raise ValueError('probe requires an original saved atom')
    if args.atom == 37 and context.grid.reference.atom_names[36] != 'H':
        raise ValueError('original H37 authority differs')
    direction = np.zeros(len(context.reference.masses))
    direction[3*(args.atom-1)+'xyz'.index(args.axis)] = 1
    observable_metadata = None
    if args.observable == 'total-force':
        probe = physical.total_force_probe(response, direction, identity='original-%s%d%s-total-force' % (
            context.grid.reference.atom_names[args.atom-1], args.atom, args.axis))
    else:
        probe, observable_metadata = physical.global_field_probe(response, context.field,
            observable=args.observable, axis=args.field_axis)
    probe.validate(response)
    preparation = physical.JointProductPreparation(response, history.initial_stage(),
        ionic_ordered=prepared['local_ordered'], ionic_minus_ordered=prepared['local_ordered'],
        field_ordered=prepared['field_ordered_covariance'] if response.gamma else np.zeros((0, 0)),
        electronic_identity=saved['original_boundary']['sha256'],
        ionic_identity=saved['preparation']['sha256'], field_identity=saved['preparation']['sha256'],
        product_identity='original-matched-product-'+context.identity)
    request = dict(initialization=inputs['initialization'], baseline=inputs['history_binding'],
        physical_qualification=inputs['physical_qualification_binding'],
        source_local_oracle=inputs['source_local_oracle_binding'], selectors=inputs['selectors'],
        times=times.tolist(), second_times=second_times.tolist(),
        frame_admissions=[binding for binding, _ in inputs['admissions']],
        observable_metadata=observable_metadata)
    identity = physical.producer_identity(context.identity, list(args.q), args.dynamics,
        source_bindings, request)
    directory = Path(directory).resolve()
    measure = abs(np.linalg.det(context.grid.reference.cell_bohr))/int(np.prod(context.grid.reference.fft_smooth))
    executor = PartitionedQResponse(response, history,
        kernel_adjoint=partial(physical.physical_kernel_adjoint, context.background),
        grid_measure=measure, budget=budget, directory=directory/'executor', identity=identity,
        dynamics=args.dynamics, cpus=args.cpus, work_directory=directory/'partitioned-work')
    frame_reuse = []
    if inputs['admissions']:
        retained = RetainedFrameWork(inputs['admissions'], history_binding=inputs['history_binding'], budget=budget)
        steps = []
        for end in (args.time, args.second_time):
            if end:
                steps.extend(time_steps(inputs['history'], end, args.subdivisions))
        requirements = frame_requirements(inputs['history'], steps,
            block_size=args.block_size, admissions=inputs['admissions'])
        retained.configure_owned(executor.frames, requirements, len(response.sectors))
        # Keep executor.frames as the genuine owned QWork for original retirement
        # APIs; only the history's read-only frame dispatcher gains admissions.
        history.frame_work = retained
        frame_reuse = retained.provenance
    provenance = dict(case=args.case, model_identity=context.identity, baseline=inputs['history_binding'],
        initialization=inputs['initialization'], history=history.records,
        physical_qualification=inputs['physical_qualification_binding'],
        source_local_numerical_proof=inputs['source_local_oracle_binding'], source_bindings=source_bindings,
        request=request, matched_aggregate=context.aggregate_binding,
        scalar=context.specification['scalar'], observable_metadata=observable_metadata,
        observable=probe.identity, all_explicit_edges=len(response.edges), fixed_original_D=True,
        frame_reuse=frame_reuse,
        mean_history_method='authenticated actual initial/final boundaries and saved midpoint stages; original full-PW exponential for missing interior frames',
        whole_spatial_reconstruction=False, numerical_accuracy_certified=False,
        scientific_scope='Actual saved selected-q runtime; no pair/noise execution claimed by construction.')
    executor.provenance = provenance
    executor.readout_source_bindings = source_bindings
    executor.readout_operator_binding = dict(request=request, model=context.identity,
        observable=args.observable, field_axis=args.field_axis, probe_identity=probe.identity)
    aliases = [state.time for state in history.priors]+[history.final.time]
    plan = source_local_plan(response, times, second_times, len(args.cpus),
        aliased_frame_times=aliases, dynamics=args.dynamics)
    needed = frame_times(times) | frame_times(second_times)
    reused_times = sorted({float(proof['expected_request']['input']['time'])
                          for _, proof in inputs['admissions']} & needed)
    plan.update(actual_moving_history_constructed=True, baseline=inputs['history_binding'],
        original_sector_count=128, original_edge_count=128,
        authenticated_retained_frame_times=reused_times,
        frame_reuse_accounting='Conservative payload bound counts retained cache times as if newly owned; only actual boundary times are labeled boundary aliases.',
        numerical_work_executed=False, production_allocation_granted=False)
    return MaterialRuntime(executor, preparation, probe, times, second_times, provenance, plan)


