#!/usr/bin/env python3
"""Bounded real-sector workload pilot, NOT a full-material evolution result."""
import argparse
import os
from pathlib import Path
import resource
from time import perf_counter

from run_streamed_material_evolution import (
    ROOT, HERE, np, explicit_binding, bound_json, preload_parsers, restore_inputs,
    build_evolution_model, slot_identity, publish, implementation, authenticate,
    DiskMean, physical_identity)
from substrate_framework.streamed_material_evolution import exponential_action
from substrate_framework.streamed_sector_executor import (
    SectorExecutor, physical_specification, KNOWN_T0_SHA256, KNOWN_T0_IDENTITY)
from threadpoolctl import threadpool_limits, threadpool_info


def memory():
    fields = {}
    for path in ('/proc/self/status', '/proc/self/smaps_rollup'):
        for line in Path(path).read_text().splitlines():
            name, _, value = line.partition(':')
            if name in ('VmRSS', 'VmHWM', 'VmSize', 'Rss', 'Pss', 'Private_Dirty', 'Shared_Clean'):
                fields[name+'_kib'] = int(value.split()[0])
    fields['ru_maxrss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('aggregate-receipt', 'manifest', 'checkpoint', 'initial-ions', 'oracle'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sector', required=True, help='Exact spin:mesh identity; no implicit selected Gamma')
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--physical-interval-ry-inverse', type=float, required=True)
    parser.add_argument('--max-h0-actions', type=int, required=True)
    parser.add_argument('--phase', choices=('action', 'moments', 'exponential', 'all'), default='all')
    parser.add_argument('--block-size', type=int, default=8)
    parser.add_argument('--reality-tolerance', type=float, default=1e-10)
    args = parser.parse_args()
    executor = SectorExecutor((args.cpu,))
    if (args.max_h0_actions < 1 or not np.isfinite(args.physical_interval_ry_inverse)
            or args.physical_interval_ry_inverse <= 0):
        parser.error('explicit positive action cap and physical interval required')
    os.sched_setaffinity(0, {args.cpu})
    bindings = {name: explicit_binding(getattr(args, name), getattr(args, name+'_sha256'))
                for name in ('aggregate_receipt', 'manifest', 'checkpoint', 'initial_ions', 'oracle')}
    oracle = bound_json(bindings['oracle'])
    if oracle.get('status') != 'passed' or oracle.get('schema') != 'P254-streamed-evolution-oracle-v1':
        raise ValueError('authenticated streamed equation/concurrency oracle required')
    for entry in oracle['source_bindings']:
        authenticate(entry)
    required = {str((ROOT/'src/substrate_framework'/name).resolve()) for name in
                ('streamed_material_state.py', 'streamed_material_evolution.py', 'streamed_sector_executor.py')}
    if not required <= {entry['path'] for entry in oracle['source_bindings']}:
        raise ValueError('oracle does not bind the actual sector consumers')
    args.output.mkdir(parents=True, exist_ok=False)
    sources = implementation((Path(__file__), HERE/'run_streamed_material_evolution.py',
                              HERE/'streamed-sector-registration.json'))
    preload_parsers()
    arrays, receipt, manifest, saved, original, grid, ordered, pairs = restore_inputs(
        bindings['aggregate_receipt'], bindings['manifest'])
    ions = bound_json(bindings['initial_ions'])
    spec = dict(aggregate=bindings['aggregate_receipt'], manifest=bindings['manifest'],
        preparation=manifest['preparation_terminal'], initial_ions=bindings['initial_ions'],
        scalar='PreparedReference joint-second-order; common Gamma coefficient only; actual per-slot H0')
    # Authenticate the exact boundary once, then only this pilot's frame and D.
    # No claim to have rehashed or evolved the other127 physical sectors.
    state = DiskMean(bindings['checkpoint']['path'], verify=False)
    known = (bindings['checkpoint']['sha256'] == KNOWN_T0_SHA256 and state.identity == KNOWN_T0_IDENTITY)
    if (not known and state.identity != physical_identity(spec)):
        raise ValueError('pilot requires the known t0 or the identical physical model boundary')
    if (state.time != 0. or not state.evidence.get('original_raw_frames_unchanged')
            or physical_specification(state.evidence['model']) != spec
            or not np.array_equal(state.u, ions['u_bohr'])
            or not np.array_equal(state.p, ions['p_canonical'])):
        raise ValueError('pilot requires the unchanged physical initial state')
    selected = [slot for slot in state.slots if slot.identity == args.sector]
    if len(selected) != 1:
        raise ValueError('explicit sector absent from complete boundary')
    slot = selected[0]
    authenticate(slot.frame)
    authenticate(slot.coefficient)
    model, background_receipt = build_evolution_model(args, arrays, manifest, saved,
        original, grid, ordered, bindings['aggregate_receipt'], ions, state.identity)
    model.validate(state)
    # Frozen aggregate mean, not a recomputed global nonlinear midpoint.
    stage = model.stage(arrays['conditional_density_bohr3'],
                        arrays['positive_gradient_ry_bohr'].reshape(-1), state.u)

    def work(_):
        observations = {}
        def timed(name, function):
            before, started = model.counts, perf_counter()
            publish(args.output/('started-'+name+'.json'),
                    dict(phase=name, counts=before, worker_pid=os.getpid()))
            try:
                return function()
            finally:
                after = model.counts
                observations[name] = dict(seconds=perf_counter()-started,
                    work_counts={key: after[key]-before[key] for key in after}, memory=memory())
                publish(args.output/('measured-'+name+'.json'), observations[name])
        started = perf_counter()
        psi, d = state.arrays(slot)
        observations['frame_mmap'] = dict(seconds=perf_counter()-started, memory=memory())
        publish(args.output/'measured-frame-mmap.json', observations['frame_mmap'])
        started = perf_counter()
        publish(args.output/'started-operator-construction.json', dict(worker_pid=os.getpid()))
        action_norm = moment_summary = exponential_norm = None
        status, failure = 'completed', None
        with model.operators(slot) as (coupling, h0):
            observations['operator_construction'] = dict(seconds=perf_counter()-started, memory=memory())
            publish(args.output/'measured-operator-construction.json', observations['operator_construction'])
            action = lambda value: model.action(coupling, h0, slot, stage, value)
            if args.phase in ('action', 'all'):
                applied = timed('full_frame_action', lambda: action(psi))
                action_norm = float(np.linalg.norm(applied))
                del applied
            if args.phase in ('moments', 'all'):
                def moments():
                    result = coupling.prepared_moments(psi, d)
                    model._counts.add('moment_reductions')
                    return result
                density, gradient = timed('full_frame_moments', moments)
                moment_summary = dict(density_norm=float(np.linalg.norm(density)),
                                      positive_gradient_norm=float(np.linalg.norm(gradient)))
                del density, gradient
            if args.phase in ('exponential', 'all'):
                try:
                    evolved = timed('full_frame_exponential', lambda: exponential_action(
                        action, psi, args.physical_interval_ry_inverse))
                    exponential_norm = float(np.linalg.norm(evolved))
                    del evolved
                except RuntimeError as error:
                    if 'H0 action cap exhausted' not in str(error):
                        raise
                    status, failure = 'workload-cap-exhausted', str(error)
        return dict(status=status, failure=failure, measurements=observations, phase=args.phase,
            full_frame_action_norm=action_norm, moments=moment_summary,
            exponential_frame_norm=exponential_norm, worker_pid=os.getpid(),
            worker_affinity=sorted(os.sched_getaffinity(0)), native_threadpools=threadpool_info())

    result, = list(executor.map(work, (None,)))
    for entry in sources:
        authenticate(entry)
    result.update(schema='P254-streamed-sector-pilot-v1',
        scope='one real full-raw-frame workload/method pilot; NOT full material evolution',
        sector=slot.metadata(), frame=slot.frame, original_D=slot.coefficient,
        inputs=bindings, implementation_sources=sources, background_receipt=background_receipt,
        physical_interval_ry_inverse=args.physical_interval_ry_inverse,
        max_h0_actions=args.max_h0_actions, work_counts=model.counts,
        potential='real saved-density kernel of aggregate conditional density minus saved density',
        limitations=['No global128-sector reduction or nonlinear midpoint executed',
            'Aggregate gradient retains declared geometric approximation',
            'No output trajectory, coherence, equilibrium or material rate result',
            'Observed single-worker RSS/PSS and timings do not certify full-stage runtime or RAM',
            'Only selected frame/D rehashed; exact whole boundary authenticated'])
    publish(args.output/'result.json', result)

if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
