#!/usr/bin/env python3
"""Bounded actual one-sector full-frame frozen-H block Krylov pilot.

Run only after the independent oracle passes. Large retained arrays must live
outside Git. An optional settled SciPy receipt is consumed only when explicit.
"""
import argparse
import json
import os
from pathlib import Path
import resource
from time import perf_counter
from types import SimpleNamespace as NS

from run_block_krylov_oracle import (
    HERE, ROOT, np, binding, authenticate, explicit, publish, sources)
from run_streamed_material_evolution import (
    bound_json, preload_parsers, restore_inputs, build_evolution_model,
    DiskMean, physical_identity, implementation)
from substrate_framework.streamed_sector_executor import (
    physical_specification, KNOWN_T0_SHA256, KNOWN_T0_IDENTITY)
from substrate_framework.streamed_material_state import save_array
from substrate_framework.material_frozen_action import FrozenMaterialAction
from substrate_framework.full_frame_krylov import FullFrameKrylov, physical_invariants, memory_bound_bytes
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


def finite_json(value):
    """Overflow is an explicitly insufficient envelope, never a JSON NaN."""
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return 'nonfinite-insufficient'
    return value


def authority(old):
    """Authenticate historical source bytes without relabeling old code current."""
    records = []
    for entry in old['source_bindings']:
        path = Path(entry['path'])
        actual = binding(path)
        if actual['sha256'] != entry['sha256']:
            archive = HERE/'source-generations'/(path.stem+'-'+entry['sha256']+path.suffix)
            actual = binding(archive)
            if actual['sha256'] != entry['sha256'] or actual['bytes'] != int(entry['bytes']):
                raise ValueError('historical scientific source bytes unavailable: '+str(path))
        elif actual['bytes'] != int(entry['bytes']):
            raise ValueError('historical scientific source size differs')
        records.append(dict(historical=entry, authenticated_bytes=actual))
    return records


def ambient_growth(frozen):
    """Conservative two-sided anti-Hermitian norm bound, no PW-square matrix.

    Fourier insertion is an isometry after reciprocal normalizations cancel.
    Species: ||beta (I_atoms tensor antiHerm(D)) beta*|| <= ||beta||F^2
    ||antiHerm(D)||2. Prepared ionic legs use their whole operator norm as a
    safe (possibly loose) anti-Hermitian bound; no small term is discarded.
    """
    factors = frozen.h0.factors
    terms = dict(kinetic=float(np.max(np.abs(np.imag(factors.kinetic_ry)))),
                 local=float(np.max(np.abs(np.imag(frozen.field)))), species=[], ionic=[])
    for species in factors.species:
        anti = (species.active_d-species.active_d.conj().T)/(2j)
        value = float(np.linalg.norm(species.beta)**2*np.linalg.norm(anti, 2))
        terms['species'].append(dict(name=species.name, bound=value))
    for pu, bd, bra in frozen.terms:
        value = float(2*np.max(np.abs(pu))*np.linalg.norm(bd)*np.linalg.norm(bra))
        terms['ionic'].append(value)
    total = terms['kinetic']+terms['local']+sum(x['bound'] for x in terms['species'])+sum(terms['ionic'])
    return total, dict(terms=terms, bound=total,
                       qualification='analytic operator-factor bound evaluated in floating point, not interval enclosed')


def comparison(args, old):
    if args.comparison_receipt is None:
        return None
    receipt = explicit(args.comparison_receipt, args.comparison_sha256)
    data = bound_json(receipt)
    if (data.get('schema') != 'P254-streamed-sector-pilot-v1'
            or data.get('status') != 'completed' or data.get('phase') not in ('exponential', 'all')
            or not data.get('exponential_frame')):
        raise ValueError('settled successful full-PW exponential with retained frame required')
    for key in ('frame', 'original_D', 'sector'):
        if data[key] != old[key]:
            raise ValueError('comparison has different exact original '+key)
    for key in ('aggregate_receipt', 'manifest', 'checkpoint', 'initial_ions'):
        if data['inputs'][key] != old['inputs'][key]:
            raise ValueError('comparison has different exact physical input '+key)
    time = float(data['physical_interval_ry_inverse'])
    if time not in args.times:
        raise ValueError('comparison time must be explicitly requested in --times')
    authenticate(data['exponential_frame'])
    # Current or archived exact implementation authority, never a pending file.
    historical = []
    for entry in data['implementation_sources']:
        path = Path(entry['path'])
        actual = binding(path)
        if actual['sha256'] != entry['sha256']:
            path = HERE/'source-generations'/(path.stem+'-'+entry['sha256']+path.suffix)
            actual = binding(path)
        if actual['sha256'] != entry['sha256'] or actual['bytes'] != int(entry['bytes']):
            raise ValueError('comparison implementation bytes differ')
        historical.append(actual)
    return dict(receipt=receipt, time=time, frame=data['exponential_frame'],
                authenticated_implementation=historical,
                qualification='independent full-PW SciPy numerical comparison, not exact or interval-certified evolution')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('material-result', 'oracle'):
        parser.add_argument('--'+name, required=True, type=Path)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--times', required=True, nargs='+', type=float)
    parser.add_argument('--max-blocks', required=True, type=int)
    parser.add_argument('--max-basis', required=True, type=int)
    parser.add_argument('--action-cap', required=True, type=int, help='new H-vector columns, not calls; reused346 excluded')
    parser.add_argument('--absolute-error-target', required=True, type=float, help='absolute full-frame Frobenius target')
    parser.add_argument('--max-memory-mib', required=True, type=float, help='numerical array workspace only; excludes frozen operator caches')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--array-root', required=True, type=Path)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--bound-cells', type=int, default=4)
    parser.add_argument('--bound-degree', type=int, default=4)
    parser.add_argument('--deflation-rtol', type=float, default=1e-13)
    parser.add_argument('--comparison-receipt', type=Path)
    parser.add_argument('--comparison-sha256')
    args = parser.parse_args()
    if bool(args.comparison_receipt) != bool(args.comparison_sha256):
        parser.error('optional comparison requires both receipt and digest')
    if (not all(np.isfinite(args.times)) or len(set(args.times)) != len(args.times)
            or not np.isfinite(args.absolute_error_target) or args.absolute_error_target <= 0
            or not np.isfinite(args.max_memory_mib) or args.max_memory_mib <= 0
            or args.max_blocks < 1 or not 346 <= args.max_basis <= 18613
            or args.action_cap < 0 or args.bound_cells < 1 or args.bound_degree < 1
            or not 0 <= args.deflation_rtol < 1):
        parser.error('invalid finite positive limits, original-column basis or times')
    args.array_root = args.array_root.resolve()
    if args.array_root.is_relative_to(ROOT) or any((p/'.git').exists() for p in (args.array_root, *args.array_root.parents)):
        parser.error('--array-root must be outside Git')
    planned = memory_bound_bytes(18613, 346, args.max_basis)
    if planned > args.max_memory_mib*1024**2:
        parser.error('planned numerical workspace exceeds --max-memory-mib')
    os.sched_setaffinity(0, {args.cpu})
    started = perf_counter()
    source = sources()
    registration = json.loads((HERE/'block-krylov-registration.json').read_text())
    material_binding = explicit(args.material_result, args.material_result_sha256)
    if material_binding['sha256'] != registration['source_authority']['fused-action-material-001/result.json']:
        raise ValueError('preregistered actual retained full346 action report required')
    old = bound_json(material_binding)
    if old.get('schema') != 'P254-fused-action-material-v1' or old.get('status') != 'passed':
        raise ValueError('passing actual frozen full346 action receipt required')
    historical_sources = authority(old)
    # Mathematical operator authority is unchanged, even when execution wrappers change.
    for name in ('material_frozen_action.py', 'qe_material_hamiltonian.py'):
        key = 'src/substrate_framework/'+name
        if binding(ROOT/key)['sha256'] != registration['source_authority'][key]:
            raise ValueError('registered physical operator changed')
    oracle_binding = explicit(args.oracle, args.oracle_sha256)
    oracle = bound_json(oracle_binding)
    consumer_path = str(Path(__file__).resolve())
    original_science = {entry['path']: entry for entry in oracle['source_bindings']
                        if entry['path'] != consumer_path}
    current_science = {entry['path']: entry for entry in source
                       if entry['path'] != consumer_path}
    if (oracle.get('schema') != 'P254-full-frame-block-krylov-oracle-v1'
            or oracle.get('status') != 'passed' or original_science != current_science
            or not oracle['checks'] or not all(oracle['checks'].values())):
        raise ValueError('passing sensitive oracle for the unchanged numerical sources required')
    oracle_source_authentication = authority(oracle)
    authenticate(oracle['arrays'])
    inputs = old['inputs']
    for entry in (*inputs.values(), old['frame'], old['original_D'], old['full_action']):
        authenticate(entry)
    optional = comparison(args, old)
    args.output.mkdir(parents=True, exist_ok=False)
    args.array_root.mkdir(parents=True, exist_ok=False)
    publish(args.output/'started.json', dict(source_bindings=source, material=material_binding,
        oracle=oracle_binding, times=args.times, max_blocks=args.max_blocks, max_basis=args.max_basis,
        action_cap_new_vector_columns=args.action_cap, absolute_error_target=args.absolute_error_target,
        planned_array_workspace_bytes=planned, memory=memory(), optional_comparison=optional,
        cpu_affinity=sorted(os.sched_getaffinity(0))))
    preload_parsers()
    arrays, receipt, manifest, saved, original, grid, ordered, pairs = restore_inputs(
        inputs['aggregate_receipt'], inputs['manifest'])
    ions = bound_json(inputs['initial_ions'])
    spec = dict(aggregate=inputs['aggregate_receipt'], manifest=inputs['manifest'],
        preparation=manifest['preparation_terminal'], initial_ions=inputs['initial_ions'],
        scalar='PreparedReference joint-second-order; common Gamma coefficient only; actual per-slot H0')
    state = DiskMean(inputs['checkpoint']['path'], verify=False)
    known = inputs['checkpoint']['sha256'] == KNOWN_T0_SHA256 and state.identity == KNOWN_T0_IDENTITY
    if not known and state.identity != physical_identity(spec):
        raise ValueError('boundary is not the original physical model')
    if (state.time != 0. or not state.evidence.get('original_raw_frames_unchanged')
            or physical_specification(state.evidence['model']) != spec
            or not np.array_equal(state.u, ions['u_bohr']) or not np.array_equal(state.p, ions['p_canonical'])):
        raise ValueError('unchanged original physical initial boundary required')
    selected = [slot for slot in state.slots if slot.identity == old['sector']['identity']]
    if len(selected) != 1:
        raise ValueError('exact original sector absent')
    slot = selected[0]
    expected_sector = dict(old['sector'], shape=tuple(old['sector']['shape']))
    if (slot.frame != old['frame'] or slot.coefficient != old['original_D']
            or tuple(slot.shape) != (18613, 346) or slot.metadata() != expected_sector):
        raise ValueError('original full346 raw frame/D, physical basis or gauge differs')
    settings = NS(block_size=8, reality_tolerance=1e-10, max_h0_actions=args.max_blocks)
    model, background = build_evolution_model(settings, arrays, manifest, saved, original,
        grid, ordered, inputs['aggregate_receipt'], ions, state.identity)
    model.validate(state)
    stage = model.stage(arrays['conditional_density_bohr3'], arrays['positive_gradient_ry_bohr'].reshape(-1), state.u)
    psi, raw_d = state.arrays(slot)
    retained = np.load(old['full_action']['path'], mmap_mode='r', allow_pickle=False)
    initial_gram = psi.conj().T @ psi
    initial_invariants = physical_invariants(psi, raw_d)
    records, expansion_records = [], []
    start_operator = perf_counter()
    with model.operators(slot) as (coupling, counted_h0):
        # The production counted wrapper is preserved, not unwrapped for H calls.
        action = model.frozen_action(coupling, counted_h0, slot, stage)
        if counted_h0.material is None:
            raise ValueError('actual material-backed counted frozen operator required')
        diagnostic_operator = FrozenMaterialAction(counted_h0.material, coupling,
            potential=stage.potential[:, slot.component], direction=stage.u)
        growth, growth_receipt = ambient_growth(diagnostic_operator)
        del diagnostic_operator
        construction = dict(seconds=perf_counter()-start_operator, memory=memory())
        krylov = FullFrameKrylov(action, psi, initial_action=retained,
            max_blocks=args.max_blocks, max_basis=args.max_basis, action_cap=args.action_cap,
            max_memory_bytes=int(args.max_memory_mib*1024**2), ambient_growth_bound=growth,
            deflation_rtol=args.deflation_rtol)
        seed_arrays = dict(reconstruction=save_array(args.array_root/'initial-R.npy', krylov.initial_r),
            initial_qr_residual=save_array(args.array_root/'initial-QR-residual.npy', krylov.initial_qr_residual))
        while True:
            index = len(records)
            st = krylov.stage()
            stage_dir = args.array_root/f'stage-{index:03d}'
            stage_dir.mkdir()
            bound_arrays = {name: save_array(stage_dir/(name+'.npy'), value)
                for name, value in dict(Q=krylov.q, HQ=krylov.hq,
                    projected=st['projected'], coefficients=st['coefficients'],
                    residual=st['residual'], initial_residual=st['initial_residual']).items()}
            time_records = []
            for time_index, time in enumerate(args.times):
                evolved, error = krylov.propagate(time, bound_cells=args.bound_cells, bound_degree=args.bound_degree)
                gram_error = float(np.linalg.norm(evolved.conj().T @ evolved-initial_gram))
                output = save_array(stage_dir/f'frame-time-{time_index:03d}.npy', evolved)
                actual_comparison = None
                if optional is not None and time == optional['time']:
                    reference = np.load(optional['frame']['path'], mmap_mode='r', allow_pickle=False)
                    if reference.shape != psi.shape:
                        raise ValueError('comparison does not retain the exact full frame shape')
                    difference = float(np.linalg.norm(evolved-reference))
                    actual_comparison = dict(absolute_frobenius=difference,
                        relative_frobenius=difference/max(float(np.linalg.norm(reference)), np.finfo(float).tiny),
                        meets_requested_target=difference <= args.absolute_error_target,
                        receipt=optional['receipt'], qualification=optional['qualification'])
                meets = (np.isfinite(error['absolute_error_envelope_estimate'])
                         and error['absolute_error_envelope_estimate'] <= args.absolute_error_target
                         and (actual_comparison is None or actual_comparison['meets_requested_target']))
                time_records.append(dict(**error, status='estimated_target_met' if meets else 'insufficient',
                    frame=output, physical_invariants=physical_invariants(evolved, raw_d),
                    gram_change_frobenius=gram_error, actual_comparison=actual_comparison))
                del evolved
            record = dict(index=index, diagnostics=st['diagnostics'], arrays=bound_arrays, times=time_records,
                          memory=memory(), work_counts=model.counts, elapsed_seconds=perf_counter()-started)
            # Durable evaluated basis receipt precedes every next H expansion.
            publish(args.output/f'stage-{index:03d}.json', finite_json(record))
            records.append(record)
            if all(item['status'] == 'estimated_target_met' for item in time_records):
                stop_reason = 'estimated_target_met'
                break
            if krylov.blocks >= args.max_blocks or krylov.size >= args.max_basis:
                stop_reason = 'basis_limit'
                break
            publish(args.output/f'before-expansion-{index:03d}.json', dict(basis=krylov.size,
                new_H_vector_columns_remaining=args.action_cap-krylov.requested_columns,
                work_counts=model.counts, memory=memory()))
            expanded, discarded = krylov.expand()
            if discarded is not None:
                expanded['discarded_residual'] = save_array(stage_dir/'discarded-expansion.npy', discarded)
            expansion_records.append(expanded)
            publish(args.output/f'expansion-{index:03d}.json', expanded)
            if expanded['status'] != 'expanded':
                stop_reason = expanded['status']
                break
    for entry in (*inputs.values(), old['frame'], old['original_D'], old['full_action']):
        authenticate(entry)
    if sources() != source:
        raise ValueError('source changed during actual consumer')
    result = dict(schema='P254-full-frame-block-krylov-material-v1',
        status='estimated_target_met' if stop_reason == 'estimated_target_met' else 'insufficient',
        stop_reason=stop_reason, scope='one exact original full-PW frozen stage, every raw column and unchanged D',
        material=material_binding, oracle=oracle_binding, source_bindings=source,
        oracle_source_authentication=oracle_source_authentication,
        consumer_delta=binding(HERE/'block-krylov-consumer-delta-001.json'),
        historical_source_authentication=historical_sources, implementation_sources=implementation((Path(__file__),)),
        inputs=inputs, sector=slot.metadata(), frame=slot.frame, original_D=slot.coefficient,
        retained_seed_action=old['full_action'], background_receipt=background,
        original_invariants=initial_invariants, seed=krylov.seed_receipt, seed_arrays=seed_arrays,
        ambient_growth=growth_receipt, memory_plan=krylov.memory_plan, construction=construction,
        stages=records, expansions=expansion_records, action_calls=krylov.calls,
        new_requested_H_vector_columns=krylov.requested_columns, reused_H_vector_columns=346,
        production_work_counts=model.counts, absolute_error_target=args.absolute_error_target,
        optional_comparison=optional, seconds=perf_counter()-started, memory=memory(),
        cpu_affinity=sorted(os.sched_getaffinity(0)), native_threadpools=threadpool_info(),
        limitations=['Exact-arithmetic residual envelope evaluated in floating point; NOT an interval certificate',
          'Seed HQ solve and H action roundoff are diagnosed, not rigorously enclosed',
          'Basis reuse only for the identical frozen Hamiltonian; changed density/ionic stage invalidates HQ',
          'Does not accelerate or replace global all-sector density/mixed-moment reductions',
          'No full128 coverage, nonlinear midpoint, coherence, long-time physical validation or speedup claim',
          'Raw original D is unchanged; no occupation refilling or physical band/row truncation',
          'Memory planning excludes frozen operator caches, native allocator overhead and disk arrays',
          'Only selected frame/D rehashed; complete original boundary and source inputs authenticated'])
    publish(args.output/'result.json', finite_json(result))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
