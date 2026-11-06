#!/usr/bin/env python3
"""Bounded original full346 column-Arnoldi LEGACY frozen-stage benchmark.

This reproduces the old block/SciPy operator, not the newer physical joint
Hamiltonian or a coupled material trajectory. All requested times share only
this identical authenticated frozen operator. Large arrays remain outside Git.
"""
import ast
import argparse
import json
import os
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace as NS

from run_column_krylov_oracle import (
    HERE, ROOT, np, binding, authenticate, explicit, publish, sources)
from run_block_krylov_material import (
    memory, finite_json, authority, ambient_growth, comparison,
    bound_json, preload_parsers, restore_inputs, build_evolution_model,
    DiskMean, physical_identity, physical_specification, KNOWN_T0_SHA256,
    KNOWN_T0_IDENTITY, save_array, FrozenMaterialAction, physical_invariants,
    threadpool_limits, threadpool_info)
from substrate_framework.column_frame_krylov import ColumnFrameKrylov, memory_bound_bytes


def operator_norm(frozen):
    """Triangle/operator-factor norm bound for seed normalization correction.

    Fourier insertion is an isometry after reciprocal normalizations cancel.
    The two directional nonlocal legs are bounded independently; no Hermitian
    simplification or numerical cutoff. Ordinary floating evaluation only.
    """
    factors = frozen.h0.factors
    terms = dict(kinetic=float(np.max(np.abs(factors.kinetic_ry))),
        local=float(np.max(np.abs(frozen.field))), species=[], ionic=[])
    for species in factors.species:
        terms['species'].append(float(np.linalg.norm(species.beta)**2*np.linalg.norm(species.active_d, 2)))
    for pu, bd, bra in frozen.terms:
        terms['ionic'].append(float(2*np.max(np.abs(pu))*np.linalg.norm(bd)*np.linalg.norm(bra)))
    total = terms['kinetic']+terms['local']+sum(terms['species'])+sum(terms['ionic'])
    return total, dict(bound=total, terms=terms, qualification='ordinary floating operator-factor norm bound')

def authenticate_frozen_composition(old, historical):
    """Require unchanged operator bytes and unchanged historical stage builder.

    Historical runtime wrappers may differ. AST equality below ignores source
    coordinates only, not expressions; no numerical similarity licenses reuse.
    """
    by_path = {entry['historical']['path']: entry for entry in historical}
    operators = ('material_frozen_action.py', 'qe_material_hamiltonian.py',
        'qe_material_dyads.py', 'qe_bloch_ionic_action.py', 'qe_material_ionic_action.py',
        'qe_bloch_qzero.py', 'qe_material_grid.py', 'qe_bloch_dyads.py')
    records = []
    for name in operators:
        path = ROOT/'src/substrate_framework'/name
        entry = by_path[str(path)]
        if binding(path)['sha256'] != entry['historical']['sha256']:
            raise ValueError('retained action operator source changed: '+name)
        records.append(dict(path=str(path), identity='exact historical source bytes'))
    for path, names in (
            (HERE/'run_streamed_material_evolution.py', ('build_evolution_model',)),
            (ROOT/'src/substrate_framework/streamed_material_evolution.py',
             ('StreamedMeanEvolution.stage', 'real_physical'))):
        historic = Path(by_path[str(path)]['authenticated_bytes']['path'])
        def definitions(source):
            tree = ast.parse(source.read_text())
            found = {}
            for node in tree.body:
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    found[node.name] = ast.dump(node, include_attributes=False)
                elif isinstance(node, ast.ClassDef):
                    for method in node.body:
                        if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef)):
                            found[node.name+'.'+method.name] = ast.dump(method, include_attributes=False)
            return found
        current, previous = definitions(path), definitions(historic)
        for name in names:
            if name not in current or current[name] != previous.get(name):
                raise ValueError('retained seed frozen-stage definition changed: '+name)
        records.append(dict(path=str(path), historical_path=str(historic), identical_definitions=list(names)))
    return records



def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('material-result', 'oracle', 'comparison-receipt'):
        parser.add_argument('--'+name, required=True, type=Path)
    for name in ('material-result-sha256', 'oracle-sha256', 'comparison-sha256'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--times', required=True, nargs='+', type=float)
    parser.add_argument('--max-dimension', required=True, type=int, help='Arnoldi dimension per original column')
    parser.add_argument('--action-cap', required=True, type=int, help='new H-vector columns; retained346 seed excluded')
    parser.add_argument('--absolute-error-target', required=True, type=float)
    parser.add_argument('--max-memory-mib', required=True, type=float, help='numerical arrays only, not total RSS')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--array-root', required=True, type=Path)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--bound-cells', type=int, default=4)
    parser.add_argument('--bound-degree', type=int, default=4)
    parser.add_argument('--breakdown-rtol', type=float, default=1e-13)
    parser.add_argument('--retain-final-basis', action='store_true', help='save per-column full Q/HQ for independent audit')
    args = parser.parse_args()
    if (not all(np.isfinite(args.times)) or len(set(args.times)) != len(args.times)
            or not np.isfinite(args.absolute_error_target) or args.absolute_error_target <= 0
            or not np.isfinite(args.max_memory_mib) or args.max_memory_mib <= 0
            or not 1 <= args.max_dimension <= 18613 or args.action_cap < 0
            or args.bound_cells < 1 or args.bound_degree < 1
            or not 0 <= args.breakdown_rtol < 1):
        parser.error('invalid finite times or hard numerical limits')
    args.array_root = args.array_root.resolve()
    if args.array_root.is_relative_to(ROOT) or any((p/'.git').exists() for p in (args.array_root, *args.array_root.parents)):
        parser.error('--array-root must be outside Git')
    args.planned = memory_bound_bytes(18613, 346, args.max_dimension)
    if args.planned > args.max_memory_mib*1024**2:
        parser.error(f'planned numerical array workspace {args.planned} exceeds supplied memory budget')
    return args


def evaluate_times(args, model, psi, raw_d, initial_gram, optional, directory):
    records = []
    for index, time in enumerate(args.times):
        evolved, error = model.propagate(time, bound_cells=args.bound_cells, bound_degree=args.bound_degree)
        actual = None
        if time == optional['time']:
            reference = np.load(optional['frame']['path'], mmap_mode='r', allow_pickle=False)
            if reference.shape != psi.shape:
                raise ValueError('comparison omitted original rows or columns')
            difference = float(np.linalg.norm(evolved-reference))
            actual = dict(absolute_frobenius=difference,
                relative_frobenius=difference/max(float(np.linalg.norm(reference)), np.finfo(float).tiny),
                meets_requested_target=difference <= args.absolute_error_target,
                receipt=optional['receipt'], qualification=optional['qualification'])
        meets = (np.isfinite(error['absolute_error_envelope_estimate'])
                 and error['absolute_error_envelope_estimate'] <= args.absolute_error_target
                 and (actual is None or actual['meets_requested_target']))
        records.append(dict(**error, status='estimated_target_met' if meets else 'insufficient',
            frame=save_array(directory/f'frame-time-{index:03d}.npy', evolved),
            physical_invariants=physical_invariants(evolved, raw_d),
            gram_change_frobenius=float(np.linalg.norm(evolved.conj().T @ evolved-initial_gram)),
            actual_comparison=actual))
    return records


def main():
    args = parse_args()
    os.sched_setaffinity(0, {args.cpu})
    started, source = perf_counter(), sources()
    registration = json.loads((HERE/'column-krylov-registration.json').read_text())
    material_binding = explicit(args.material_result, args.material_result_sha256)
    if material_binding['sha256'] != registration['source_authority']['fused-action-material-001/result.json']:
        raise ValueError('registered legacy retained full346 action receipt required')
    old = bound_json(material_binding)
    if old.get('schema') != 'P254-fused-action-material-v1' or old.get('status') != 'passed':
        raise ValueError('passing legacy frozen action material receipt required')
    historical = authority(old)
    frozen_authentication = authenticate_frozen_composition(old, historical)
    for name in ('material_frozen_action.py', 'qe_material_hamiltonian.py'):
        key = 'src/substrate_framework/'+name
        if binding(ROOT/key)['sha256'] != registration['source_authority'][key]:
            raise ValueError('registered legacy operator source changed')
    oracle_binding = explicit(args.oracle, args.oracle_sha256)
    oracle = bound_json(oracle_binding)
    if (oracle.get('schema') != 'P254-column-frame-krylov-oracle-v1'
            or oracle.get('status') != 'passed' or oracle['source_bindings'] != source
            or not oracle['checks'] or not all(oracle['checks'].values())):
        raise ValueError('passing independent oracle bound to every unchanged source required')
    authenticate(oracle['arrays'])
    inputs = old['inputs']
    for entry in (*inputs.values(), old['frame'], old['original_D'], old['full_action']):
        authenticate(entry)
    for key, name in (('frame', 'original_frame'), ('original_D', 'original_D'),
                      ('full_action', 'fused-action-material-001/new-full346-action.npy')):
        if old[key]['sha256'] != registration['source_authority'][name]:
            raise ValueError('registered exact original '+key+' differs')
    optional = comparison(args, old)
    args.output.mkdir(parents=True, exist_ok=False)
    args.array_root.mkdir(parents=True, exist_ok=False)
    publish(args.output/'started.json', dict(source_bindings=source, material=material_binding,
        oracle=oracle_binding, times=args.times, max_dimension_per_column=args.max_dimension,
        action_cap_new_vector_columns=args.action_cap, absolute_error_target=args.absolute_error_target,
        planned_array_workspace_bytes=args.planned, memory=memory(), scipy_comparison=optional,
        cpu_affinity=sorted(os.sched_getaffinity(0)), scope='legacy identical frozen stage only'))
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
        raise ValueError('boundary differs from authenticated legacy model')
    if (state.time != 0. or not state.evidence.get('original_raw_frames_unchanged')
            or physical_specification(state.evidence['model']) != spec
            or not np.array_equal(state.u, ions['u_bohr']) or not np.array_equal(state.p, ions['p_canonical'])):
        raise ValueError('unchanged original initial boundary required')
    selected = [slot for slot in state.slots if slot.identity == old['sector']['identity']]
    if len(selected) != 1:
        raise ValueError('exact original sector absent')
    slot = selected[0]
    if (slot.frame != old['frame'] or slot.coefficient != old['original_D']
            or tuple(slot.shape) != (18613, 346)
            or slot.metadata() != dict(old['sector'], shape=tuple(old['sector']['shape']))):
        raise ValueError('full346 original raw frame/D or gauge differs')
    settings = NS(block_size=8, reality_tolerance=1e-10, max_h0_actions=args.max_dimension+1)
    model, background = build_evolution_model(settings, arrays, manifest, saved, original,
        grid, ordered, inputs['aggregate_receipt'], ions, state.identity)
    model.validate(state)
    stage = model.stage(arrays['conditional_density_bohr3'], arrays['positive_gradient_ry_bohr'].reshape(-1), state.u)
    psi, raw_d = state.arrays(slot)
    retained = np.load(old['full_action']['path'], mmap_mode='r', allow_pickle=False)
    initial_gram, initial_invariants = psi.conj().T @ psi, physical_invariants(psi, raw_d)
    records, expansions, basis_arrays = [], [], []
    with model.operators(slot) as (coupling, counted_h0):
        action = model.frozen_action(coupling, counted_h0, slot, stage)
        if counted_h0.material is None:
            raise ValueError('counted actual material operator required')
        diagnostic = FrozenMaterialAction(counted_h0.material, coupling,
            potential=stage.potential[:, slot.component], direction=stage.u)
        growth, growth_receipt = ambient_growth(diagnostic)
        norm, norm_receipt = operator_norm(diagnostic)
        del diagnostic
        identity = material_binding['sha256']
        krylov = ColumnFrameKrylov(action, psi, initial_action=retained,
            operator_identity=identity, initial_action_identity=identity,
            max_dimension=args.max_dimension, action_cap=args.action_cap,
            max_memory_bytes=int(args.max_memory_mib*1024**2), ambient_growth_bound=growth,
            operator_norm_bound=norm, breakdown_rtol=args.breakdown_rtol)
        while True:
            index = len(records)
            st = krylov.stage()
            directory = args.array_root/f'stage-{index:03d}'
            directory.mkdir()
            projected = {f'column-{j:03d}': t for j, t in enumerate(st['projected'])}
            np.savez(directory/'projected.npz', **projected)
            time_records = evaluate_times(args, krylov, psi, raw_d, initial_gram, optional, directory)
            record = dict(index=index, diagnostics=st['diagnostics'], times=time_records,
                projected=binding(directory/'projected.npz'), memory=memory(),
                work_counts=dict(model.counts), elapsed_seconds=perf_counter()-started)
            # An evaluated full-frame result for every requested time is durable
            # before another H batch. Insufficiency never relaxes the target.
            publish(args.output/f'stage-{index:03d}.json', finite_json(record))
            records.append(record)
            if all(item['status'] == 'estimated_target_met' for item in time_records):
                stop_reason = 'estimated_target_met'
                break
            publish(args.output/f'before-expansion-{index:03d}.json', dict(dimensions=krylov.sizes.tolist(),
                new_H_vector_columns_remaining=args.action_cap-krylov.requested_columns,
                work_counts=dict(model.counts), memory=memory()))
            expansion = krylov.expand()
            expansions.append(expansion)
            publish(args.output/f'expansion-{index:03d}.json', finite_json(expansion))
            if expansion['status'] != 'expanded':
                stop_reason = expansion['status']
                break
        if args.retain_final_basis:
            for j, k in enumerate(krylov.sizes):
                basis_arrays.append(dict(column=j, dimension=int(k),
                    Q=save_array(args.array_root/f'final-Q-{j:03d}.npy', krylov.q[j, :k].T),
                    HQ=save_array(args.array_root/f'final-HQ-{j:03d}.npy', krylov.hq[j, :k].T)))
    for entry in (*inputs.values(), old['frame'], old['original_D'], old['full_action'], optional['frame']):
        authenticate(entry)
    if sources() != source:
        raise ValueError('source changed during material benchmark')
    result = dict(schema='P254-column-frame-krylov-material-v1',
        status='estimated_target_met' if stop_reason == 'estimated_target_met' else 'insufficient',
        stop_reason=stop_reason, scope='legacy one-sector identical frozen stage, NOT physical joint dynamics',
        material=material_binding, oracle=oracle_binding, source_bindings=source,
        historical_source_authentication=historical, inputs=inputs, sector=slot.metadata(),
        frozen_composition_authentication=frozen_authentication,
        frame=slot.frame, original_D=slot.coefficient, retained_seed_action=old['full_action'],
        background_receipt=background, original_invariants=initial_invariants,
        seed=krylov.seed_receipt, ambient_growth=growth_receipt, operator_norm=norm_receipt,
        memory_plan=krylov.memory_plan, stages=records, expansions=expansions,
        final_basis_arrays=basis_arrays, action_calls=krylov.calls,
        new_requested_H_vector_columns=krylov.requested_columns, reused_H_vector_columns=346,
        production_work_counts=model.counts, absolute_error_target=args.absolute_error_target,
        scipy_comparison=optional, seconds=perf_counter()-started, memory=memory(),
        cpu_affinity=sorted(os.sched_getaffinity(0)), native_threadpools=threadpool_info(),
        limitations=['Legacy common-Gamma frozen-stage comparison; not the newer physical all-source joint operator',
            'Every original row/column and complex raw D retained; no Gram repair or saved-band truncation',
            'Full-row residual includes decomposition, initial reconstruction and breakdown; no interval certificate',
            'Seed scaling defect uses an analytic operator-norm bound; H-action/expm/norm roundoff not enclosed',
            'Reusing seed or basis at a changed ionic/electronic stage is invalid',
            'Memory cap is numerical arrays, not process RSS; excludes material operator caches and allocators',
            'No claim of coupled trajectory, all128-slot evolution, physical loss, long-time accuracy or speedup'])
    publish(args.output/'result.json', finite_json(result))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
