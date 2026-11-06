#!/usr/bin/env python3
"""Verify detached native metadata and consume the actual whole moving history.

Retains existing numerical q/field proofs when their code is unchanged, with
only the native reader's metadata copies and its explicit consumer binding new.
Checks actual saved data, omitted-copy mutations and the same four-GiB full
history intake that previously failed. No q response is run.
"""
import argparse
import ast
from dataclasses import fields, replace
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from run_joint_q_material import bindings, selected_initialization
from matched_joint_case import load_case
from substrate_framework.material_control_reference import bound_json
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.partitioned_q_response import partition_plan
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.streamed_material_state import file_binding


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    for name in ('initialization', 'history', 'retained-oracle', 'prior-native-source', 'prior-q-source'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--time', type=float, required=True)
    parser.add_argument('--second-time', type=float, required=True)
    parser.add_argument('--subdivisions', type=int, required=True)
    parser.add_argument('--workers', type=int, choices=range(1, 6), required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    return parser.parse_args()


def checked_binding(path, digest):
    binding = file_binding(path)
    if binding['sha256'] != digest:
        raise ValueError('explicit input digest differs: '+str(path))
    return binding


def unchanged_numerics(previous, current, component):
    old, new = ast.parse(previous), ast.parse(current)
    if component == 'native':
        function = next(node for node in new.body if isinstance(node, ast.FunctionDef) and node.name == 'load_wavefunctions')
        copies, freezes, retained = set(), set(), []
        for node in function.body:
            if (isinstance(node, ast.Assign) and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name) and node.targets[0].id in ('reciprocal', 'miller')):
                call = node.value
                if (not isinstance(call, ast.Call) or not isinstance(call.func, ast.Attribute)
                        or call.func.attr != 'copy' or call.args
                        or ast.dump(ast.Tuple(elts=call.keywords, ctx=ast.Load())) != ast.dump(ast.Tuple(elts=[ast.keyword(arg='order', value=ast.Constant(value='F'))], ctx=ast.Load()))):
                    return False
                copies.add(node.targets[0].id)
                node.value = call.func.value
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                call = node.value
                if (isinstance(call.func, ast.Attribute) and isinstance(call.func.value, ast.Name)
                        and call.func.value.id in ('reciprocal', 'miller') and call.func.attr == 'setflags'):
                    if call.args or len(call.keywords) != 1 or call.keywords[0].arg != 'write' or not isinstance(call.keywords[0].value, ast.Constant) or call.keywords[0].value.value is not False:
                        return False
                    freezes.add(call.func.value.id)
                    continue
            retained.append(node)
        if copies != {'reciprocal', 'miller'} or freezes != copies:
            return False
        function.body = retained
    elif component == 'q-consumer':
        assignment = next(node for node in new.body if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'MODULES' for t in node.targets))
        additions = [node for node in assignment.value.elts if isinstance(node, ast.Constant) and node.value == 'qe_saved_native.py']
        if len(additions) != 1:
            return False
        assignment.value.elts.remove(additions[0])
    else:
        raise ValueError('unknown source comparison')
    return ast.dump(old, include_attributes=False) == ast.dump(new, include_attributes=False)


def retained_metadata_bytes(basis):
    roots = {}
    for name in ('reciprocal_bohr_inv', 'miller'):
        root = getattr(basis, name)
        while True:
            if isinstance(root, np.ndarray) and root.base is not None:
                root = root.base
            elif isinstance(root, memoryview):
                root = root.obj
            else:
                break
        if isinstance(root, np.ndarray):
            size = root.nbytes
        elif isinstance(root, (bytes, bytearray)):
            size = len(root)
        else:
            raise TypeError('unrecognized retained metadata backing owner')
        roots[id(root)] = size
    return sum(roots.values())


def same_saved_values(a, b):
    for field in fields(a):
        left, right = getattr(a, field.name), getattr(b, field.name)
        if isinstance(left, np.ndarray):
            if (left.dtype != right.dtype or left.shape != right.shape
                    or left.tobytes(order='A') != right.tobytes(order='A')):
                return False
        elif left != right:
            return False
    return True


def execute(args, budget):
    retained_binding = checked_binding(args.retained_oracle, args.retained_oracle_sha256)
    retained = bound_json(retained_binding)
    prior_binding = checked_binding(args.prior_native_source, args.prior_native_source_sha256)
    prior_q_binding = checked_binding(args.prior_q_source, args.prior_q_source_sha256)
    native_path, q_path = ROOT/'src/substrate_framework/qe_saved_native.py', HERE/'run_joint_q_material.py'
    current_sources = bindings()
    old_sources = {item['path']: item for item in retained['source_bindings']}
    changed = [item['path'] for item in current_sources if old_sources.get(item['path'], {}).get('sha256') != item['sha256']]
    if (retained['schema'] != 'P254-joint-partitioned-q-oracle-v1' or retained['status'] != 'passed'
            or retained['extension'] != 'physical-global-field-readouts'
            or set(changed) != {str(native_path), str(q_path)}
            or prior_q_binding['sha256'] != old_sources[str(q_path)]['sha256']):
        raise ValueError('retained numerical proof plus only native ownership and explicit reader inventory changes required')
    checks = []
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        if not passed:
            raise AssertionError(name)
    check('native_decoding_unchanged_outside_owned_readonly_metadata_copies',
          unchanged_numerics(args.prior_native_source.read_bytes(), native_path.read_bytes(), 'native'))
    check('q_numerics_unchanged_except_explicit_native_reader_inventory',
          unchanged_numerics(args.prior_q_source.read_bytes(), q_path.read_bytes(), 'q-consumer'))
    initialization, saved = selected_initialization(args)
    hb = checked_binding(args.history, args.history_sha256)
    context = load_case(args.case)
    if saved['model_identity'] != context.identity or saved['specification'] != context.specification:
        raise ValueError('actual saved case identity/specification differs')
    spec = importlib.util.spec_from_file_location('p254_prior_native_reader', args.prior_native_source)
    prior_reader = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = prior_reader
    spec.loader.exec_module(prior_reader)
    selected_source = dict(context.ordered[0].source_binding)
    original = prior_reader.load_wavefunctions(selected_source)
    detached = load_wavefunctions(selected_source)
    small_bytes = sum(getattr(original, key).nbytes for key in ('reciprocal_bohr_inv', 'miller'))
    old_bytes, new_bytes = retained_metadata_bytes(original), retained_metadata_bytes(detached)
    check('actual_native_values_dtypes_bytes_and_authorities_preserved', same_saved_values(original, detached))
    check('detached_metadata_remains_readonly', all(not getattr(detached, key).flags.writeable for key in ('reciprocal_bohr_inv', 'miller')))
    check('retained_metadata_fits_its_actual_data_bytes', new_bytes <= small_bytes)
    check('no_detachment_mutation_rejected', old_bytes > small_bytes)
    for key in ('reciprocal_bohr_inv', 'miller'):
        mutation = replace(detached, **{key: getattr(original, key)})
        check('omitted_'+key+'_copy_rejected', retained_metadata_bytes(mutation) > small_bytes)
    del mutation, detached, original
    print('Actual ownership/value checks passed; consuming all original history sectors.', flush=True)
    history = context.make_history(hb, budget, block_size=8)
    if file_binding(history.initial.path) != saved['joint_boundary']:
        raise ValueError('actual history does not start at the selected initialization')
    response = context.bind_response(history, np.zeros(3))
    times = history.time_grid(args.time, args.subdivisions)
    second_times = history.time_grid(args.second_time, args.subdivisions)
    aliases = [state.time for state in history.priors]+[history.final.time]
    plans = {mode: partition_plan(response, times, second_times, args.workers,
        mode=mode, aliased_frame_times=aliases, dynamics='screened-joint') for mode in ('forward', 'pair', 'all')}
    budget.check_processes(os.getpid(), [])
    check('whole_actual_history_and_q_plan_complete_without_h0_actions', budget.actions == 0)
    stat = os.statvfs(args.output)
    return dict(schema='P254-joint-partitioned-q-oracle-v1', status='passed',
        extension='physical-global-field-readouts', material_result=False,
        source_bindings=current_sources, retained_proof=retained_binding,
        retained_numerical_checks_not_rerun=dict(base=retained['retained_check_count'], field=retained['new_check_count']),
        proof_sources=[file_binding(Path(__file__).resolve()), prior_binding, prior_q_binding],
        native_source=selected_source,
        original_native_parser_provenance='Pre-repair bytes captured before this correction, not retrospectively represented as part of the earlier numerical-proof source inventory. Current real consumer now binds its native reader explicitly.',
        checks=checks, new_check_count=len(checks), case=args.case, initialization=initialization,
        history=hb, model_identity=context.identity, explicit_sector_count=len(history.bases),
        ownership=dict(original_backing_bytes=old_bytes, detached_backing_bytes=new_bytes, actual_metadata_bytes=small_bytes),
        q=[0, 0, 0], time_grid=times.tolist(), second_time_grid=second_times.tolist(),
        resource_plans=plans, filesystem_free_bytes=stat.f_bavail*stat.f_frsize,
        numerical_accuracy_certified=False,
        scope='Existing numerical q/field proofs retained by exact unchanged-code and actual data identity; not re-executed or relabelled. New checks cover native-buffer ownership, omitted-copy sensitivity and the genuine full-material history/planner under the same memory cap that previously failed. No propagated q response, spatial result or physical accuracy certificate.')


def main():
    args = parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1, max_seconds=args.max_seconds, max_bytes=1048576,
        max_memory_bytes=args.max_memory_bytes, worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(524288)
    result = dict(schema='P254-joint-history-memory-check-v1', status='failed')
    started = perf_counter()
    try:
        result = execute(args, budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(seconds=perf_counter()-started, budget=budget.snapshot(),
        process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        cpu_affinity=sorted(os.sched_getaffinity(0)))
    result['budget']['time_cap'] = 'Cooperative intake checks and external systemd hard wall supervisor; no worker watchdog in this standalone consumer.'
    result['budget']['memory_accounting'] = 'Parent RSS sampled after full intake; whole-invocation RSS peak reported separately; external systemd hard memory limit.'
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
