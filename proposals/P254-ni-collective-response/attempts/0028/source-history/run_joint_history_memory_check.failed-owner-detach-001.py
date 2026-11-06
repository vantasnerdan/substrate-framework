#!/usr/bin/env python3
"""Verify detached native metadata and consume the actual whole moving history.

Retains existing numerical q/field proofs only when the source change is exactly
the metadata-detachment wrapper and every other producer source is unchanged.
The new checks use actual material basis data, omitted-copy mutations and the
same four-GiB full-history intake that previously failed. No q response is run.
"""
import argparse
import ast
from dataclasses import fields, replace
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
from substrate_framework.material_joint_response import _detached_basis
from substrate_framework.partitioned_q_response import partition_plan
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.streamed_material_state import file_binding


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    for name in ('initialization', 'history', 'retained-oracle', 'prior-history-source'):
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


def unchanged_numerics(previous, current):
    old, new = ast.parse(previous), ast.parse(current)
    helpers = [node for node in new.body if isinstance(node, ast.FunctionDef) and node.name == '_detached_basis']
    if len(helpers) != 1:
        raise ValueError('one explicit metadata-detachment helper required')
    new.body.remove(helpers[0])
    for node in new.body:
        if isinstance(node, ast.ImportFrom) and node.module == 'dataclasses':
            node.names = [alias for alias in node.names if alias.name != 'replace']
    class RemoveWrapper(ast.NodeTransformer):
        count = 0
        def visit_Call(self, node):
            node = self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == '_detached_basis':
                if len(node.args) != 1 or node.keywords:
                    raise ValueError('unexpected detachment call')
                self.count += 1
                return node.args[0]
            return node
    normalization = RemoveWrapper()
    new = normalization.visit(new)
    return normalization.count == 1 and ast.dump(old, include_attributes=False) == ast.dump(new, include_attributes=False)


def retained_metadata_bytes(basis):
    roots = {}
    for name in ('reciprocal_bohr_inv', 'source_miller'):
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


def same_basis_values(a, b):
    for field in fields(a):
        left, right = getattr(a, field.name), getattr(b, field.name)
        if isinstance(left, np.ndarray):
            if (left.dtype != right.dtype or left.shape != right.shape
                    or left.tobytes(order='A') != right.tobytes(order='A')):
                return False
        elif left is not right:
            return False
    return True


def execute(args, budget):
    retained_binding = checked_binding(args.retained_oracle, args.retained_oracle_sha256)
    retained = bound_json(retained_binding)
    prior_binding = checked_binding(args.prior_history_source, args.prior_history_source_sha256)
    history_path = ROOT/'src/substrate_framework/material_joint_response.py'
    current_sources = bindings()
    old_sources = {item['path']: item for item in retained['source_bindings']}
    changed = [item['path'] for item in current_sources if old_sources.get(item['path'], {}).get('sha256') != item['sha256']]
    if (retained['schema'] != 'P254-joint-partitioned-q-oracle-v1' or retained['status'] != 'passed'
            or retained['extension'] != 'physical-global-field-readouts'
            or changed != [str(history_path)]
            or prior_binding['sha256'] != old_sources[str(history_path)]['sha256']):
        raise ValueError('exact retained numerical proof and only the declared history source change required')
    checks = []
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        if not passed:
            raise AssertionError(name)
    check('all_previous_numerical_code_unchanged_outside_detachment',
          unchanged_numerics(args.prior_history_source.read_bytes(), history_path.read_bytes()))
    initialization, saved = selected_initialization(args)
    hb = checked_binding(args.history, args.history_sha256)
    context = load_case(args.case)
    if saved['model_identity'] != context.identity:
        raise ValueError('actual saved case identity differs')
    frame = context.grid.load_frame(context.ordered[0])
    original, detached = frame.basis, _detached_basis(frame.basis)
    small_bytes = sum(getattr(original, key).nbytes for key in ('reciprocal_bohr_inv', 'source_miller'))
    old_bytes, new_bytes = retained_metadata_bytes(original), retained_metadata_bytes(detached)
    check('actual_basis_values_and_authorities_preserved', same_basis_values(original, detached))
    check('detached_metadata_remains_readonly', all(not getattr(detached, key).flags.writeable for key in ('reciprocal_bohr_inv', 'source_miller')))
    check('retained_metadata_fits_its_actual_data_bytes', new_bytes <= small_bytes)
    check('no_detachment_mutation_rejected', old_bytes > small_bytes)
    for key in ('reciprocal_bohr_inv', 'source_miller'):
        mutation = replace(detached, **{key: getattr(original, key)})
        check('omitted_'+key+'_copy_rejected', retained_metadata_bytes(mutation) > small_bytes)
    del mutation, detached, original, frame
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
        proof_sources=[file_binding(Path(__file__).resolve()), prior_binding],
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
