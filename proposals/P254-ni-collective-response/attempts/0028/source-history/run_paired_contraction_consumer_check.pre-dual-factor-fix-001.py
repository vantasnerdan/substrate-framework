#!/usr/bin/env python3
"""Exercise old/new real q-edge consumers on the genuine moving material history.

Retains settled numerical/component proofs. Only the selected edge is evaluated;
its original full PW rows, every saved column and raw D are kept. This is a
nonzero complex operator comparison, not a propagated material response.
"""
import argparse
import ast
from copy import deepcopy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from run_joint_q_material import bindings, selected_initialization, bind_saved_model
from matched_joint_case import load_case
from matched_joint_state import preparation_arrays
from substrate_framework.material_control_reference import bound_json, bound_arrays
from substrate_framework.partitioned_q_response import PartitionedQResponse
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.streamed_material_state import authenticate, file_binding


def checked(path, expected):
    binding = file_binding(path)
    if binding['sha256'] != expected:
        raise ValueError('explicit input digest differs: '+str(path))
    return binding, bound_json(binding)


def dump(tree):
    return ast.dump(tree, include_attributes=False)


def unchanged_other_edge_code(before, after):
    old, new = ast.parse(before), ast.parse(after)
    for tree in (old, new):
        owner = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'PartitionedQResponse')
        owner.body = [n for n in owner.body if not isinstance(n, ast.FunctionDef) or n.name != 'independent_edge']
    return dump(old) == dump(new)


def unchanged_physical_stage(before, after):
    old, new = ast.parse(before), ast.parse(after)
    new.body = [n for n in new.body if not isinstance(n, ast.ImportFrom)
                or n.module not in ('paired_finiteq_moments', 'paired_field_moments')]
    owner = next(n for n in new.body if isinstance(n, ast.ClassDef) and n.name == 'MaterialJointHistory')
    owner.body = [n for n in owner.body if not isinstance(n, ast.FunctionDef)
                  or n.name not in ('paired_moments_and_drive', 'paired_field_moments_and_drive')]
    class RestoreAxes(ast.NodeTransformer):
        count = 0
        def visit_Call(self, node):
            node = self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'CartesianPhysicalFieldAxes':
                if len(node.args) != 1 or node.keywords:
                    raise ValueError('unexpected physical axes construction')
                self.count += 1
                return ast.Attribute(value=ast.Call(func=ast.Attribute(value=ast.Name(id='CartesianFieldActions', ctx=ast.Load()), attr='from_vertex', ctx=ast.Load()), args=node.args, keywords=[]), attr='axes', ctx=ast.Load())
            return node
    normalize = RestoreAxes()
    new = normalize.visit(new)
    return normalize.count == 1 and dump(old) == dump(new)


def retained_sources(args, check):
    retained_binding, retained = checked(args.retained_oracle, args.retained_oracle_sha256)
    component_bindings, components = [], []
    for name in ('matter', 'field'):
        binding, value = checked(getattr(args, name+'_oracle'), getattr(args, name+'_oracle_sha256'))
        if value['status'] != 'passed' or any(c['passed'] is not True for c in value['checks']):
            raise ValueError('complete passed physical component proof required')
        component_bindings.append(binding); components.append(value)
    current = bindings()
    old = {b['path']: b for b in retained['source_bindings']}
    now = {b['path']: b for b in current}
    modified = {str(ROOT/'src/substrate_framework'/name) for name in ('material_joint_response.py', 'partitioned_q_response.py')}
    modified.add(str(HERE/'run_joint_q_material.py'))
    added = {str(ROOT/'src/substrate_framework'/name) for name in ('paired_finiteq_moments.py', 'paired_field_moments.py')}
    added.update(str(HERE/name) for name in ('run_paired_finiteq_moments_check.py', 'run_paired_field_moments_check.py', 'run_paired_contraction_consumer_check.py'))
    check('all_other_production_sources_retained_exactly', set(now)-set(old) == added
          and set(old).issubset(now)
          and {p for p in old if old[p]['sha256'] != now[p]['sha256']} == modified)
    snapshots = []
    for module, validator in (('material_joint_response.py', unchanged_physical_stage), ('partitioned_q_response.py', unchanged_other_edge_code)):
        path = ROOT/'src/substrate_framework'/module
        prior = HERE/'source-history'/(module.removesuffix('.py')+'.pre-paired-edges-001.py')
        binding = file_binding(prior)
        check(module+'_original_source_preserved', binding['sha256'] == old[str(path)]['sha256'])
        check(module+'_unchanged_outside_declared_native_hooks', validator(prior.read_bytes(), path.read_bytes()))
        snapshots.append(binding)
    previous_field = HERE/'source-history/paired_field_moments.pre-native-axes-001.py'
    field_path = ROOT/'src/substrate_framework/paired_field_moments.py'
    previous_ast, current_ast = ast.parse(previous_field.read_bytes()), ast.parse(field_path.read_bytes())
    current_ast.body = [n for n in current_ast.body if not isinstance(n, ast.ClassDef) or n.name != 'CartesianPhysicalFieldAxes']
    check('executed_field_numerics_unchanged_when_adding_lazy_axes', dump(previous_ast) == dump(current_ast))
    for component in components:
        for source in component['source_bindings']:
            if source['path'] == str(field_path):
                check('executed_field_source_bound_to_exact_preserved_bytes', file_binding(previous_field)['sha256'] == source['sha256'])
            else:
                authenticate(source)
    compatibility = dict(kind='same-physical-scalar-new-execution',
        before=retained['source_compatibility']['before'],
        after=now[str(ROOT/'src/substrate_framework/material_joint_response.py')])
    return retained_binding, retained, component_bindings, components, current, compatibility, snapshots


def digest_values(values):
    return [hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest() for value in values]


def load_previous_edge_class(binding):
    path = authenticate(binding)
    name = 'substrate_framework._p254_previous_partitioned_q_response'
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module.PartitionedQResponse


def execute(args, budget):
    checks, timings = [], []
    def check(name, passed, **details):
        checks.append(dict(name=name, passed=bool(passed), **details))
        if not passed:
            raise AssertionError(name)
    retained_binding, retained, component_bindings, components, sources, compatibility, snapshots = retained_sources(args, check)
    initialization, saved = selected_initialization(args)
    history_binding, _ = checked(args.history, args.history_sha256)
    context = bind_saved_model(load_case(args.case), saved, compatibility)
    preparation = preparation_arrays(context)
    original_preparation = bound_arrays(saved['preparation'])
    check('canonical_model_and_every_original_preparation_array_preserved',
          context.identity == saved['model_identity'] and context.specification == saved['specification']
          and set(preparation) == set(original_preparation)
          and all(np.array_equal(v, original_preparation[k]) for k, v in preparation.items()))
    rejected = False
    wrong = deepcopy(compatibility); wrong['before']['sha256'] = '0'*64
    try:
        bind_saved_model(load_case(args.case), saved, wrong)
    except ValueError:
        rejected = True
    check('wrong_source_continuation_remains_rejected', rejected)
    history = context.make_history(history_binding, budget, block_size=8)
    check('genuine_history_starts_at_original_initialization', file_binding(history.initial.path) == saved['joint_boundary'])
    response = context.bind_response(history, np.zeros(3))
    link = response.edges[args.edge_index]
    history.set_interval(0., history.final.time)
    psi = history.frame(link.source, 0.)
    target = history.frame(link.target, 0.)
    if psi.shape != target.shape:
        raise ValueError('this whole-edge smoke selects original Gamma only')
    ket = np.array(psi, order='F', copy=True)*(.013+.017j)
    ket[::3] *= .8+.2j
    bra = np.array(psi.conj().T, order='F', copy=True)*(.019-.011j)
    bra[:, 1::4] *= -.2+.3j
    coordinate = np.array([.03+.02j, -.01+.04j, .02-.03j])
    x = SimpleNamespace(ket={args.edge_index: ket}, bra={args.edge_index: bra},
        u=np.arange(1, response.ndof+1)*(.0001+.0002j),
        p=np.arange(1, response.ndof+1)*(-.0002+.0001j),
        coordinate=coordinate, momentum=np.roll(coordinate, 1)*(.7-.3j))
    inputs = (psi, target, ket, bra, x.u, x.p, x.coordinate, x.momentum,
              response.sectors[link.source].d, response.sectors[link.target].d)
    original_bytes = digest_values(inputs)
    previous_class = load_previous_edge_class(next(b for b in snapshots if 'partitioned_q_response' in b['path']))
    old = object.__new__(previous_class); current = object.__new__(PartitionedQResponse)
    for executor in (old, current):
        executor.r, executor.edge_stage, executor.dynamics = response, history, 'screened-joint'
    for dual in (False, True):
        budget.check_processes(os.getpid(), [])
        label = 'dual' if dual else 'forward'
        started = perf_counter(); expected = old.independent_edge(args.edge_index, 0., x, dual); old_seconds = perf_counter()-started
        budget.check_processes(os.getpid(), [])
        started = perf_counter(); actual = current.independent_edge(args.edge_index, 0., x, dual); new_seconds = perf_counter()-started
        budget.check_processes(os.getpid(), [])
        check(label+'_same_physical_output_components', set(expected) == set(actual))
        for name in expected:
            a, b = actual[name], expected[name]
            error = float(np.max(np.abs(a-b), initial=0.))
            check(label+'_actual_full_edge_'+name, a.shape == b.shape and a.dtype == b.dtype
                  and np.isfinite(a).all() and error <= args.absolute_target,
                  max_abs_difference=error, absolute_target=args.absolute_target)
        # A skipped physical term cannot hide behind a zero/unexcited test state.
        zero_error = float(np.max(np.abs(expected['ket']), initial=0.))
        check(label+'_zero_output_mutation_rejected', zero_error > args.absolute_target,
              zero_mutation_max_abs_difference=zero_error)
        timings.append(dict(direction=label, original_edge_seconds=old_seconds,
                            paired_edge_seconds=new_seconds, original_PW_rows=psi.shape[0], original_columns=psi.shape[1]))
        del expected, actual
        (args.output/'progress.json').write_text(json.dumps(dict(status='partial-until-both-directions',
            source_bindings=sources, checks=checks, timings=timings), indent=2)+'\n')
        print('Actual full-column '+label+' comparison passed.', flush=True)
    check('all_original_material_and_probe_bytes_preserved', original_bytes == digest_values(inputs))
    return dict(schema='P254-joint-partitioned-q-oracle-v1', status='passed',
        extension='physical-global-field-readouts', source_bindings=sources,
        source_compatibility=compatibility, retained_proof=retained_binding,
        retained_numerical_checks_not_rerun=retained['retained_numerical_checks_not_rerun'],
        retained_mapping_checks_not_rerun=retained['new_check_count'],
        retained_history_checks_not_rerun=retained['retained_history_ownership_checks_not_rerun'],
        component_proofs=component_bindings, component_checks_not_rerun=[len(v['checks']) for v in components],
        proof_sources=[file_binding(Path(__file__).resolve()), *snapshots,
                      file_binding(HERE/'source-history/paired_field_moments.pre-native-axes-001.py')],
        checks=checks, new_check_count=len(checks), initialization=initialization,
        history=history_binding, model_identity=context.identity, explicit_history_sectors=len(history.bases),
        selected_edge=args.edge_index, full_PW_rows=psi.shape[0], all_original_columns=psi.shape[1],
        q=[0., 0., 0.], stage_time=0., timings=timings, material_result=False,
        numerical_accuracy_certified=False,
        scope='Actual complete moving-history intake and selected original full-PW/all-column nonzero complex forward/dual consumer comparison. Primitive proofs retained, not replayed. Old source-bound edge algorithm is loaded from its exact preserved bytes. No other original q/RK4/storage code or physical-stage arithmetic changes. Timings are this selected edge only; no complete q/spatial/accuracy/ranking result.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    for name in ('initialization', 'history', 'retained-oracle', 'matter-oracle', 'field-oracle'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--edge-index', type=int, default=0)
    parser.add_argument('--absolute-target', type=float, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=64, max_seconds=args.max_seconds, max_bytes=2097152,
        max_memory_bytes=args.max_memory_bytes, worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(1048576)
    result = dict(schema='P254-paired-contraction-consumer-v1', status='failed')
    started = perf_counter()
    try:
        result = execute(args, budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(seconds=perf_counter()-started, budget=budget.snapshot(),
        process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result['budget']['time_cap'] = 'Cooperative checks and external hard systemd wall cap; standalone consumer has no worker watchdog.'
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
