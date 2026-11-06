#!/usr/bin/env python3
"""Finite behavioral oracle for committed partitioned selected readout.

Run only after integration, on explicitly allocated CPUs6-10. Reuses the real
moving finite q fixture, original raw frames/D, ordered product sources, mapped
checkpoints and QWork. The independent comparator contracts full rectangular
densities against literal matrices, not the production factor readout. No
material/Native trajectory, additional qualification gate, or propagation is
inferred from this finite calculation.
"""
import argparse
from copy import copy
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from run_joint_response_oracle import FiniteSystem, MatrixAction, pack_legs
from run_joint_q_oracle import mapped
from run_joint_partitioned_q_oracle import DiskFiniteHistory
from partitioned_material_readout import partitioned_selected_readout
from run_joint_q_readout import _execution_snapshot, _execution_summary
from substrate_framework.joint_ordered_sources import (
    JointProbe, JointProductPreparation, linear_combination)
from substrate_framework.joint_response_observables import total_force_probe
from substrate_framework.partitioned_joint_slots import array_identity
from substrate_framework.partitioned_q_response import PartitionedQResponse
from substrate_framework.partitioned_q_storage import QBudget, durable_publish
from substrate_framework.partitioned_q_work import QWork, arrays
from substrate_framework.streamed_joint_response import MappedResponse
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-joint-q-readout-oracle-v1'
SMALL = ('u', 'p', 'coordinate', 'momentum')
TOLERANCE = 2e-11  # Existing finite q algebra/discrete-adjoint tolerance.
MUTANT_TOLERANCE = 2e-8  # Existing finite q mutation convention, not relaxed.


def encode(value):
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=float(value.real), imag=float(value.imag))
    if isinstance(value, np.generic):
        return value.item()
    return value


class Evidence:
    def __init__(self):
        self.checks = []

    def require(self, name, passed, **details):
        self.checks.append(dict(name=name, passed=bool(passed), **details))

    def compare(self, name, actual, expected, *, mutation=False, tolerance=None):
        tolerance = (MUTANT_TOLERANCE if mutation else TOLERANCE) if tolerance is None else tolerance
        a, b = np.asarray(actual), np.asarray(expected)
        finite = bool(np.isfinite(a).all() and np.isfinite(b).all())
        same = a.shape == b.shape
        error = float(np.max(np.abs(a-b), initial=0)) if same and finite else None
        close = same and finite and bool(np.allclose(a, b, atol=tolerance, rtol=tolerance))
        self.require(name, finite and same and (not close if mutation else close),
                     actual=actual, expected=expected, max_abs=error,
                     mutation=mutation, atol=tolerance, rtol=tolerance)

    def refuse(self, name, action):
        try:
            value = action()
        except (ValueError, TypeError) as error:
            self.require(name, True, rejected=True, error=type(error).__name__+': '+str(error))
        else:
            self.require(name, False, rejected=False, returned=value)


def source_bindings():
    """Bind actual imported repository dependencies plus both runtime entrypoints."""
    paths = {HERE/'run_joint_q_readout_oracle.py', HERE/'partitioned_material_readout.py',
             HERE/'run_joint_q_readout.py', HERE/'mixed_material_driver_support.py'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and (path.is_relative_to(ROOT/'src') or path.parent == HERE):
                paths.add(path)
    return [file_binding(path) for path in sorted(paths)]


def make_probe(system, response, identity, *, matrices=None, small=None):
    if matrices is None:
        matrices = [system.random((system.sizes[e.target], system.sizes[e.source]))
                    for e in response.edges]
    if small is None:
        small = [system.random((size,)) for size in
                 (response.ndof, response.ndof, response.nfield, response.nfield)]
    operators = [MatrixAction(matrix) for matrix in matrices]
    probe = JointProbe(tuple(op.apply for op in operators), tuple(op.adjoint_apply for op in operators),
                       *small, identity, 'finite complex observable')
    return probe, matrices


def operator_binding(probe, matrices):
    return dict(construction='literal finite rectangular observable matrices; original units',
                matrices=[array_identity(matrix) for matrix in matrices],
                small={key: array_identity(getattr(probe, key)) for key in SMALL})


def full_density_reference(response, stage, state, probe, matrices):
    """Independent scalar convention: sum w Tr(P† X) + c† y, no D here."""
    densities = [state.ket[i] @ stage.frames[edge.source].conj().T
                 + stage.frames[edge.target] @ state.bra[i]
                 for i, edge in enumerate(response.edges)]
    canonical = {key: sum((complex(c).conjugate()*complex(y)
                          for c, y in zip(getattr(probe, key), getattr(state, key))), 0j)
                 for key in SMALL}
    unweighted = [complex(np.trace(matrix.conj().T @ density))
                  for matrix, density in zip(matrices, densities)]
    electronic = [response.sectors[edge.source].weight*value
                  for edge, value in zip(response.edges, unweighted)]
    total = sum(canonical.values(), 0j)
    for term in electronic:
        total += term
    return dict(total=complex(total), canonical=canonical, edge_values=electronic,
                unweighted=unweighted, densities=densities)


def executor_for(system, response, root, budget, cpus, sources, label):
    executor = PartitionedQResponse(response, DiskFiniteHistory(system, response, budget),
        directory=root/label, kernel_adjoint=lambda n, **kw: system.k*n,
        grid_measure=1., budget=budget, identity='finite-current-readout-consumer',
        cpus=cpus, work_directory=root/'stage-work')
    executor.readout_source_bindings = sources
    return executor


def owned_receipts(work):
    key, directory, count = work.last
    records = work.collect(directory, key, count)
    return dict(request=file_binding(directory/'request.json'), key=key,
                directory=str(directory), records=records,
                complete_bindings=[file_binding(directory/('item-%04d' % i)/'complete.json')
                                   for i in range(count)])


def immutable_work_bindings(work):
    # Replay must retain the exact scientific outputs AND original worker costs,
    # rather than merely returning the same number after recomputation.
    return {str(path.relative_to(work.directory)): file_binding(path)
            for path in sorted(work.directory.rglob('*')) if path.is_file()}


def publish_state(response, tangent, path, time, budget):
    writable = mapped(response, path, budget, pack_legs(tangent))
    binding = writable.publish(time=time, identity='finite-historical-committed-producer',
                               budget=budget, kind='forward',
                               provenance=dict(finite_only=True, construction='original ordered product source'))
    state, record = MappedResponse.restore(response, binding, budget,
                                          'finite-historical-committed-producer')
    return state, binding, record


def exercise(gamma, root, budget, cpus, sources, evidence, findings):
    root.mkdir(parents=True, exist_ok=False)
    system = FiniteSystem()
    response = system.response(gamma)
    tag = 'Gamma' if gamma else 'nonGamma'
    if gamma:
        # Each Gamma edge preserves its own sector. Unequal original declared
        # weights expose source-index/weight permutations as well as omission.
        response.sectors = tuple(replace(sector, weight=weight) for sector, weight in
                                 zip(response.sectors, (.19, .31, .47)))
    time = .0173
    stage = system.stage(time)
    evidence.require(tag+'/nonorthogonal_original_frames',
        max(np.linalg.norm(p.conj().T @ p-np.eye(p.shape[1])) for p in system.frames) > 1e-4)
    ranks = [int(np.linalg.matrix_rank(d)) for d in system.ds]
    evidence.require(tag+'/singular_complex_original_D',
        any(rank < d.shape[0] for rank, d in zip(ranks, system.ds))
        and max(np.max(np.abs(d.imag)) for d in system.ds) > 1e-4,
        ranks=ranks, original_columns=list(system.ranks))
    evidence.require(tag+'/nonunit_weights', all(s.weight != 1. for s in response.sectors),
                     weights=[s.weight for s in response.sectors])
    evidence.require(tag+'/genuine_moving_frames',
        max(np.linalg.norm(a-b) for a, b in zip(stage.frames, system.frames)) > 1e-4)
    right, right_matrices = make_probe(system, response, tag+'-ordered-right')
    probe, matrices = make_probe(system, response, tag+'-selected-complex')
    preparation = JointProductPreparation(response, system.stage(0),
        ionic_ordered=system.wi, ionic_minus_ordered=system.wi,
        field_ordered=system.wf if gamma else np.zeros((0, 0)),
        electronic_identity='original-raw-singular-D', ionic_identity='original-oscillator-W',
        field_identity='original-global-field-W', product_identity='original-product')
    executor = executor_for(system, response, root, budget, cpus, sources, 'executor')
    work = QWork(root/'readout-work', budget, cpus, 'finite-current-scalar-work')
    findings.update(q=response.q_fractional, time=time, scalar_reference=[],
                    qwork_receipts=[], mutant_sensitivity=[], replay=[], source_refusals=[],
                    finite_state_construction='AB/BA original product-source legs at saved moving-frame time; no propagation claim')
    states = []
    for ordering in ('AB', 'BA'):
        ordered = preparation.source(right, ordering=ordering)
        # Independent ordered full-density formulas verify a genuine D-bearing
        # source, not an arbitrary matrix echoed through a worker.
        for index, edge in enumerate(response.edges):
            ps, pt = system.frames[edge.source], system.frames[edge.target]
            rs = ps @ system.ds[edge.source] @ ps.conj().T
            rt = pt @ system.ds[edge.target] @ pt.conj().T
            b = right_matrices[index]
            weight = response.ncells*response.sectors[edge.source].weight
            expected = weight*((np.eye(len(rt))-rt) @ b @ rs if ordering == 'AB'
                               else rt @ b @ (np.eye(len(rs))-rs))
            actual = ordered.electronic.ket[index] @ ps.conj().T+pt @ ordered.electronic.bra[index]
            evidence.compare(tag+'/'+ordering+'/independent_product_density_%d' % index, actual, expected)
        tangent = linear_combination((.7+.2j, ordered.electronic),
                                     (-.4+.3j, ordered.ionic), (.6-.5j, ordered.field))
        state, binding, record = publish_state(response, tangent, root/ordering, time, budget)
        states.append((state, binding, record))
        protected = [binding, *record['legs']]
        original_inputs = [file_binding(item['path']) for item in protected]
        evidence.require(tag+'/'+ordering+'/readonly_original_leg_mapping',
                         all(not value.flags.writeable for value in state.ket+state.bra))
        reference = full_density_reference(response, stage, state, probe, matrices)
        executor.readout_operator_binding = operator_binding(probe, matrices)
        check_execution_costs = gamma and ordering == 'AB'
        journals_before = _execution_snapshot(work) if check_execution_costs else None
        actual = partitioned_selected_readout(executor, state, probe, time, work, state_binding=binding)
        evidence.compare(tag+'/'+ordering+'/full_density_trace_and_canonical_scalar', actual, reference['total'])
        evidence.require(tag+'/'+ordering+'/nontrivial_complex_electronic',
                         abs(sum(reference['edge_values']).imag) > 1e-6)
        for key in SMALL if gamma else ('u', 'p'):
            value = reference['canonical'][key]
            evidence.require(tag+'/'+ordering+'/complex_nonzero_'+key, abs(value.real) > 1e-6 and abs(value.imag) > 1e-6,
                             scalar=value)
        receipt = owned_receipts(work)
        if check_execution_costs:
            allocation = _execution_summary(work, journals_before)
            evidence.require('current_finite_readout_costs_are_actual_new_work',
                allocation['current_items'] == set(range(len(response.edges)))
                and allocation['current_completed_worker_cpu_seconds'] > 0,
                allocated_indices=sorted(allocation['current_items']),
                current_worker_cpu=allocation['current_completed_worker_cpu_seconds'])
        for index, row in enumerate(receipt['records']):
            values = arrays(row)
            valid = set(values) == {'value'} and values['value'].shape == (1,) and values['value'].dtype == np.complex128
            evidence.require(tag+'/'+ordering+'/scalar_only_receipt_%d' % index, valid,
                             arrays=row['arrays'], worker_cpu=row['cpu'], worker_pid=row['pid'])
            if valid:
                evidence.compare(tag+'/'+ordering+'/edge_full_density_%d' % index,
                                 values['value'][0], reference['edge_values'][index])
            evidence.require(tag+'/'+ordering+'/actual_pinned_worker_cost_%d' % index,
                row['pid'] != os.getpid() and row['cpu'] in cpus and row['affinity'] == [row['cpu']]
                and row['blas_threads'] == 1 and row['cpu_seconds'] > 0
                and row['elapsed_seconds'] > 0
                and row['process_peak_rss_bytes'] > 0, receipt=row)
        findings['qwork_receipts'].append(dict(ordering=ordering, **receipt))
        findings['scalar_reference'].append(dict(ordering=ordering, checkpoint=binding,
            producer=record['identity'], consumer=executor.identity, actual=actual,
            expected=reference['total'], canonical=reference['canonical'], edge_values=reference['edge_values'],
            density_shapes=[list(x.shape) for x in reference['densities']]))
        # Fresh objects reuse actual content-bound receipts and their original
        # immutable costs. No echo callback or invocation-count surrogate.
        before = immutable_work_bindings(work)
        fresh = executor_for(system, response, root, budget, cpus, sources, 'replay-'+ordering)
        fresh.readout_operator_binding = operator_binding(probe, matrices)
        replay_work = QWork(work.directory, budget, cpus, work.identity)
        replay_journals_before = _execution_snapshot(replay_work) if check_execution_costs else None
        replay = partitioned_selected_readout(fresh, state, probe, time, replay_work, state_binding=binding)
        evidence.compare(tag+'/'+ordering+'/fresh_object_replay_scalar', replay, actual, tolerance=0.)
        if check_execution_costs:
            replay_costs = _execution_summary(replay_work, replay_journals_before)
            evidence.require('cached_replay_keeps_historical_costs_out_of_current_work',
                replay_costs['current_items'] == set()
                and replay_costs['current_completed_worker_cpu_seconds'] == 0
                and replay_costs['receipt_worker_cpu_seconds'] == allocation['receipt_worker_cpu_seconds'],
                current_worker_cpu=replay_costs['current_completed_worker_cpu_seconds'],
                retained_worker_cpu=replay_costs['receipt_worker_cpu_seconds'])
        after = immutable_work_bindings(replay_work)
        evidence.require(tag+'/'+ordering+'/replay_preserves_actual_outputs_and_costs', before == after,
                         prior_bindings=before, replay_bindings=after)
        findings['replay'].append(dict(ordering=ordering, value=replay, immutable_bindings=after))
        evidence.require(tag+'/'+ordering+'/original_checkpoint_and_legs_unchanged',
                         [file_binding(item['path']) for item in protected] == original_inputs)
        # Numerically plausible wrong contractions, not tests of source spelling.
        small = sum(reference['canonical'].values(), 0j)
        densities = reference['densities']
        mutants = dict(
            omitted_weight=small+sum(reference['unweighted']),
            duplicated_weight=small+sum(response.sectors[e.source].weight*v
                                        for e, v in zip(response.edges, reference['edge_values'])),
            extra_cell_factor=small+response.ncells*sum(reference['edge_values']),
            real_projection=complex(reference['total'].real),
            unconjugated_probe=small+sum(response.sectors[e.source].weight*np.trace(p.T @ x)
                                        for e, p, x in zip(response.edges, matrices, densities)),
            conjugated_bra=small+sum(response.sectors[e.source].weight*np.trace(p.conj().T @
                (state.ket[i] @ stage.frames[e.source].conj().T+stage.frames[e.target] @ state.bra[i].conj()))
                for i, (e, p) in enumerate(zip(response.edges, matrices))),
            unconjugated_canonical=sum(sum(c*y for c, y in zip(getattr(probe, key), getattr(state, key)))
                                      for key in SMALL)+sum(reference['edge_values']),
            extra_original_D=small+sum(response.sectors[e.source].weight*np.trace(p.conj().T @
                (state.ket[i] @ response.sectors[e.source].d @ stage.frames[e.source].conj().T
                 +stage.frames[e.target] @ response.sectors[e.target].d @ state.bra[i]))
                for i, (e, p) in enumerate(zip(response.edges, matrices))),
            initial_not_saved_frames=full_density_reference(response, system.stage(0), state, probe, matrices)['total'])
        if gamma:
            permuted = [densities[2], densities[1], densities[0]]
            mutants['wrong_source_order'] = small+sum(response.sectors[e.source].weight*np.trace(p.conj().T @ x)
                                                     for e, p, x in zip(response.edges, matrices, permuted))
        for name, value in mutants.items():
            evidence.compare(tag+'/'+ordering+'/sensitive_'+name, value, reference['total'], mutation=True)
        findings['mutant_sensitivity'].append(dict(ordering=ordering, expected=reference['total'], mutants=mutants))
    first, binding, metadata = states[0]
    second, other_binding, _ = states[1]
    refs = [row['expected'] for row in findings['scalar_reference']]
    evidence.compare(tag+'/ordered_AB_not_BA', refs[0], refs[1], mutation=True)
    # Keep probe's declared identity/small coefficients unchanged, but bind and
    # evaluate genuinely changed literal operators: stale source cache must fail.
    changed = [matrix*(.73+.21j) for matrix in matrices]
    changed_probe, _ = make_probe(system, response, probe.identity, matrices=changed,
                                 small=[getattr(probe, key) for key in SMALL])
    executor.readout_operator_binding = operator_binding(changed_probe, changed)
    changed_reference = full_density_reference(response, stage, first, changed_probe, changed)
    changed_value = partitioned_selected_readout(executor, first, changed_probe, time, work, state_binding=binding)
    evidence.compare(tag+'/changed_bound_operator_correct_scalar', changed_value, changed_reference['total'])
    evidence.compare(tag+'/changed_bound_operator_not_stale_cache', changed_value, refs[0], mutation=True)
    findings['operator_rebinding'] = dict(actual=changed_value, expected=changed_reference['total'],
                                         prior=refs[0], binding=executor.readout_operator_binding,
                                         receipts=owned_receipts(work))
    executor.readout_operator_binding = operator_binding(probe, matrices)
    direction = np.array([1., 0., 0.])
    force = total_force_probe(response, direction, identity=tag+'-selected-total-force')
    force_matrices = [-np.einsum('a,aij->ij', direction, e.coupling.g) for e in response.edges]
    executor.readout_operator_binding = operator_binding(force, force_matrices)
    force_value = partitioned_selected_readout(executor, first, force, time, work, state_binding=binding)
    grad = sum(response.sectors[e.source].weight*np.array([np.trace(g.conj().T @ x) for g in e.coupling.g])
               for e, x in zip(response.edges, full_density_reference(response, stage, first, probe, matrices)['densities']))
    expected_force = -response.curvature @ first.u-grad
    if gamma:
        expected_force -= response.sqrt_n*response.mixed @ first.coordinate
    evidence.compare(tag+'/actual_selected_total_force_probe', force_value, np.vdot(direction, expected_force))
    findings['selected_total_force'] = dict(actual=force_value, expected=np.vdot(direction, expected_force),
                                           receipts=owned_receipts(work))
    executor.readout_operator_binding = operator_binding(probe, matrices)
    def consume(state=first, source=binding, at=time):
        return partitioned_selected_readout(executor, state, probe, at, work, state_binding=source)
    start = len(evidence.checks)
    evidence.refuse(tag+'/refuse_wrong_checkpoint_digest', lambda: consume(source=dict(binding, sha256='0'*64)))
    evidence.refuse(tag+'/refuse_wrong_checkpoint_bytes', lambda: consume(source=dict(binding, bytes=binding['bytes']+1)))
    evidence.refuse(tag+'/refuse_different_authentic_committed_state', lambda: consume(source=other_binding))
    evidence.refuse(tag+'/refuse_wrong_endpoint_time', lambda: consume(at=time+.001))
    mismatched = copy(first)
    mismatched.u = first.u+(.13+.07j)
    evidence.refuse(tag+'/refuse_canonical_checkpoint_mismatch', lambda: consume(state=mismatched))
    for name, patch in (('q', dict(q=[.17, .23, .31])), ('complete', dict(complete=False)),
                        ('kind', dict(kind='dual'))):
        invalid = durable_publish(root/('wrong-'+name+'.json'), dict(metadata, **patch), budget)
        evidence.refuse(tag+'/refuse_authenticated_wrong_'+name, lambda invalid=invalid: consume(source=invalid))
    executor.readout_source_bindings = [dict(sources[0], sha256='0'*64), *sources[1:]]
    try:
        evidence.refuse(tag+'/refuse_current_source_digest_mismatch', consume)
    finally:
        executor.readout_source_bindings = sources
    executor.readout_source_bindings = []
    try:
        evidence.refuse(tag+'/refuse_missing_current_source_authority', consume)
    finally:
        executor.readout_source_bindings = sources
    executor.readout_operator_binding = {}
    try:
        evidence.refuse(tag+'/refuse_missing_operator_construction_authority', consume)
    finally:
        executor.readout_operator_binding = operator_binding(probe, matrices)
    findings['source_refusals'] = evidence.checks[start:]
    # A refusal cannot replace/poison already valid scalar receipts.
    recovery = consume()
    evidence.compare(tag+'/valid_source_replay_after_refusals', recovery,
                     findings['scalar_reference'][0]['actual'], tolerance=0.)
    findings['after_refusals'] = dict(value=recovery, receipts=owned_receipts(work))
    for index, (state, committed, record) in enumerate(states):
        try:
            for original in (committed, *record['legs']):
                authenticate(original)
            evidence.require(tag+'/all_controls_preserve_original_state_bytes_%d' % index, True)
        except (ValueError, OSError) as error:
            evidence.require(tag+'/all_controls_preserve_original_state_bytes_%d' % index,
                             False, error=str(error))
        for key in SMALL:
            expected = np.asarray([complex(*z) for z in record['small'][key]])
            evidence.compare(tag+'/all_controls_preserve_canonical_%d_%s' % (index, key),
                             getattr(state, key), expected, tolerance=0.)
    return findings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', type=int, nargs='+', choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, default=1800.)
    parser.add_argument('--max-actions', type=int, default=1024)
    parser.add_argument('--max-output-bytes', type=int, default=268435456)
    parser.add_argument('--max-memory-bytes', type=int, default=17179869184)
    parser.add_argument('--worker-memory-bytes', type=int, default=4294967296)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('oracle output must be new')
    args.work.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, {args.cpus[0]})
    budget = QBudget(max_actions=args.max_actions, max_seconds=args.max_seconds,
                     max_bytes=args.max_output_bytes, max_memory_bytes=args.max_memory_bytes,
                     worker_memory_bytes=args.worker_memory_bytes)
    evidence = Evidence()
    sources = source_bindings()
    report = dict(schema=SCHEMA, status='failed', source_bindings=sources,
        checks=evidence.checks, findings={}, cpus=args.cpus, finite_only=True,
        material_result=False, tolerance=TOLERANCE, mutant_tolerance=MUTANT_TOLERANCE,
        reference='independent sum of weighted full-density Tr(P dagger X) and explicit conjugate canonical coefficients',
        limitations=['Finite Gamma and directed nonGamma three-sector fixture only; no material/Native numerical claim.',
                     'Original product-source legs are read at saved moving frames; no new trajectory is propagated.',
                     'No claim of full-grid/64q, finite amplitude, energy, ranking, or whole-objective closure.',
                     'No interrupted in-flight-worker kill is injected; genuine completed QWork replay and source refusals are exercised.'])
    for gamma, tag in ((True, 'Gamma'), (False, 'nonGamma')):
        report['findings'][tag] = {}
        try:
            exercise(gamma, args.work/tag, budget, args.cpus, sources, evidence, report['findings'][tag])
        except Exception as error:
            evidence.require(tag+'/execution_completed', False,
                             error=type(error).__name__+': '+str(error), traceback=traceback.format_exc())
    try:
        for binding in sources:
            authenticate(binding)
        evidence.require('all_bound_implementation_and_fixture_sources_unchanged', True)
    except (ValueError, OSError) as error:
        evidence.require('all_bound_implementation_and_fixture_sources_unchanged', False, error=str(error))
    report['budget'] = budget.snapshot()
    report['check_count'] = len(evidence.checks)
    if report['findings'].keys() == {'Gamma', 'nonGamma'} and all(row['passed'] for row in evidence.checks):
        report['status'] = 'passed'
    with args.output.open('x') as stream:
        json.dump(encode(report), stream, indent=2, allow_nan=False)
        stream.write('\n')
    if report['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
