#!/usr/bin/env python3
"""Exercise context-only q pairs against independent finite matrix/Fock oracles."""
import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import sympy as sp

from run_joint_response_oracle import (FiniteSystem, DenseTangent, FockOracle,
    check, reject, CHECKS, pack_legs, unpack_legs)
from run_joint_partitioned_q_oracle import DiskFiniteHistory
from run_joint_q_material import encode
from substrate_framework.partitioned_q_response import PartitionedQResponse
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.source_local_q_response import SourceLocalQResponse, prepare_pullback_frames
from substrate_framework.source_local_q_covariance import SourceLocalTwoSidedCovariance
from substrate_framework.source_local_q_state import next_coefficients
from substrate_framework.source_local_q_plan import source_local_plan
from substrate_framework.streamed_joint_covariance import initial_contraction
from substrate_framework.source_local_q_edge import dual_moments, dual_free
from substrate_framework.streamed_material_state import file_binding


def symbolic_reverse():
    h = sp.symbols('h', real=True)
    a = [sp.Matrix([[1+sp.I, 2], [-1, 0]]),
         sp.Matrix([[0, -sp.I], [2, 1]]), sp.Matrix([[2, 1+sp.I], [0, -1]])]
    eye, seed = sp.eye(2), sp.Matrix([1+sp.I, 2-sp.I])
    k0 = a[0]
    k1 = a[1]*(eye+h*k0/2)
    k2 = a[1]*(eye+h*k1/2)
    k3 = a[2]*(eye+h*k2)
    forward = eye+h*(k0+2*k1+2*k2+k3)/6
    current, out = h*seed/6, seed
    bad_current, bad_out = current, seed
    for stage, matrix in enumerate((a[2], a[1], a[1], a[0])):
        derivative, bad_derivative = matrix.H*current, matrix.H*bad_current
        out, bad_out = out+derivative, bad_out+bad_derivative
        if stage < 3:
            x, y = next_coefficients(stage, h)
            current = x*seed+y*derivative
            bad_current = x*seed+(h/2 if stage == 0 else y)*bad_derivative
    if any(sp.simplify(z) != 0 for z in out-forward.H*seed):
        raise AssertionError('noncommuting exact reverse graph differs')
    CHECKS.append(dict(name='symbolic_noncommuting_reverse_graph', passed=True, exact=True))
    witness = np.asarray((bad_out-forward.H*seed).subs(h, sp.Rational(1,5)), complex)
    check('wrong_first_reverse_chain_coefficient_sensitive', witness, np.zeros((2,1)), mutation=True)


def generator(response, system, time):
    size = len(pack_legs(response.zero_tangent()))
    return np.column_stack([pack_legs(response.rhs(system.stage(time), unpack_legs(response, column)))
                            for column in np.eye(size, dtype=complex)])


def forward_matrix(response, system, time, dt):
    a, b, c = (generator(response, system, t) for t in (time, time+dt/2, time+dt))
    eye = np.eye(len(a), dtype=complex)
    k0 = a
    k1 = b @ (eye+dt*k0/2)
    k2 = b @ (eye+dt*k1/2)
    k3 = c @ (eye+dt*k2)
    return eye+dt*(k0+2*k1+2*k2+k3)/6


def terminal(response, system, probe, time):
    stage = system.stage(time)
    ket, bra = [], []
    for i, edge in enumerate(response.edges):
        weight = response.sectors[edge.source].weight
        ket.append(weight*probe.actions[i](stage.frames[edge.source]))
        bra.append(weight*probe.adjoints[i](stage.frames[edge.target]).conj().T)
    return response.tangent(ket, bra, probe.u, probe.p, probe.coordinate, probe.momentum)


def exact_dual(response, system, probe, times):
    value = pack_legs(terminal(response, system, probe, times[-1]))
    for a, b in reversed(list(zip(times[:-1], times[1:]))):
        value = forward_matrix(response, system, a, b-a).conj().T @ value
    return unpack_legs(response, value)


def split_edge_checks(system, response, executor, probe, budget, tag):
    """Compare both consumed projections to the unchanged complete generator."""
    state = terminal(response, system, probe, 0.)
    executor.prepare_time(0.)
    for index in range(len(response.edges)):
        reference = executor.independent_edge(index, 0., state, True)
        before = budget.snapshot()['hamiltonian_actions']
        moments = dual_moments(executor, index, 0., state)
        check(tag+'_moment_only_no_discarded_H0_'+str(index),
              budget.snapshot()['hamiltonian_actions']-before, 0, tolerance=0.)
        before = budget.snapshot()['hamiltonian_actions']
        free = dual_free(executor, index, 0., state)
        check(tag+'_free_only_original_two_H0_'+str(index),
              budget.snapshot()['hamiltonian_actions']-before, 2, tolerance=0.)
        for name in ('n', 'g', 'ell'):
            check(tag+'_moment_only_'+name+'_'+str(index),
                  moments[name], reference[name], tolerance=2e-11)
        for name, value in zip(('ket', 'bra'), free):
            check(tag+'_free_only_'+name+'_'+str(index),
                  value, reference[name], tolerance=2e-11)


def symbolic_split_cost():
    count = sp.symbols('count', integer=True, positive=True)
    previous = sp.symbols('previous', integer=True, nonnegative=True)
    free = sp.summation(16*previous+6, (previous, 0, count-1))+4*count
    if sp.simplify(free-(8*count**2+2*count)) != 0:
        raise AssertionError('exact free-only replay count differs')
    if sp.simplify(free+1-(8*count**2+2*count)) == 0:
        raise AssertionError('one discarded H0 operation mutation was not detected')
    CHECKS.append(dict(name='symbolic_split_replay_cost', passed=True, exact=True,
        free_generator_calls=str(free), moment_only_calls=str(4*count),
        one_extra_generator_sensitive=True))


def contraction_restart_check(preparation, states, executor, directory, grids, reference):
    from substrate_framework import source_local_q_contraction as contraction
    original, cpus = contraction._electronic_source, executor.edge_work.cpus
    last = len(executor.r.edges)-1
    def interrupted(*args):
        value = original(*args)
        if args[-1] == last:
            raise RuntimeError('injected interruption after charged final-source computation')
        return value
    for state in states:
        state.release()
    executor.edge_work.cpus = (cpus[0],)
    contraction._electronic_source = interrupted
    try:
        try:
            contraction.source_contraction(preparation, *states, executor, directory)
        except RuntimeError:
            pass
        else:
            raise AssertionError('interrupted final source was published')
    finally:
        contraction._electronic_source = original
        executor.edge_work.cpus = cpus
    prefix = {path: path.read_bytes() for path in directory.glob('*/generation-*/item-*/complete.json')}
    check('contraction_interruption_preserves_completed_prefix', len(prefix), last, tolerance=0.)
    before = executor.budget.snapshot()['hamiltonian_actions']
    recovered = contraction.source_contraction(preparation, *states, executor, directory)
    check('contraction_resume_charges_only_uncommitted_source',
        executor.budget.snapshot()['hamiltonian_actions']-before,
        8*sum(len(grid)-1 for grid in grids), tolerance=0.)
    if any(path.read_bytes() != raw for path, raw in prefix.items()):
        raise AssertionError('contraction resume rewrote completed source receipts')
    for ordering in ('AB', 'BA'):
        check('contraction_resume_original_'+ordering,
              recovered[ordering]['total'], reference[ordering]['total'], tolerance=2e-11)


def exercise(system, gamma, unequal, directory, budget, cpus):
    r = system.response(gamma)
    tag = ('Gamma' if gamma else 'directed-rectangular-q')+('-unequal-weights' if unequal else '')
    if unequal:
        r.sectors = tuple(replace(sector, weight=weight)
                          for sector, weight in zip(r.sectors, (.2, .5, .3)))
    dense = DenseTangent(system, r)
    left = dense.probe(system.random((dense.size,)), tag+'-A')
    right = dense.probe(system.random((dense.size,)), tag+'-B')
    history = DiskFiniteHistory(system, r, budget)
    executor = PartitionedQResponse(r, history, kernel_adjoint=lambda n, **kw: system.k*n,
        grid_measure=1., budget=budget, directory=directory/'executor', identity=tag,
        cpus=cpus, work_directory=directory/'physical-work')
    split_edge_checks(system, r, executor, left, budget, tag)
    from substrate_framework.joint_ordered_sources import JointProductPreparation
    preparation = JointProductPreparation(r, system.stage(0), ionic_ordered=system.wi,
        ionic_minus_ordered=system.wi, field_ordered=system.wf if gamma else np.zeros((0,0)),
        electronic_identity='original-raw-singular-D', ionic_identity='original-ionic-W',
        field_identity='original-global-field-W', product_identity='original-product')
    grids = (np.array([0., .003, .009]), np.array([0., .004]))
    prepare_pullback_frames(executor, grids)
    states, expected = [], []
    for side, probe, grid in zip(('left', 'right'), (left, right), grids):
        state = SourceLocalQResponse(executor, directory/'pair'/side).pullback(probe, grid)
        reference = exact_dual(r, system, probe, grid)
        check(tag+'_'+side+'_full_rectangular_discrete_dual', pack_legs(state), pack_legs(reference), tolerance=2e-11)
        states.append(state)
        expected.append(reference)
    reference = initial_contraction(preparation, *expected)
    if gamma and not unequal:
        contraction_restart_check(preparation, states, executor,
            directory/'interrupted-contraction', grids, reference)
    covariance = SourceLocalTwoSidedCovariance(preparation, executor)
    actual = covariance.pair(left, right, *grids, directory/'pair')
    for ordering in ('AB', 'BA'):
        for origin in ('electronic', 'ionic', 'field', 'total'):
            check(tag+'_'+ordering+'_'+origin, actual[ordering][origin], reference[ordering][origin], tolerance=2e-11)
    for item in ('commutator', 'symmetrized'):
        check(tag+'_'+item, actual[item], reference[item], tolerance=2e-11)
    check(tag+'_raw_compatibility_preserved',
        [[x['left'], x['right']] for x in actual['raw_frame_compatibility']],
        [[x['left'], x['right']] for x in reference['raw_frame_compatibility']], tolerance=2e-11)
    b0 = terminal(r, system, right, 0.)
    check(tag+'_genuine_second_pullback_sensitive', initial_contraction(preparation, expected[0], b0)['AB']['total'],
          actual['AB']['total'], mutation=True)
    if not unequal:
        fock = FockOracle(system)
        pa, pb = dense.pullback(grids[0][-1], left), dense.pullback(grids[1][-1], right)
        for ordering in ('AB', 'BA'):
            for origin, value in fock.parts(r, pa, pb, ordering).items():
                check(tag+'_independent_Fock_'+ordering+'_'+origin,
                      actual[ordering][origin], value, tolerance=2e-8)
    prepare_pullback_frames(executor, grids)
    restored = SourceLocalQResponse(executor, directory/'pair'/'left').restore(left, grids[0], actual['left_context'])
    check(tag+'_restored_source_context', pack_legs(restored), pack_legs(expected[0]), tolerance=2e-11)
    last = len(r.edges)-1
    check(tag+'_last_edge_after_sequence_end', restored.edge(last)[0],
          expected[0].ket[last], tolerance=2e-11)
    try:
        restored.edge(len(r.edges))
    except IndexError:
        pass
    else:
        raise AssertionError('invalid edge did not raise IndexError')
    check(tag+'_cache_after_invalid_index', restored.edge(last)[1],
          expected[0].bra[last], tolerance=2e-11)
    reject(tag+'_foreign_probe_context_refused', lambda: SourceLocalQResponse(
        executor, directory/'pair'/'left').restore(right, grids[0], actual['left_context']))
    record = json.loads(Path(actual['left_context']['path']).read_text())
    context = record['intervals'][0]['contexts'][0]
    path = Path(context['path'])
    original = path.read_bytes()
    path.write_bytes(original+b' ')
    try:
        reject(tag+'_changed_context_digest_refused', lambda: SourceLocalQResponse(
            executor, directory/'pair'/'left').restore(left, grids[0], actual['left_context']))
    finally:
        path.write_bytes(original)
    payloads = list((directory/'pair').glob('**/*.npy'))
    for path in payloads:
        if path.name == 'potential.npy':
            continue
        value = np.load(path, mmap_mode='r', allow_pickle=False)
        if (path.name not in ('ordered.npy', 'compatibility.npy') or value.shape != (2,)
                or value.nbytes > 32):
            raise AssertionError('source-local propagation retained a nonscalar electronic payload')
    CHECKS.append(dict(name=tag+'_actual_payload_inventory', passed=True,
        payloads=[str(p.relative_to(directory)) for p in payloads],
        result_binding=actual['context_pair_binding']))
    return actual


def field_smoke(system, directory, budget, cpus):
    from substrate_framework.joint_ordered_sources import JointProbe, JointProductPreparation
    r = system.response(True)
    zero = r.zero_tangent()
    actions = tuple(lambda value: np.zeros_like(value) for edge in r.edges)
    probe = JointProbe(actions, actions, zero.u, zero.p, zero.coordinate,
        np.ones(r.nfield), 'original-Gamma-field-momentum', 'finite Maxwell canonical momentum')
    executor = PartitionedQResponse(r, DiskFiniteHistory(system, r, budget),
        kernel_adjoint=lambda n, **kw: system.k*n, grid_measure=1., budget=budget,
        directory=directory/'executor', identity='field-only-Gamma', cpus=cpus,
        work_directory=directory/'physical-work')
    preparation = JointProductPreparation(r, system.stage(0), ionic_ordered=system.wi,
        ionic_minus_ordered=system.wi, field_ordered=system.wf,
        electronic_identity='original-raw-singular-D', ionic_identity='original-ionic-W',
        field_identity='original-global-field-W', product_identity='original-product')
    pair = SourceLocalTwoSidedCovariance(preparation, executor)
    initial = pair.pair(probe, probe, [0.], [0.], directory/'initial')
    check('zero_time_field_only_electronic', initial['AB']['electronic'], 0., tolerance=0.)
    check('zero_time_field_only_ionic', initial['AB']['ionic'], 0., tolerance=0.)
    check('zero_time_field_original_preparation', initial['AB']['field'],
          system.wf[r.nfield, r.nfield], tolerance=0.)
    grids = ([0., .009], [0., .004])
    before = budget.snapshot()
    plan = source_local_plan(r, *grids, len(cpus))
    actual = pair.pair(probe, probe, *grids, directory/'moving')
    after = budget.snapshot()
    check('measured_bare_actions_match_context_replay_plan',
          after['hamiltonian_actions']-before['hamiltonian_actions'],
          plan['bare_hamiltonian_actions'], tolerance=0.)
    check('measured_bare_columns_match_context_replay_plan',
          after['action_columns']-before['action_columns'],
          plan['bare_hamiltonian_columns'], tolerance=0.)
    check('stock_four_stage_action_forecast_is_insufficient',
          2*len(r.edges)*4*2, after['hamiltonian_actions']-before['hamiltonian_actions'],
          mutation=True, tolerance=0.)
    payload = sum(np.load(path, mmap_mode='r', allow_pickle=False).nbytes
                  for path in (directory/'moving').glob('**/*.npy'))
    check('actual_final_context_payload_matches_schedule', payload,
          plan['final_payload_upper_bound'], tolerance=0.)
    CHECKS.append(dict(name='measured_context_only_pair_cost', passed=True,
        hamiltonian_actions=after['hamiltonian_actions']-before['hamiltonian_actions'],
        action_columns=after['action_columns']-before['action_columns'],
        payload_bytes=payload, plan=plan))
    reference = initial_contraction(preparation, *[
        exact_dual(r, system, probe, grid) for grid in grids])
    for ordering in ('AB', 'BA'):
        for origin in ('electronic', 'ionic', 'field', 'total'):
            check('field_only_seed_'+ordering+'_'+origin,
                  actual[ordering][origin], reference[ordering][origin], tolerance=2e-11)
    check('field_only_seed_develops_electron_origin', actual['AB']['electronic'],
          0., mutation=True, tolerance=1e-14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', nargs='+', type=int, default=[6])
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpus[0]})
    args.work.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1000000, max_seconds=300, max_bytes=128000000,
                     max_memory_bytes=4000000000, worker_memory_bytes=1000000000)
    result = dict(schema='P254-source-local-q-oracle-v1', status='failed', checks=CHECKS,
        sources=[file_binding(path) for path in [
            *[ROOT/'src/substrate_framework'/name for name in (
                'source_local_q_state.py', 'source_local_q_context.py',
                'source_local_q_response.py', 'source_local_q_covariance.py', 'source_local_q_plan.py',
                'source_local_q_edge.py', 'source_local_q_contraction.py',
                'partitioned_q_response.py', 'streamed_joint_response.py',
                'streamed_joint_covariance.py', 'partitioned_q_work.py',
                'partitioned_q_storage.py', 'joint_ordered_sources.py')],
            *[HERE/name for name in ('run_source_local_q_oracle.py',
                'run_joint_response_oracle.py', 'run_joint_q_oracle.py',
                'run_joint_partitioned_q_oracle.py')]]],
        scientific_scope='Finite moving examples only; not an actual material-noise result or a production disk-fit claim')
    try:
        symbolic_reverse()
        symbolic_split_cost()
        system = FiniteSystem()
        for gamma, unequal, name in ((True, False, 'Gamma'), (False, False, 'directed'), (True, True, 'weighted')):
            exercise(system, gamma, unequal, args.work/name, budget, args.cpus)
        field_smoke(system, args.work/'field-only', budget, args.cpus)
        result['status'] = 'passed'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    result['measured_file_bytes'] = sum(p.stat().st_size for p in args.work.rglob('*') if p.is_file())
    with args.output.open('x') as stream:
        json.dump(encode(result), stream, indent=2, allow_nan=False)
        stream.write('\n')
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
