#!/usr/bin/env python3
"""Independent dense SciPy/SymPy column-Arnoldi oracle; never material dynamics."""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter

from run_block_krylov_oracle import (
    HERE, ROOT, np, sp, expm, binding, authenticate, explicit, publish,
    symbolic_proof, threadpool_limits)
from substrate_framework.column_frame_krylov import ColumnFrameKrylov, memory_bound_bytes
from substrate_framework.column_krylov_bounds import column_envelope
from substrate_framework.full_frame_krylov import physical_invariants


def source_paths():
    return [HERE/name for name in (
        'run_column_krylov_oracle.py', 'run_column_krylov_material.py',
        'run_block_krylov_oracle.py', 'run_block_krylov_material.py',
        'run_streamed_material_evolution.py')] + [ROOT/'src/substrate_framework'/name for name in (
        'column_frame_krylov.py', 'column_krylov_bounds.py', 'full_frame_krylov.py',
        'material_frozen_action.py', 'qe_material_hamiltonian.py',
        'streamed_material_evolution.py', 'streamed_sector_executor.py',
        'streamed_material_state.py', 'qe_material_dyads.py', 'qe_bloch_ionic_action.py',
        'qe_material_ionic_action.py', 'qe_bloch_qzero.py', 'qe_material_grid.py',
        'qe_bloch_dyads.py')]


def sources():
    current = [binding(path) for path in source_paths()]
    registration = json.loads((HERE/'column-krylov-registration.json').read_text())
    if current != registration['registered_source_bindings']:
        raise ValueError('column preregistration no longer binds exact scientific source bytes')
    return current+[binding(HERE/'column-krylov-registration.json')]


def make(h, frame, *, dimension=None, cap=None, action=None, **options):
    return ColumnFrameKrylov(action or (lambda x: h @ x), frame,
        initial_action=h @ frame, initial_action_identity='dense-frozen-H-and-frame',
        operator_identity='dense-frozen-H-and-frame',
        max_dimension=h.shape[0] if dimension is None else dimension,
        action_cap=h.shape[0]*frame.shape[1] if cap is None else cap,
        max_memory_bytes=64*1024**2,
        ambient_growth_bound=float(np.linalg.norm((h-h.conj().T)/(2j), 2)),
        operator_norm_bound=float(np.linalg.norm(h, 2)), **options)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    source, started = sources(), perf_counter()
    checks, symbolic = symbolic_proof()
    time_symbol, shift = sp.symbols('t mu', real=True)
    hs = sp.Matrix([[1, sp.I], [-sp.I, 2]])
    # Independent exact shift compensation, including a complex nonorthogonal frame.
    raw = sp.Matrix([[2, 1+sp.I], [sp.I, 3]])
    shifted = sp.exp(-sp.I*time_symbol*shift)*(-sp.I*time_symbol*(hs-shift*sp.eye(2))).exp()*raw
    checks['sympy_real_shift_compensation'] = (shifted-(-sp.I*time_symbol*hs).exp()*raw).applyfunc(sp.simplify) == sp.zeros(2)
    rng = np.random.default_rng(225280181)
    a = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
    h = (a+a.conj().T)/4+3*np.eye(6)
    frame = rng.normal(size=(6, 4))+1j*rng.normal(size=(6, 4))
    frame[:, 2] = (.3+.4j)*frame[:, 0]+.7*frame[:, 1]
    frame[:, 3] = 0
    d = np.array([[.8, .1j, .13, .03j], [-.1j, .6, .07j, .02],
                  [.13, -.07j, .4, .04j], [-.03j, .02, -.04j, .2]])
    times = [0., .13, -.21, .37]
    calls = []
    def action(x):
        calls.append(x.shape[1])
        return h @ x
    model = make(h, frame, action=action)
    checks['seed_action_reused'] = calls == []
    arrays = dict(H=h, frame=frame, raw_D=d)
    records, errors = [], []
    for dimension in range(1, 7):
        rows = []
        stage = model.stage()
        for index, time in enumerate(times):
            actual, estimate = model.propagate(time)
            exact = expm(-1j*time*h) @ frame
            error = float(np.linalg.norm(actual-exact))
            checks[f'bound_dim{dimension}_time{index}'] = error <= estimate['absolute_error_envelope_estimate']+2e-10
            checks[f'raw_shape_dim{dimension}_time{index}'] = actual.shape == frame.shape and np.all(actual[:, 3] == 0)
            arrays[f'dim{dimension}_time{index}'] = actual
            arrays[f'exact_time{index}'] = exact
            rows.append(dict(actual_dense_error=error, **estimate))
        records.append(dict(diagnostics=stage['diagnostics'], times=rows))
        errors.append(rows[-1]['actual_dense_error'])
        if dimension < 6:
            model.expand()
    checks['full_dimension_dense_match'] = max(r['actual_dense_error'] for r in records[-1]['times']) < 2e-10
    checks['initial_space_not_falsely_closed'] = errors[0] > 1e-3
    checks['refinement_reduces_error'] = errors[-1] < errors[2] < errors[0]
    checks['batched_action_accounting'] = calls == [3]*5 and model.requested_columns == 15
    checks['no_cross_column_projected_matrix'] = all(t.shape == (6, 6) for t in model.stage()['projected'][:3])
    final, _ = model.propagate(times[-1])
    rho = final @ d @ final.conj().T
    wrong_sign = expm(1j*times[-1]*h) @ frame
    missing = final.copy()
    missing[:, 2] = 0
    normalized = frame.copy()
    normalized[:, :3] /= np.linalg.norm(normalized[:, :3], axis=0)
    no_phase = final*np.exp(1j*times[-1]*model.shifts)[None, :]
    mutations = dict(omit_original_column=float(np.linalg.norm(missing-final)),
        wrong_sign=float(np.linalg.norm(wrong_sign-final)),
        omit_shift_phase=float(np.linalg.norm(no_phase-final)),
        normalize_raw_columns=float(np.linalg.norm(expm(-1j*times[-1]*h) @ normalized-final)),
        diagonalize_D=float(np.linalg.norm(final @ np.diag(np.diag(d)) @ final.conj().T-rho)),
        false_convergence_if_residual_zeroed=errors[0])
    for name, discrepancy in mutations.items():
        checks['mutation_'+name] = discrepancy > 1e-6
    initial_inv, final_inv = physical_invariants(frame, d), physical_invariants(final, d)
    for key in ('number', 'purity'):
        checks['unchanged_raw_D_'+key] = abs(initial_inv[key]['real']-final_inv[key]['real']) < 2e-9
    nh = h+1j*np.diag([.5, -.4, .2, .1, -.3, .6])
    nonhermitian = make(nh, frame)
    for index, time in enumerate((.31, -.29)):
        approx, estimate = nonhermitian.propagate(time)
        checks[f'nonhermitian_signed_bound{index}'] = np.linalg.norm(approx-expm(-1j*time*nh) @ frame) <= estimate['absolute_error_envelope_estimate']+2e-10
        checks[f'nonhermitian_growth{index}'] = max(r.get('projected_two_sided_growth', 0.) for r in estimate['columns']) > 1e-3
    for _ in range(5):
        nonhermitian.expand()
    approx, estimate = nonhermitian.propagate(.31)
    checks['nonhermitian_full_dense'] = np.linalg.norm(approx-expm(-.31j*nh) @ frame) < 2e-10
    q, t, residual = nonhermitian._decomposition(0)
    hermitized = np.exp(-.31j*nonhermitian.shifts[0])*q @ (expm(-.31j*(t+t.conj().T)/2)[:, 0]*nonhermitian.beta[0])
    mutations['silent_hermitization'] = float(np.linalg.norm(hermitized-approx[:, 0]))
    checks['mutation_silent_hermitization'] = mutations['silent_hermitization'] > 1e-6
    # Deliberately nonorthogonal Q and non-Galerkin T: independently form the
    # whole decomposition defect and nonzero initial residual, not Ritz surrogates.
    q = frame[:, :2]/4
    t = np.array([[.2+.1j, .3j], [.1, -.4-.2j]])
    e = nh @ q-q @ t
    x = q[:, 0]+.01*frame[:, 2]
    estimate = column_envelope(t, e, 1., np.linalg.norm(x-q[:, 0]), .19,
        ambient_growth=float(np.linalg.norm((nh-nh.conj().T)/(2j), 2)))
    defect_error = float(np.linalg.norm(expm(-.19j*nh) @ x-q @ expm(-.19j*t)[:, 0]))
    checks['whole_decomposition_and_initial_residual_enclosed'] = defect_error <= estimate['absolute_error_envelope_estimate']+2e-10
    checks['initial_reconstruction_not_omitted'] = estimate['initial_error_term'] > 1e-4
    defect_case = dict(actual_error=defect_error, **estimate)
    invariant = make(np.diag(np.arange(6.)), np.eye(6, 2, dtype=complex))
    stop = invariant.expand()
    state, estimate = invariant.propagate(.7)
    checks['exact_breakdown_dense'] = np.linalg.norm(state-expm(-.7j*np.diag(np.arange(6.)))[:, :2]) < 2e-10
    checks['exact_breakdown_stops'] = stop['status'] == 'numerical_breakdown'
    premature = make(h, frame, breakdown_rtol=.999999, shifts=np.zeros(frame.shape[1]))
    premature_stop = premature.expand()
    checks['deliberate_breakdown_exercised'] = premature_stop['accepted_columns'] == 0
    _, estimate = premature.propagate(.37)
    checks['breakdown_not_false_convergence'] = estimate['absolute_error_envelope_estimate'] > 1e-4
    limited = make(h, frame, cap=0)
    checks['work_cap_before_call'] = limited.expand()['status'] == 'action_limit' and limited.calls == []
    partial = make(h, frame, cap=2)
    partial.expand()
    checks['partial_batch_preserves_all_columns'] = partial.requested_columns == 2 and partial.propagate(.13)[0].shape == frame.shape
    try:
        ColumnFrameKrylov(lambda x: h @ x, frame, initial_action=h @ frame,
            operator_identity='changed', initial_action_identity='old', max_dimension=2,
            action_cap=0, max_memory_bytes=64*1024**2, ambient_growth_bound=0., operator_norm_bound=10.)
    except ValueError:
        checks['changed_operator_seed_rejected'] = True
    else:
        checks['changed_operator_seed_rejected'] = False
    try:
        ColumnFrameKrylov(lambda x: h @ x, frame, operator_identity='dense', max_dimension=2,
            action_cap=100, max_memory_bytes=1, ambient_growth_bound=0., operator_norm_bound=10.)
    except MemoryError:
        checks['memory_cap_before_allocation'] = True
    else:
        checks['memory_cap_before_allocation'] = False
    np.savez(args.output/'oracle-arrays.npz', **arrays)
    result = dict(schema='P254-column-frame-krylov-oracle-v1',
        status='passed' if all(checks.values()) else 'failed',
        checks={key: bool(value) for key, value in checks.items()}, symbolic=symbolic,
        stages=records, mutations=mutations, finite_defect_case=defect_case,
        invariants=dict(initial=initial_inv, final=final_inv), source_bindings=source,
        arrays=binding(args.output/'oracle-arrays.npz'), seconds=perf_counter()-started,
        limitations=['Independent finite dense oracle, not material execution',
                     'Ordinary floating bounds, not interval certificates or H-action rounding enclosures',
                     'No joint physical dynamics, long-time accuracy or speedup claim'])
    if sources() != source:
        raise ValueError('sources changed during oracle')
    publish(args.output/'result.json', result)
    if result['status'] != 'passed':
        raise AssertionError('column oracle failed; retained result identifies failing checks')


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
