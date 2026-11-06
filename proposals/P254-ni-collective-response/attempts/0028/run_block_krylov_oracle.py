#!/usr/bin/env python3
"""Independent dense/SymPy proof and sensitive mutations; no material run."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from scipy.linalg import expm
import sympy as sp
from threadpoolctl import threadpool_limits
from substrate_framework.full_frame_krylov import FullFrameKrylov, physical_invariants


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return dict(path=str(path), sha256=digest.hexdigest(), bytes=path.stat().st_size)


def authenticate(entry):
    if binding(entry['path']) != dict(entry, bytes=int(entry['bytes'])):
        raise ValueError('input binding differs: '+entry['path'])


def explicit(path, digest):
    result = binding(path)
    if result['sha256'] != digest:
        raise ValueError('explicit digest mismatch: '+str(path))
    return result


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def sources():
    paths = [HERE/name for name in ('block-krylov-registration.json',
             'run_block_krylov_oracle.py', 'run_block_krylov_material.py')]
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'full_frame_krylov.py', 'material_frozen_action.py', 'qe_material_hamiltonian.py',
        'streamed_material_evolution.py', 'streamed_sector_executor.py')]
    paths.append(HERE/'run_streamed_material_evolution.py')
    return [binding(p) for p in paths]


def symbolic_proof():
    t = sp.symbols('t', real=True)
    h = sp.Matrix([[0, sp.I], [-sp.I, 0]])
    frame = sp.Matrix([[2, 1], [0, 3]])
    d = sp.Matrix([[1, sp.I/3], [-sp.I/3, 2]])
    direct = (-sp.I*t*h).exp()
    formula = sp.cos(t)*sp.eye(2)-sp.I*sp.sin(t)*h
    simplify = lambda x: x.applyfunc(sp.simplify)
    state = direct*frame
    rho = state*d*state.conjugate().T
    initial_rho = frame*d*frame.conjugate().T
    q = sp.Matrix([1, 0])
    # Exact Duhamel identity for a saved-span approximation T=0, E=Hq.
    s = sp.symbols('s', real=True)
    integral = sp.integrate((-sp.I*(t-s)*h).exp()*(-sp.I*h*q), (s, 0, t))
    checks = dict(exponential=simplify(direct-formula) == sp.zeros(2),
        differential_equation=simplify(sp.diff(state, t)+sp.I*h*state) == sp.zeros(2),
        initial_state=simplify(state.subs(t, 0)-frame) == sp.zeros(2),
        duhamel_complement=simplify(direct*q-q-integral) == sp.zeros(2, 1),
        raw_D_noncommuting=simplify(h*initial_rho-initial_rho*h) != sp.zeros(2),
        number_conserved=sp.simplify(sp.trace(rho)-sp.trace(initial_rho)) == 0,
        purity_conserved=sp.simplify(sp.trace(rho*rho)-sp.trace(initial_rho*initial_rho)) == 0)
    return checks, dict(H=str(h), frame=str(frame), raw_D=str(d),
                        exponential=str(simplify(direct)), duhamel=str(simplify(integral)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    source = sources()
    start = perf_counter()
    checks, symbolic = symbolic_proof()
    rng = np.random.default_rng(22528071)
    a = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
    h = (a+a.conj().T)/3
    frame = np.zeros((6, 2), dtype=complex)
    frame[:2] = [[1.7, .4+.3j], [.1j, .8]]
    d = np.array([[.9, .17+.08j], [.17-.08j, .4]])
    times = [.0, .13, -.21, .37]
    calls = []
    def action(x):
        calls.append(x.shape[1])
        return h @ x
    def make(matrix, initial=frame, cap=4, blocks=3, basis=6, counted=None):
        return FullFrameKrylov(counted or (lambda x: matrix @ x), initial,
            initial_action=matrix @ initial, max_blocks=blocks, max_basis=basis,
            action_cap=cap, max_memory_bytes=64*1024**2,
            ambient_growth_bound=float(np.linalg.norm((matrix-matrix.conj().T)/(2j), 2)))
    model = make(h, counted=action)
    checks['seed_action_reused'] = calls == []
    records, arrays = [], dict(H=h, frame=frame, raw_D=d)
    errors = []
    for stage_index in range(3):
        stage = model.stage()
        rows = []
        for index, time in enumerate(times):
            actual, estimate = model.propagate(time)
            exact = expm(-1j*time*h) @ frame
            error = float(np.linalg.norm(actual-exact))
            checks[f'bound_stage{stage_index}_time{index}'] = error <= estimate['absolute_error_envelope_estimate']+2e-10
            arrays[f'stage{stage_index}_time{index}'] = actual
            arrays[f'exact_time{index}'] = exact
            rows.append(dict(error=error, **estimate))
        errors.append(rows[-1]['error'])
        records.append(dict(diagnostics=stage['diagnostics'], times=rows))
        arrays[f'residual{stage_index}'] = stage['residual']
        if stage_index < 2:
            model.expand()
    checks['full_ambient_dense_match'] = max(row['error'] for row in records[-1]['times']) < 2e-10
    checks['ambient_complement_not_closed_initially'] = errors[0] > 1e-3
    checks['refinement_sensitive'] = errors[2] < errors[1] < errors[0]
    checks['all_original_columns_retained'] = model.stage()['coefficients'].shape == (6, 2)
    checks['exact_action_column_accounting'] = calls == [2, 2] and model.requested_columns == 4
    final, _ = model.propagate(times[-1])
    initial_inv, final_inv = physical_invariants(frame, d), physical_invariants(final, d)
    checks['raw_D_number'] = abs(initial_inv['number']['real']-final_inv['number']['real']) < 2e-10
    checks['raw_D_purity'] = abs(initial_inv['purity']['real']-final_inv['purity']['real']) < 2e-10
    rho = final @ d @ final.conj().T
    q = np.linalg.qr(frame, mode='reduced')[0]
    mutated = dict(
        omission=final[:, :1] @ d[:1, :1] @ final[:, :1].conj().T,
        sign=(expm(1j*times[-1]*h) @ frame) @ d @ (expm(1j*times[-1]*h) @ frame).conj().T,
        complement=arrays['stage0_time3'] @ d @ arrays['stage0_time3'].conj().T,
        normalization=(expm(-1j*times[-1]*h) @ q) @ d @ (expm(-1j*times[-1]*h) @ q).conj().T,
        raw_D_diagonal=final @ np.diag(np.diag(d)) @ final.conj().T)
    mutations = {name: float(np.linalg.norm(value-rho)) for name, value in mutated.items()}
    for name, error in mutations.items():
        checks['mutation_'+name] = error > 1e-6
    invariant_h = np.diag(np.arange(6.)+.4).astype(complex)
    invariant_h[:2, :2] = [[.4, .2j], [-.2j, 1.4]]
    invariant = make(invariant_h)
    state, estimate = invariant.propagate(.7)
    checks['invariant_subspace_exact'] = np.linalg.norm(state-expm(-.7j*invariant_h) @ frame) < 2e-10
    checks['invariant_defect_small'] = invariant.stage()['diagnostics']['residual_frobenius'] < 2e-10
    closed, discarded = invariant.expand()
    checks['invariant_numerical_closure'] = closed['accepted_columns'] == 0
    nh = h+1j*np.diag([.5, -.4, .2, .1, -.3, .6])
    nonhermitian = make(nh)
    for index, time in enumerate((.31, -.29)):
        approx, estimate = nonhermitian.propagate(time)
        exact = expm(-1j*time*nh) @ frame
        checks[f'nonhermitian_bound{index}'] = np.linalg.norm(approx-exact) <= estimate['absolute_error_envelope_estimate']+2e-10
        checks[f'projected_growth_retained{index}'] = estimate['projected_two_sided_growth'] > .1
    nonhermitian.expand()
    nonhermitian.expand()
    approx, estimate = nonhermitian.propagate(.31)
    checks['nonhermitian_full_dense'] = np.linalg.norm(approx-expm(-.31j*nh) @ frame) < 2e-10
    st = nonhermitian.stage()
    hermitized = nonhermitian.q @ expm(-.31j*(st['projected']+st['projected'].conj().T)/2) @ st['coefficients']
    mutations['hermitization'] = float(np.linalg.norm(hermitized-approx))
    checks['mutation_hermitization'] = mutations['hermitization'] > 1e-6
    blocked_calls = []
    limited = make(h, cap=0, counted=lambda x: blocked_calls.append(x.shape[1]) or h @ x)
    stopped, _ = limited.expand()
    checks['cap_stops_before_action'] = stopped['status'] == 'action_limit' and blocked_calls == []
    checks['insufficient_target_not_relaxed'] = records[0]['times'][-1]['absolute_error_envelope_estimate'] > 1e-12
    np.savez(args.output/'oracle-arrays.npz', **arrays)
    result = dict(schema='P254-full-frame-block-krylov-oracle-v1',
                  status='passed' if all(checks.values()) else 'failed',
                  checks={key: bool(value) for key, value in checks.items()},
                  symbolic=symbolic, stages=records, mutations=mutations,
                  invariants=dict(initial=initial_inv, final=final_inv),
                  arrays=binding(args.output/'oracle-arrays.npz'), source_bindings=source,
                  seconds=perf_counter()-start,
                  limitations=['Finite independent proof only; no material execution',
                    'Floating-point analytic bound estimates, not interval certificates',
                    'No nonlinear/changing-H, all-sector, long-time or speedup claim'])
    if sources() != source:
        raise ValueError('source changed during proof')
    publish(args.output/'result.json', result)
    if result['status'] != 'passed':
        raise AssertionError('sensitive block Krylov proof failed; see retained receipts')


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
