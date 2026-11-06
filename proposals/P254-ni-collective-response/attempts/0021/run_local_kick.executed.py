#!/usr/bin/env python3
"""Actual full-PW ionic/electronic local-kick pilot; physical-time checkpointed."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback

import numpy as np
from scipy.integrate import solve_ivp
from substrate_framework.pw_coupled_dynamics import PWCoupledDynamics
from substrate_framework.pw_electronic_dynamics import PWElectronicDynamics
from substrate_framework.qe_ionic_mass import ionic_mass_diagonal
from substrate_framework.qe_pw_response_frames import build_pw_response_models
from substrate_framework.spin_pbe_kernel import apply_spin_pbe_kernel
from time_pilot_checkpoints import (HERE, binding, save_json, implementation, save_sample,
                                    recover_samples, cell_trajectory)
from pilot_inputs import load_inputs


def now():
    return datetime.now(timezone.utc).isoformat()


def interrupted(signum, frame):
    raise InterruptedError(f'signal {signum}; completed full physical-time states retained')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    target = args.output_dir.resolve()
    if target.exists() and not args.resume:
        parser.error('existing output requires explicit same-generation --resume')
    target.mkdir(parents=True, exist_ok=args.resume)
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, interrupted)
    settings = json.loads((HERE/'registration.json').read_bytes())['actual_pilot']
    times = np.linspace(*settings['time_window_tau_ry'], settings['output_samples'])
    scales = {'rtol': 1e-4, 'displacement_atol_bohr': 1e-12, 'velocity_atol_bohr_per_tau': 1e-12,
              'sqrt_weighted_occupation_amplitude_atol': 1e-10,
              'norm': 'SciPy standard real RMS; weighted absolute electronic scale, not per-component error guarantee'}
    identity = {'implementation': implementation(), 'registration': binding(HERE/'registration.json'),
                'input_inventory': binding(HERE.parent/'0020/actual-input-inventory.json'), 'scales': scales}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    started = time.time()
    allocation = {'started_unix': started, 'deadline_unix': started+settings['total_wall_cap_seconds'], 'sectors': {}}
    records, admitted = {}, False
    report = {'started_at': now(), 'status': 'running', 'settings': settings, 'scales': scales,
              'preregistration_commit': '2cb8638a', 'command': [sys.executable, *sys.argv],
              'environment': {'cpu_affinity': [str(x) for x in sorted(os.sched_getaffinity(0))],
                              'BLAS_threads': os.environ.get('OPENBLAS_NUM_THREADS')},
              'recovery': 'Full spectral state at completed quarter-tau outputs; adaptive integrator history is not saved. Work since last committed output can be lost. Original wall/action allocation never resets.'}
    def event(kind, **values):
        record = {'at': now(), 'event': kind, 'elapsed_seconds': time.time()-allocation['started_unix'], **values}
        with (target/'telemetry.jsonl').open('a') as stream:
            stream.write(json.dumps(record, allow_nan=False)+'\n')
            stream.flush()
        print(json.dumps(record, allow_nan=False), flush=True)
    def checkpoint():
        report['completed_samples'] = list(records.values())
        report['elapsed_seconds'] = time.time()-allocation['started_unix']
        save_json(target/'result.json', report)
    def budget(point=None):
        current = time.time()
        if current >= allocation['deadline_unix']:
            raise TimeoutError('600s fixed cumulative allocation exhausted')
        if not any(x['tau_ry'] > 0 for x in records.values()) and current-allocation['started_unix'] > settings['first_nonzero_time_state_due_seconds']:
            raise TimeoutError('no completed nonzero physical-time output by150s; method decision required')
        if point is not None and current-allocation['sectors'][point]['started_unix'] > settings['per_sector_wall_cap_seconds']:
            raise TimeoutError('240s fixed sector allocation exhausted: '+point)
    try:
        if args.resume:
            previous = json.loads((target/'identity.json').read_bytes())
            if previous['fingerprint'] != fingerprint:
                raise ValueError('source/input/scales differ; original result remains untouched')
            allocation = json.loads((target/'allocation.json').read_bytes())
            records = recover_samples(target, fingerprint, times)
            report['reference'] = json.loads((target/'result.json').read_bytes()).get('reference')
        else:
            save_json(target/'identity.json', {'fingerprint': fingerprint, **identity})
            save_json(target/'allocation.json', allocation)
        admitted = True
        event('pilot_started', output_dir=str(target), recovered_samples=str(len(records)))
        bundles, frozen, input_receipts = load_inputs()
        report['actual_inputs'] = input_receipts
        masses = ionic_mass_diagonal(bundles['gamma'])
        for point, expected in (('gamma', [0., 0., 0.]), ('half', [.5, 0., 0.])):
            bg = bundles[point].background
            qfrac = bg.array('at').T @ bg.array('q_cart_2pi_alat')
            if not np.allclose(qfrac, expected, atol=2e-12, rtol=0):
                raise ValueError('actual q does not occupy declared complete two-cell grid')
            if not np.array_equal(ionic_mass_diagonal(bundles[point]), masses):
                raise ValueError('source sectors have different actual ionic inertia')
        report['masses_ry'] = masses.tolist()
        report['mechanical_scope'] = 'Per-q ionic kinetic+raw frozen quadratic and instantaneous electronic power only; no actual total electronic/orbit energy or equilibrium covariance is claimed. q-to-cell transform has normalization1/2, no extra atom phase. Electronic density arrays are q-periodic envelopes, not fabricated real-space cell densities.'
        kick = np.zeros((2, 1, 1, 6), complex)
        kick[0, 0, 0, 4] = 1e-4
        qkick = np.fft.fftn(kick, axes=(0, 1, 2))[:, 0, 0]
        models = None
        if len(records) < 2*len(times):
            budget()
            frame_start = time.monotonic()
            models = build_pw_response_models(bundles['gamma'], bundles['half'], bundles['minus'])
            report['reference'] = models['gamma'].metadata
            event('frames_constructed', seconds=time.monotonic()-frame_start)
        checkpoint()
        for qindex, point in enumerate(('gamma', 'half')):
            present = [index for index in range(len(times)) if f'{point}-t{index:02d}' in records]
            if len(present) == len(times):
                event('completed_sector_reused', point=point)
                continue
            if point not in allocation['sectors']:
                allocation['sectors'][point] = {'started_unix': time.time()}
                save_json(target/'allocation.json', allocation)
            budget(point)
            electronic = PWElectronicDynamics(models[point])
            coupled = PWCoupledDynamics(electronic, lambda n: apply_spin_pbe_kernel(bundles[point], n),
                                        frozen[point].cartesian_matrix, masses)
            actions = 0
            action_path = target/(point+'-rhs-actions.jsonl')
            if action_path.exists():
                actions = max((int(json.loads(line)['action']) for line in action_path.read_text().splitlines()), default=0)
            if present:
                index = present[-1]
                with np.load(records[f'{point}-t{index:02d}']['archive']['path'], allow_pickle=False) as stored:
                    state = stored['state'].copy()
                if state.shape != (coupled.state_size,):
                    raise ValueError('saved full state differs from actual current model')
            else:
                index, state = 0, coupled.zero_state()
                state[6:12] = masses*qkick[qindex]
                observable, mechanical = coupled.mechanical_observables(state)
                records[f'{point}-t00'] = save_sample(target, point, 0, 0., state, observable, mechanical, fingerprint, now())
                checkpoint()
            event('sector_started', point=point, state_size=str(coupled.state_size), resumed_tau=float(times[index]))
            def rhs(tau, y):
                nonlocal actions
                budget(point)
                if actions >= settings['max_rhs_actions_per_sector']:
                    raise TimeoutError('400 cumulative RHS action cap reached: '+point)
                actions += 1
                with action_path.open('a') as stream:
                    stream.write(json.dumps({'action': str(actions), 'tau_ry': float(tau), 'at': now()})+'\n')
                    stream.flush()
                return coupled.derivative_real(tau, y)
            for next_index in range(index+1, len(times)):
                solution = solve_ivp(rhs, (times[next_index-1], times[next_index]), state.view(float),
                    t_eval=[times[next_index]], method='DOP853', rtol=scales['rtol'],
                    atol=coupled.absolute_tolerances())
                if not solution.success or solution.y.shape[1] != 1:
                    raise RuntimeError(solution.message)
                state = np.ascontiguousarray(solution.y[:, -1]).view(complex)
                observable, mechanical = coupled.mechanical_observables(state)
                name = f'{point}-t{next_index:02d}'
                first_nonzero = not any(x['tau_ry'] > 0 for x in records.values())
                records[name] = save_sample(target, point, next_index, times[next_index], state, observable, mechanical, fingerprint, now())
                checkpoint()
                event('physical_time_completed', point=point, tau_ry=float(times[next_index]), rhs_actions=str(actions),
                      displacement_norm=float(np.linalg.norm(state[:6])), density_norm=float(np.linalg.norm(observable.density)))
                if first_nonzero and time.time()-allocation['started_unix'] > settings['first_nonzero_time_state_due_seconds']:
                    raise TimeoutError('first nonzero-time output committed after150s; retained, further work stopped')
                budget(point)
            event('sector_completed', point=point, rhs_actions=str(actions))
        report['trajectory_archive'] = cell_trajectory(target, records, masses, times)
        report['status'] = 'complete'
    except Exception as error:
        report['status'] = 'stopped_requires_root_decision'
        report['failure'] = {'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()}
    finally:
        report['finished_at'] = now()
        if admitted:
            checkpoint()
        else:
            save_json(target/'resume-refusal.json', report)
        if report['status'] == 'complete':
            (target/'.gate_done_local_kick').write_text(report['finished_at']+'\n')
        event('pilot_finished', status=report['status'], samples=str(len(records)))
        print('[dft-done]' if report['status'] == 'complete' else '[dft-fail]', 'p254-local-kick-0021', flush=True)
    return 0 if report['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
