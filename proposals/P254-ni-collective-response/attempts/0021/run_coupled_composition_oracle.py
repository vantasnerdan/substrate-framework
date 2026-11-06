#!/usr/bin/env python3
"""Narrow new composition proof against original-basis coupled EOM; no Native run."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import sympy as sp

import run_electronic_time_oracle as original
from substrate_framework.pw_coupled_dynamics import PWCoupledDynamics
from substrate_framework.pw_electronic_dynamics import PWElectronicDynamics


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return {"path": str(path), "bytes": str(len(data)), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('exclusive output required')
    started = time.monotonic()
    checks = []
    def compare(name, value, reference):
        error = float(np.linalg.norm(value-reference)/max(1., np.linalg.norm(reference)))
        checks.append({'name': name, 'scaled_error': error, 'passed': error <= 1e-10})
    fixture = original.load_fixture()
    for step in (0, 1):
        f = fixture.fixture(step)
        electronic = PWElectronicDynamics(f.model)
        spin_kernel = np.array([[.3, .08], [.08, .24]])
        kernel = lambda n: n @ spin_kernel.T
        frozen = .7*np.eye(6)+.02*np.ones((6, 6))
        masses = np.array([7., 7., 7., 2., 2., 2.])
        model = PWCoupledDynamics(electronic, kernel, frozen, masses)
        rng = np.random.default_rng(25421+step)
        initial = .02*fixture.complex_random(rng, (model.state_size,))
        delta = original.original_state(f, electronic, initial[12:], fixture)
        full_initial = np.r_[initial[:12], delta.ravel()]
        shape = delta.shape
        def direct_rhs(tau, state):
            u, p, rho = state[:6], state[6:12], state[12:].reshape(shape)
            obs = original.original_observables(f, rho, fixture)
            hs, _, forcing = original.original_operators(f, kernel(obs.density), u, fixture)
            drho = original.original_derivative(hs, rho, forcing)
            return np.r_[p/masses, -frozen @ u-obs.ionic_gradient, drho.ravel()]
        def mapped(state):
            return np.r_[state[:12], original.original_state(f, electronic, state[12:], fixture).ravel()]
        compare(f'q{step}/combined_original_rhs', mapped(model.derivative(initial)), direct_rhs(0, full_initial))
        generator = np.column_stack([direct_rhs(0, column) for column in np.eye(full_initial.size, dtype=complex)])
        times = np.array([0., .2, .7])
        solution = solve_ivp(model.derivative_real, (0., .7), initial.view(float), t_eval=times,
                             method='DOP853', rtol=2e-13, atol=2e-14)
        if not solution.success or solution.y.shape[1] != len(times):
            raise RuntimeError(solution.message)
        for index, tau in enumerate(times):
            actual = np.ascontiguousarray(solution.y[:, index]).view(complex)
            compare(f'q{step}/real_IVP_original_exponential/{index}', mapped(actual), expm(tau*generator) @ full_initial)
        gradient = electronic.observables(initial[12:]).ionic_gradient
        wrong = model.derivative(initial)
        wrong[6:12] += 2*gradient
        separation = float(np.linalg.norm(mapped(wrong)-direct_rhs(0, full_initial)))
        checks.append({'name': f'q{step}/wrong_force_sign_exposed', 'separation': separation, 'passed': separation > 1e-8})
        obs, mechanical = model.mechanical_observables(initial)
        velocity, acceleration = initial[6:12]/masses, model.derivative(initial)[6:12]/masses
        direct_power = np.vdot(initial[6:12], acceleration).real + np.vdot(velocity, frozen @ initial[:6]).real
        compare(f'q{step}/ionic_power_identity', np.asarray(direct_power), np.asarray(mechanical['electronic_power_into_ions_ry_per_tau']))
    transform = sp.Matrix([[1, 1], [1, -1]])
    checks.append({'name': 'symbolic/two_cell_inverse_no_atomic_phase', 'passed': transform.T*transform/2 == sp.eye(2)})
    kick = np.zeros((2, 1, 1, 6), complex)
    kick[0, 0, 0, 4] = 1e-4
    qkick = np.fft.fftn(kick, axes=(0, 1, 2))
    compare('Fourier/local_kick_direct_sum', qkick[:, 0, 0], np.asarray(transform, complex) @ kick[:, 0, 0])
    compare('Fourier/local_kick_inverse', np.fft.ifftn(qkick, axes=(0, 1, 2)), kick)
    wrong = qkick.copy()
    wrong[1] *= np.exp(1j*np.pi/3)
    checks.append({'name': 'Fourier/extra_atomic_phase_exposed', 'passed': bool(np.linalg.norm(np.fft.ifftn(wrong, axes=(0, 1, 2))-kick) > 1e-8)})
    import substrate_framework.pw_coupled_dynamics as module
    report = {'preregistration': '2cb8638a', 'scope': 'new composed full-PW/inertial EOM and interleaved real IVP; finite mathematical fixture, not actual material time comparison',
              'checks': checks, 'passed': all(x['passed'] for x in checks), 'wall_seconds': time.monotonic()-started,
              'source_bindings': [binding(x) for x in (Path(__file__), Path(module.__file__), Path(original.__file__), original.FIXTURE_PATH)]}
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'output': str(args.output), 'passed': report['passed'], 'wall_seconds': report['wall_seconds']}))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
