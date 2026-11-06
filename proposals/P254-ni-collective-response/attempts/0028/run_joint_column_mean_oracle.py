#!/usr/bin/env python3
"""Independent checks for the new bounded controller and physical field bound.

The unchanged Arnoldi/envelope mathematics retains its earlier complete oracle.
"""
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.linalg import expm
import sympy as sp
from joint_column_mean import ColumnLimits, bounded_columns, field_norm_bound
from substrate_framework.material_physical_operators import PhysicalFieldAction
from substrate_framework.qe_material_hamiltonian import SpeciesProjectors
from substrate_framework.streamed_material_state import authenticate, file_binding

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sources():
    paths = [Path(__file__), HERE/'joint_column_mean.py', HERE/'run_column_krylov_material.py',
             HERE/'run_block_krylov_material.py']
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'column_frame_krylov.py', 'column_krylov_bounds.py', 'material_physical_operators.py',
        'material_frozen_action.py', 'qe_material_hamiltonian.py', 'partitioned_joint_mean.py',
        'partitioned_joint_slots.py', 'streamed_joint_mean.py', 'streamed_joint_reduction.py')]
    return [file_binding(path) for path in paths]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--registration', type=Path, required=True)
    p.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    args = p.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    registration = json.loads(args.registration.read_text())
    current = sources()
    expected = registration['source_bindings']
    if {b['path']:b['sha256'] for b in current} != {b['path']:b['sha256'] for b in expected}:
        raise ValueError('registered changed sources differ')
    for binding in expected:
        authenticate(binding)
    checks, values = {}, {}
    beta = np.array([[1.], [0.]], complex)
    species = SpeciesProjectors('independent-rank-one', (1,), beta, np.ones((1, 1)))
    field = PhysicalFieldAction(np.zeros((2, 3)), np.ones(3), ((species, beta.copy()),), 2)
    dense = field.apply(np.eye(2, dtype=complex))
    symbolic = sp.Matrix([1, 0])*sp.Matrix([[1, 0]])*2
    checks['independent_symbolic_two_leg_field'] = np.array_equal(dense, np.array(symbolic, complex))
    bound = field_norm_bound(field)
    checks['both_physical_legs_bound'] = np.linalg.norm(dense, 2) <= bound
    checks['removed_leg_mutation_detected'] = np.linalg.norm(dense, 2) > bound/2
    rng = np.random.default_rng(22528023)
    a = rng.normal(size=(7, 7))+1j*rng.normal(size=(7, 7))
    h = (a+a.conj().T)/2
    frame = rng.normal(size=(7, 3))+1j*rng.normal(size=(7, 3))
    frame[:, 2] = (.2+.3j)*frame[:, 0]+.7*frame[:, 1]
    for time in (.013, -.029):
        actual, receipt = bounded_columns(lambda x: h @ x, frame, time,
            limits=ColumnLimits(7, 8*1024**2, 1e-10), norm_bound=float(np.linalg.norm(h, 2)),
            growth_bound=0., identity='independent-dense-H-and-raw-frame')
        expected_frame = expm(-1j*time*h) @ frame
        error = float(np.linalg.norm(actual-expected_frame))
        checks['controlled_original_frame_'+str(time)] = error <= 1e-10
        checks['original_gram_'+str(time)] = np.linalg.norm(actual.conj().T @ actual-frame.conj().T @ frame) <= 1e-9
        values[str(time)] = dict(error=error, actual_dimensions=receipt['dimensions'],
                                estimated_error=receipt['envelope']['absolute_error_envelope_estimate'])
    try:
        bounded_columns(lambda x: h @ x, frame, 1., limits=ColumnLimits(1, 8*1024**2, 1e-14),
            norm_bound=float(np.linalg.norm(h, 2)), growth_bound=0., identity='insufficient-one-dimensional-action')
    except RuntimeError:
        checks['insufficient_allocation_is_not_a_result'] = True
    else:
        checks['insufficient_allocation_is_not_a_result'] = False
    checks = {name: bool(value) for name, value in checks.items()}
    result = dict(schema='P254-joint-column-controller-oracle-v1', status='passed' if all(checks.values()) else 'failed',
        source_bindings=current, registration=file_binding(args.registration), checks=checks, values=values,
        scope='Changed bounded controller and physical field bound only; prior full Arnoldi and joint-midpoint proofs retained')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
