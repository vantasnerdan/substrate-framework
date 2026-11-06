#!/usr/bin/env python3
"""Independent small dense algebra oracle; execute only with Root authorization."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace as NS
from time import perf_counter
import resource

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from substrate_framework.material_operator_covariance import (
    UnitaryPWMap, fourier_witness, covariance_bound, action_witness)
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_material_hamiltonian import PWFactors, SpeciesProjectors, MaterialHamiltonian
from threadpoolctl import threadpool_limits, threadpool_info


def binding(path):
    path = Path(path).resolve()
    digest = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return dict(path=str(path), bytes=str(path.stat().st_size), sha256=digest.hexdigest())


def authenticate(record):
    actual = binding(record['path'])
    if actual['sha256'] != record['sha256'] or int(actual['bytes']) != int(record['bytes']):
        raise ValueError('raw input bytes differ: '+record['path'])
    return actual


def explicit(path, digest):
    result = binding(path)
    if result['sha256'] != digest:
        raise ValueError('explicit digest differs: '+str(path))
    return result


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def sources():
    paths = [Path(__file__), HERE/'operator-covariance-registration.json',
             ROOT/'src/substrate_framework/material_operator_covariance.py',
             ROOT/'src/substrate_framework/qe_material_grid.py',
             ROOT/'src/substrate_framework/qe_material_hamiltonian.py']
    return [binding(path) for path in paths]


def dense_fft(fft):
    # Explicit exp(2 pi i m.x), independent of production FFT and bound algebra.
    x = np.array([(i/3, j/3, k/3) for k in range(3) for j in range(3) for i in range(3)])
    return np.exp(2j*np.pi*x@fft.miller)/np.sqrt(len(x))


def dense_operator(factors, potential):
    e = dense_fft(factors.fft)
    matrix = np.diag(factors.kinetic_ry.astype(complex)) + e.conj().T@(potential[:, None]*e)
    for species in factors.species:
        nh = len(species.active_d)
        for j in range(len(species.atom_indices)):
            b = species.beta[:, j*nh:(j+1)*nh]
            matrix += b@species.active_d@b.conj().T
    return matrix


def require(checks, name, condition):
    checks[name] = bool(condition)
    if not condition:
        raise AssertionError(name)


def run():
    rng = np.random.default_rng(2540028)
    shape = (3, 3, 3)
    m = np.array([(i, j, k) for k in range(3) for j in range(3) for i in range(3)]).T
    s = np.diag([-1, 1, 1]).astype(int)
    f = np.array([1/3, 0., 0.])
    indices = m.copy()
    pull = np.array([int((-i-1) % 3 + 3*j+9*k) for i, j, k in indices.T])
    rows = rng.permutation(27)
    phase = np.exp(-2j*np.pi*(m[:, rows].T@f))
    mapping = UnitaryPWMap(rows, phase, 'original-dense-source', 'nontrivial-unitary-target')
    u = np.zeros((27, 27), complex)
    for t in range(27):
        u[t, rows[t]] = phase[t]
    sf = PhysicalRhoFFT.create(shape, m, np.eye(3), 1.)
    tf = PhysicalRhoFFT.create(shape, s@m[:, rows], np.eye(3), 1.)
    b = (rng.normal(size=(27, 4))+1j*rng.normal(size=(27, 4)))/9
    d = np.diag([0.7, 1.4])
    angle = 0.37
    q = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    dt = q.T@d@q
    bt = np.concatenate([u@b[:, 2:4]@q, u@b[:, :2]@q], axis=1)
    kinetic = np.sum(m*m, axis=0).astype(float)
    v = rng.normal(size=27)
    checks, records = {}, []
    require(checks, 'independent_discrete_group_involution', np.array_equal(pull[pull], np.arange(27)))
    require(checks, 'independent_unitary_rows', np.linalg.norm(u.conj().T@u-np.eye(27)) < 2e-11)
    probes = rng.normal(size=(27, 5))+1j*rng.normal(size=(27, 5))
    for case in ('exact', 'kinetic', 'potential', 'geometry', 'conjugacy', 'fft-phase', 'intertwiner', 'combined'):
        target_t, target_v, target_b, target_d = kinetic[rows].copy(), v[pull].copy(), bt.copy(), dt.copy()
        current_map = mapping
        declared_q = q.copy()
        if case in ('intertwiner', 'combined'):
            declared_q[0, 1] += 0.05
        if case in ('kinetic', 'combined'):
            target_t[3] += 0.21
        if case in ('potential', 'combined'):
            target_v[7] += 0.31
        if case in ('geometry', 'combined'):
            # Actual phase response to a nonzero atom displacement, never repaired.
            target_b[:, :2] *= np.exp(-1j*(tf.miller.T@np.array([0.013, -0.007, 0.003])))[:, None]
        if case in ('conjugacy', 'combined'):
            target_d[0, 0] += 0.09
        if case in ('fft-phase', 'combined'):
            defect_phase = phase.copy(); defect_phase[4] *= np.exp(0.11j)
            current_map = UnitaryPWMap(rows, defect_phase, mapping.source_identity, mapping.target_identity)
            # Keep the nonlocal term exact under the changed row phase except
            # when the combined case explicitly includes geometry/D defects.
            if case == 'fft-phase':
                target_b = current_map.apply(np.concatenate([b[:, 2:4]@q, b[:, :2]@q], axis=1))
        source = PWFactors(sf, kinetic, (SpeciesProjectors('X', (1, 2), b, d),), np.zeros(3), None, ({'oracle': True},))
        target = PWFactors(tf, target_t, (SpeciesProjectors('X', (1, 2), target_b, target_d),), np.zeros(3), None, ({'oracle': True},))
        fw = fourier_witness(sf, tf, current_map, rotation=s, translation=f,
             reciprocal_shift=np.zeros(3, dtype=int), source_k_fractional=np.zeros(3), fft_pullback=pull)
        bound = covariance_bound(source, target, current_map, v, target_v,
              fft_pullback=pull, fourier=fw, atom_pullback=np.array([1, 0]), angular_intertwiners={'X': declared_q})
        hs, ht = dense_operator(source, v), dense_operator(target, target_v)
        actual_u = np.zeros_like(u)
        for t in range(27):
            actual_u[t, rows[t]] = current_map.phase[t]
        error = ht-actual_u@hs@actual_u.conj().T
        actual_norm = float(np.linalg.norm(error, 2))
        require(checks, case+'_dense_bound', actual_norm <= bound['total_bound_ry']+2e-11*(1+actual_norm))
        expected_fourier = dense_fft(tf)@actual_u-dense_fft(sf)[pull]
        require(checks, case+'_independent_fourier_defect', np.isclose(np.linalg.norm(expected_fourier), fw['embedding_defect_frobenius'], atol=2e-11, rtol=2e-11))
        source_h = MaterialHamiltonian(NS(fft=sf, inputs=NS(sources=()), vrs=np.column_stack([v, v])), source, 1)
        target_h = MaterialHamiltonian(NS(fft=tf, inputs=NS(sources=()), vrs=np.column_stack([target_v, target_v])), target, 1)
        require(checks, case+'_independent_production_action', np.allclose(source_h.apply(probes), hs@probes, atol=2e-11, rtol=2e-11) and np.allclose(target_h.apply(probes), ht@probes, atol=2e-11, rtol=2e-11))
        witness = action_witness(source_h.apply, target_h.apply, current_map, probes, block_size=2)
        direct = np.linalg.norm(ht@actual_u@probes-actual_u@hs@probes, axis=0)
        require(checks, case+'_independent_action_residual', np.allclose(direct, witness['residual_norms_ry'], atol=2e-11, rtol=2e-11))
        if case == 'exact':
            require(checks, 'exact_group_symmetry', actual_norm < 2e-11 and bound['total_bound_ry'] < 2e-11)
        elif case == 'intertwiner':
            require(checks, 'unfitted_intertwiner_defects_retained', actual_norm < 2e-11 and bound['total_bound_ry'] > 1e-6 and all(a['d_conjugacy_defect_norm2'] > 1e-6 and a['projector_residual_frobenius'] > 1e-6 for a in bound['atoms']))
        elif case != 'combined':
            term = {'kinetic': 'kinetic_bound_ry', 'potential': 'local_bound_ry',
                    'geometry': 'nonlocal_bound_ry', 'conjugacy': 'nonlocal_bound_ry',
                    'fft-phase': 'local_fourier_bound_ry'}[case]
            require(checks, case+'_omitted_term_mutation_detected', actual_norm > bound['total_bound_ry']-bound[term]+1e-6)
        else:
            require(checks, 'combined_all_defects_retained', bound['kinetic_bound_ry'] > 0 and bound['local_potential_defect_ry'] > 0 and bound['local_fourier_bound_ry'] > 0 and all(a['d_conjugacy_defect_norm2'] > 0 and a['projector_residual_frobenius'] > 0 for a in bound['atoms']))
        records.append(dict(case=case, dense_norm_ry=actual_norm, bound=bound, action=witness))
    mutations = {
        'duplicate_rows': lambda: UnitaryPWMap(np.zeros(27, dtype=int), phase, 's', 't'),
        'nonunit_phase': lambda: UnitaryPWMap(rows, phase*1.1, 's', 't'),
        'antiunitary': lambda: UnitaryPWMap(rows, phase, 's', 't', antiunitary=True),
        'wrong_fft_pullback': lambda: fourier_witness(sf, tf, mapping, rotation=s, translation=f, reciprocal_shift=np.zeros(3, dtype=int), source_k_fractional=np.zeros(3), fft_pullback=np.arange(27)),
        'duplicate_fft_pullback': lambda: fourier_witness(sf, tf, mapping, rotation=s, translation=f, reciprocal_shift=np.zeros(3, dtype=int), source_k_fractional=np.zeros(3), fft_pullback=np.zeros(27, dtype=int)),
    }
    for name, function in mutations.items():
        rejected = False
        try:
            function()
        except ValueError:
            rejected = True
        require(checks, name+'_rejected', rejected)
    return checks, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    before = sources()
    started = perf_counter()
    with threadpool_limits(limits=1):
        checks, records = run()
        pools = threadpool_info()
    if sources() != before:
        raise ValueError('source bytes changed during oracle execution')
    publish(args.output, dict(schema='P254-operator-covariance-oracle-v1', status='passed',
        checks=checks, cases=records, source_bindings=before, seconds=perf_counter()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        cpu_affinity=sorted(os.sched_getaffinity(0)), threadpools=pools,
        qualification='independent numerical dense algebra and mutations, not interval certification or actual-material symmetry'))


if __name__ == '__main__':
    main()
