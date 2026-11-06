#!/usr/bin/env python3
"""Concurrency/checkpoint delta oracle; reuse the bound original scientific result."""
import argparse
import json
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from substrate_framework.streamed_material_state import DiskMean


def state_check(directory, frames, coefficients):
    records = [(dict(identity=str(i), component=i % 2, weight=.5, basis_identity='basis-'+str(i),
                     gauge_identity='raw-'+str(i), source_identity='source-'+str(i)), psi, d)
               for i, (psi, d) in enumerate(zip(frames, coefficients))]
    state = DiskMean.create(directory, records, identity='dense-two-spin-full-coverage',
                            time=0., u=[.07], p=[-.03])
    restored = DiskMean(state.path, expected_identity=state.identity)
    for s, psi, d in zip(restored.slots, frames, coefficients):
        a, b = restored.arrays(s)
        assert np.array_equal(a, psi) and np.array_equal(b, d)
        assert isinstance(a, np.memmap) and isinstance(b, np.memmap)
    return restored


def concurrency_checks(args, reference, factory, kwargs, frames, ds):
    """Delta only: reuse the authenticated scientific oracle, do not rerun it."""
    from contextlib import contextmanager
    from dataclasses import replace
    from hashlib import sha256
    import multiprocessing as mp
    import os
    import time
    from threadpoolctl import threadpool_info
    from substrate_framework.streamed_material_evolution import StreamedMeanEvolution
    from substrate_framework.streamed_material_state import file_binding, authenticate
    from substrate_framework.streamed_sector_executor import (
        SectorExecutor, adopt_initial_boundary, physical_identity)
    prior_binding = file_binding(args.prior_oracle)
    if prior_binding['sha256'] != args.prior_oracle_sha256:
        raise ValueError('prior scientific oracle digest differs')
    prior = json.loads(args.prior_oracle.read_text())
    if prior['status'] != 'passed' or prior['schema'] != 'P254-streamed-evolution-oracle-v1':
        raise ValueError('passed prior independent scientific oracle required')
    preserved = []
    for entry in prior['source_bindings']:
        try:
            authenticate(entry)
            preserved.append(entry)
        except ValueError:
            original = Path(entry['path'])
            archive = HERE/'source-generations'/(original.stem+'-'+entry['sha256']+original.suffix)
            archived = file_binding(archive)
            assert archived['sha256'] == entry['sha256'] and archived['bytes'] == int(entry['bytes'])
            preserved.append(archived)
    executor = SectorExecutor(args.worker_cpus)
    if len(executor.cpus) < 2:
        raise ValueError('parallel delta oracle requires at least two explicitly allocated CPUs')
    def probe(_):
        return dict(affinity=sorted(os.sched_getaffinity(0)), pools=threadpool_info())
    probes = list(executor.map(probe, range(len(executor.cpus))))
    assert all(len(row['affinity']) == 1 and row['affinity'][0] in executor.cpus
               and all(pool['num_threads'] == 1 for pool in row['pools']) for row in probes)
    active = mp.get_context('fork').Value('i', 0)
    completed = mp.get_context('fork').Value('i', 0)
    actual_h0 = mp.get_context('fork').Value('i', 0)
    completion_order = mp.get_context('fork').Array('i', 4)
    @contextmanager
    def delayed_factory(slot):
        with active.get_lock():
            active.value += 1
        try:
            # Different finish order, same deterministic physical reduction.
            time.sleep(.01*(3-int(slot.identity)))
            with factory(slot) as (coupling, h0):
                def counted(v):
                    with actual_h0.get_lock():
                        actual_h0.value += 1
                    return h0(v)
                yield coupling, counted
        finally:
            with completed.get_lock():
                completion_order[completed.value % 4] = int(slot.identity)
                completed.value += 1
            with active.get_lock():
                active.value -= 1
    old_kernel = reference.kernel
    barrier_observations = []
    def checked_kernel(n):
        assert active.value == 0 and completed.value % 4 == 0
        barrier_observations.append(completed.value)
        return old_kernel(n)
    parallel_reference = replace(reference, kernel=checked_kernel)
    with tempfile.TemporaryDirectory(prefix='P254-sector-delta-') as temporary:
        base = Path(temporary)
        state = state_check(base/'initial', frames, ds)
        serial = StreamedMeanEvolution(reference, factory, **kwargs)
        parallel = StreamedMeanEvolution(parallel_reference, delayed_factory, executor=executor, **kwargs)
        serial_stage = serial.reduce(state, energy=True)
        parallel_stage = parallel.reduce(state, energy=True)
        assert np.array_equal(serial_stage.density, parallel_stage.density)
        assert np.array_equal(serial_stage.gradient, parallel_stage.gradient)
        assert serial_stage.energy == parallel_stage.energy
        observed_completion_order = list(completion_order[:])
        assert observed_completion_order != list(range(4))
        serial_final = serial.step(state, base/'serial', dt=.075, max_iterations=40, tolerance=2e-13)
        parallel_final = parallel.step(state, base/'parallel', dt=.075, max_iterations=40, tolerance=2e-13)
        equation_error = max(float(np.linalg.norm(serial_final.arrays(s)[0]-parallel_final.arrays(t)[0]))
                             for s, t in zip(serial_final.slots, parallel_final.slots))
        assert equation_error < 2e-12
        assert np.linalg.norm(serial_final.u-parallel_final.u) < 2e-12
        assert np.linalg.norm(serial_final.p-parallel_final.p) < 2e-12
        assert all(a.coefficient == b.coefficient for a, b in zip(state.slots, parallel_final.slots))
        payload = json.loads(state.path.read_text())
        payload['slots'].reverse()
        reverse_path = base/'reverse.json'
        reverse_path.write_text(json.dumps(payload))
        reverse = DiskMean(reverse_path)
        reversed_model = StreamedMeanEvolution(parallel_reference, delayed_factory, executor=executor,
            **dict(kwargs, expected_slots=tuple(reversed(kwargs['expected_slots']))))
        reverse_final = reversed_model.step(reverse, base/'reverse-step', dt=.075,
            max_iterations=40, tolerance=2e-13)
        reversed_frames = {s.identity: reverse_final.arrays(s)[0] for s in reverse_final.slots}
        permutation_error = max(float(np.linalg.norm(serial_final.arrays(s)[0]-reversed_frames[s.identity]))
                                for s in serial_final.slots)
        assert permutation_error < 2e-12
        assert np.linalg.norm(reverse_final.u-serial_final.u) < 2e-12
        before_calls = actual_h0.value
        capped = StreamedMeanEvolution(reference, delayed_factory, executor=executor,
                                      max_h0_actions=3, **kwargs)
        try:
            capped.step(state, base/'capped', dt=.075, max_iterations=40, tolerance=2e-13)
        except RuntimeError as error:
            assert 'action cap exhausted' in str(error)
        else:
            raise AssertionError('shared global action cap was bypassed')
        assert capped.counts['h0_calls'] == 3 and actual_h0.value-before_calls == 3
        assert not (base/'capped'/'boundary.json').exists()
        def failing_frame(slot):
            if slot.identity == '2':
                raise RuntimeError('injected worker frame failure')
            return state.arrays(slot)[0]
        try:
            executor.successor(state, base/'failed-publication', failing_frame,
                time=.1, u=state.u, p=state.p, evidence={}, counts=lambda: {})
        except RuntimeError as error:
            assert 'injected worker' in str(error)
        else:
            raise AssertionError('failed sector published a complete boundary')
        assert not (base/'failed-publication'/'boundary.json').exists()
        assert (base/'failed-publication'/'frame-0000.npy').exists()
        # Independent digest-bound finite authority exercises the same adoption
        # checks; production CLI has no override for its fixed known t0 digest.
        specification = dict(aggregate={'sha256': 'aggregate'}, manifest={'sha256': 'manifest'},
            preparation={'sha256': 'preparation'}, initial_ions={'sha256': 'ions'}, scalar='finite scalar',
            implementation=['old execution'], oracle={'sha256': 'old oracle'})
        payload = json.loads(state.path.read_text())
        payload['model_identity'] = sha256(json.dumps(specification, sort_keys=True).encode()).hexdigest()
        payload['evidence'] = dict(model=specification, original_raw_frames_unchanged=True)
        legacy_path = base/'legacy.json'
        legacy_path.write_text(json.dumps(payload))
        legacy_binding = file_binding(legacy_path)
        authority_path = base/'old-model.json'
        authority_path.write_text(json.dumps(dict(identity=payload['model_identity'],
            specification=specification, initial_boundary=legacy_binding)))
        authority_binding = file_binding(authority_path)
        changed_execution = dict(specification, implementation=['new execution'], oracle={'sha256': 'new oracle'})
        adoption_args = dict(specification=changed_execution, u=state.u, p=state.p,
            expected_slots=kwargs['expected_slots'], execution={'executor': 'isolated'},
            authority_binding=authority_binding)
        adopted = adopt_initial_boundary(legacy_binding, base/'adopted', **adoption_args)
        assert adopted.identity == physical_identity(specification)
        assert all(a.frame == b.frame and a.coefficient == b.coefficient
                   and a.metadata() == b.metadata() for a, b in zip(state.slots, adopted.slots))
        assert file_binding(legacy_path) == legacy_binding
        refusals = []
        for field in ('aggregate', 'manifest', 'preparation', 'initial_ions', 'scalar'):
            altered = dict(changed_execution, **{field: 'changed'})
            try:
                adopt_initial_boundary(legacy_binding, base/('changed-'+field),
                    **dict(adoption_args, specification=altered))
            except ValueError:
                refusals.append(field)
            else:
                raise AssertionError('changed physical specification adopted: '+field)
        for label, options in (('ions', dict(u=state.u+.01)),
                               ('coverage', dict(expected_slots=tuple(reversed(kwargs['expected_slots']))))):
            try:
                adopt_initial_boundary(legacy_binding, base/label, **dict(adoption_args, **options))
            except ValueError:
                refusals.append(label)
            else:
                raise AssertionError('changed initial state adopted: '+label)
        for label, change in (('nonzero-time', lambda p: p.update(time_ry_inverse=.1)),
                              ('changed-source', lambda p: p['slots'][0].update(source_identity='changed')),
                              ('changed-D', lambda p: p['slots'][0].update(coefficient=p['slots'][1]['coefficient']))):
            changed = json.loads(legacy_path.read_text())
            change(changed)
            target = base/(label+'.json')
            target.write_text(json.dumps(changed))
            try:
                adopt_initial_boundary(file_binding(target), base/(label+'-adopt'), **adoption_args)
            except ValueError:
                refusals.append(label)
            else:
                raise AssertionError('changed bound bytes adopted: '+label)
        # Production route must reject even an internally valid unrecognized t0.
        try:
            adopt_initial_boundary(legacy_binding, base/'unknown-production',
                **{k: v for k, v in adoption_args.items() if k != 'authority_binding'})
        except ValueError:
            refusals.append('unknown-production-t0')
        else:
            raise AssertionError('production adoption accepted an unrecognized t0')
    report = dict(schema='P254-streamed-evolution-oracle-v1', status='passed',
        scope='concurrency/checkpoint delta only; scientific checks reused without reexecution',
        prior_scientific_oracle=prior_binding, preserved_scientific_sources=preserved,
        serial_parallel_error=equation_error, permutation_error=permutation_error,
        completion_order=observed_completion_order, complete_stage_barriers=barrier_observations,
        exact_global_cap=3, failure_nonpublication=True, copy_free_adoption=True,
        adoption_refusals=refusals, worker_probes=probes,
        source_bindings=[file_binding(p) for p in (Path(__file__),
            ROOT/'src/substrate_framework/streamed_material_state.py',
            ROOT/'src/substrate_framework/streamed_material_evolution.py',
            ROOT/'src/substrate_framework/streamed_sector_executor.py',
            HERE/'streamed-evolution-registration.json', HERE/'streamed-sector-registration.json')])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    from contextlib import contextmanager
    from substrate_framework.prepared_reference import PreparedReference
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prior-oracle', type=Path, required=True)
    parser.add_argument('--prior-oracle-sha256', required=True)
    parser.add_argument('--worker-cpus', type=int, nargs='+', required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(73211)
    frames, ds, hs, gs = [], [], [], []
    for j in range(4):
        frames.append((rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2)))*.21)
        z = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
        ds.append(z @ z.conj().T*.12)
        z = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        hs.append((z+z.conj().T)*.22)
        z = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        gs.append((z+z.conj().T)*.14)
    k = np.array([[.7, -.11], [-.11, .51]])
    nbar = np.full((3, 2), .09)
    mass, curvature, ionic = 2.3, .8, .13

    class Coupling:
        npw, nr, ndof = 3, 3, 1
        grid_measure_bohr3, cell_volume_bohr3 = 1., 3.
        def __init__(self, index):
            self.index = index
        def moments(self, left, right, coefficient):
            weighted = left @ coefficient
            return np.sum(weighted*right.conj(), axis=1), np.array([
                np.sum(right.conj()*(gs[self.index] @ weighted))])
        def potential_action(self, field, vectors):
            return field[:, None]*vectors
        def ionic_action(self, direction, vectors):
            return direction[0]*(gs[self.index] @ vectors)

    @contextmanager
    def factory(slot):
        i = int(slot.identity)
        yield Coupling(i), lambda x: hs[i] @ x

    reference = PreparedReference(nbar, np.array([ionic]), np.zeros_like(nbar), np.zeros(1),
        np.array([[curvature]]), np.array([mass]), np.zeros(1), lambda n: n @ k, 2,
        'dense-scalar', 'dense-ionic', 'four-sectors-two-spin', 'explicit-zero',
        'one-coordinate', 'fixed-complex-D')
    kwargs = dict(model_identity='dense-two-spin-full-coverage',
                  expected_slots=tuple(str(i) for i in range(4)), grid_measure=1.)
    concurrency_checks(args, reference, factory, kwargs, frames, ds)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()