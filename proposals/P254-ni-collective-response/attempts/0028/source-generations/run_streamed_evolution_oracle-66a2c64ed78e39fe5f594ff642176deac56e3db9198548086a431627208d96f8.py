#!/usr/bin/env python3
"""Independent finite dense oracle; Root executes, never a material result."""
import argparse
import json
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
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


def main():
    from contextlib import contextmanager
    from substrate_framework.prepared_reference import PreparedReference
    from substrate_framework.streamed_material_evolution import StreamedMeanEvolution, exponential_action
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
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

    def moments(xs):
        density, gradient = np.zeros((3, 2), complex), 0j
        for i, (x, d) in enumerate(zip(xs, ds)):
            rho = x @ d @ x.conj().T
            density[:, i % 2] += .5*np.diag(rho)
            gradient += .5*np.trace(rho @ gs[i])
        return density, gradient

    def scalar(xs, u, p):
        n, g = moments(xs)
        delta = n-nbar
        energy = sum(.5*np.trace(x @ d @ x.conj().T @ h)
                     for x, d, h in zip(xs, ds, hs))
        return (energy + .5*np.sum(delta*(delta @ k))
                + u*(g-ionic) + .5*curvature*u*u + p*p/(2*mass)).real

    def pack(xs, u, p):
        return np.concatenate([x.ravel() for x in xs]+[np.array([u, p])])

    def unpack(y):
        return [y[6*i:6*(i+1)].reshape(3, 2) for i in range(4)], y[-2], y[-1]

    def rhs(t, y):
        xs, u, p = unpack(y)
        n, g = moments(xs)
        v = (n-nbar) @ k
        dx = [-1j*(hs[i]+np.diag(v[:, i % 2])+u*gs[i]) @ x for i, x in enumerate(xs)]
        return pack(dx, p/mass, ionic-g-curvature*u)

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
    model = StreamedMeanEvolution(reference, factory, **kwargs)
    initial = pack(frames, .07, -.03)
    interval = .3
    exact = solve_ivp(rhs, (0., interval), initial, method='DOP853',
                      rtol=2e-13, atol=2e-14)
    assert exact.success
    oracle_final = exact.y[:, -1]
    eps = 2e-6
    fd_force = -(scalar(frames, .07+eps, -.03)-scalar(frames, .07-eps, -.03))/(2*eps)
    assert abs(fd_force-rhs(0., initial)[-1]) < 1e-9
    direction = rng.normal(size=(3, 2))
    delta = moments(frames)[0]-nbar
    quadratic = lambda x: .5*np.sum(x*(x @ k))
    fd_kernel = (quadratic(delta+eps*direction)-quadratic(delta-eps*direction))/(2*eps)
    assert abs(fd_kernel-np.sum(direction*(delta @ k))) < 1e-9
    frozen_error = np.linalg.norm(exponential_action(lambda x: hs[0] @ x, frames[0], .23)
                                  - expm(-.23j*hs[0]) @ frames[0])
    assert frozen_error < 2e-12
    with tempfile.TemporaryDirectory(prefix='P254-stream-oracle-') as temporary:
        base = Path(temporary)
        state = state_check(base/'initial', frames, ds)
        evaluation = model.reduce(state, energy=True)
        assert abs(evaluation.energy-scalar(frames, .07, -.03)) < 1e-12
        assert abs(evaluation.force[0]-fd_force) < 1e-9
        assert abs(evaluation.number-sum(.5*np.trace(x @ d @ x.conj().T)
                   for x, d in zip(frames, ds))) < 1e-12
        capped = StreamedMeanEvolution(reference, factory, max_h0_actions=1, **kwargs)
        try:
            capped.step(state, base/'capped', dt=.03, max_iterations=40, tolerance=2e-13)
        except RuntimeError as error:
            assert 'action cap' in str(error)
        else:
            raise AssertionError('H0 work cap was not enforced')
        assert not (base/'capped'/'boundary.json').exists()
        try:
            model.step(state, base/'unconverged', dt=.3, max_iterations=1, tolerance=1e-15)
        except RuntimeError as error:
            assert 'failed to converge' in str(error)
        else:
            raise AssertionError('nonconverged stage was published')
        assert not (base/'unconverged'/'boundary.json').exists()
        errors, energies, gram_errors, number_errors = [], [], [], []
        def run(consumer, steps, name):
            final, history = consumer.evolve(state, base/name, physical_interval=interval,
                steps=steps, max_steps=32, max_iterations=40, tolerance=2e-13,
                refinement_label=name)
            xs = [np.array(final.arrays(s)[0]) for s in final.slots]
            return final, xs, history
        for count in (4, 8, 16):
            final, xs, history = run(model, count, 'n'+str(count))
            errors.append(float(np.linalg.norm(pack(xs, final.u[0], final.p[0])-oracle_final)))
            energies.append(float(max(abs(row['energy_ry']-evaluation.energy) for row in history)))
            gram_errors.append(float(max(np.linalg.norm(x.conj().T @ x-y.conj().T @ y)
                                         for x, y in zip(xs, frames))))
            number_errors.append(float(abs(sum(.5*np.trace(x @ d @ x.conj().T)
                                      for x, d in zip(xs, ds))-evaluation.number)))
        assert errors[0]/errors[1] > 3.3 and errors[1]/errors[2] > 3.3
        assert energies[0]/energies[1] > 2.5 and energies[1]/energies[2] > 2.5
        assert max(gram_errors+number_errors) < 2e-10

        class WrongForce(StreamedMeanEvolution):
            def force(self, gradient, u):
                return self.reference.ionic_force+gradient-self.reference.curvature @ u
        class WrongKernel(StreamedMeanEvolution):
            @staticmethod
            def action(c, h, s, stage, v):
                return h(v)-c.potential_action(stage.potential[:, s.component], v)+c.ionic_action(stage.u, v)
        class MissingWeight(StreamedMeanEvolution):
            def stage(self, density, gradient, u, **options):
                return super().stage(2*density, 2*gradient, u, **options)
        class StaleStage(StreamedMeanEvolution):
            def reduce(self, state, **options):
                if options.get('frozen_stage') is not None:
                    options['frozen_stage'] = None
                return super().reduce(state, **options)
        mutations = {}
        for cls in (WrongForce, WrongKernel, MissingWeight, StaleStage):
            mutated = cls(reference, factory, **kwargs)
            final, xs, _ = run(mutated, 16, cls.__name__)
            error = float(np.linalg.norm(pack(xs, final.u[0], final.p[0])-oracle_final))
            assert error > 20*errors[-1], (cls.__name__, error, errors[-1])
            mutations[cls.__name__] = dict(consumer_error=error, detected=True)
        payload = json.loads(state.path.read_text())
        payload['slots'].reverse()
        reversed_path = base/'reversed.json'
        reversed_path.write_text(json.dumps(payload))
        reversed_state = DiskMean(reversed_path)
        reverse_model = StreamedMeanEvolution(reference, factory, **dict(
            kwargs, expected_slots=tuple(reversed(kwargs['expected_slots']))))
        reverse_final = reverse_model.step(reversed_state, base/'reverse-step', dt=interval/4,
                                           max_iterations=40, tolerance=2e-13)
        forward_final = DiskMean(base/'n4'/'step-000001'/'boundary.json')
        reverse_frames = {s.identity: reverse_final.arrays(s)[0] for s in reverse_final.slots}
        permutation_error = max(float(np.linalg.norm(forward_final.arrays(s)[0]-reverse_frames[s.identity]))
                                for s in forward_final.slots)
        assert permutation_error < 2e-12
        assert np.linalg.norm(reverse_final.u-forward_final.u) < 2e-12
        assert model.counts['h0_calls'] > 0 and model.counts['moment_reductions'] > 0
        assert final.evidence['work_counts']['h0_calls'] > 0
    report = dict(schema='P254-streamed-evolution-oracle-v1', status='passed',
        scope='independent finite dense oracle only; no actual material execution',
        physical_interval_ry_inverse=interval, steps=[4, 8, 16], state_errors=errors,
        energy_errors=energies, gram_errors=gram_errors, number_errors=number_errors,
        frozen_exponential_error=frozen_error, stream_permutation_error=permutation_error,
        mutations=mutations, bounded_action_rejection=True, nonlinear_rejection=True,
        observed_work_counts=model.counts)
    from substrate_framework.streamed_material_state import file_binding
    report['source_bindings'] = [file_binding(p) for p in (Path(__file__),
        ROOT/'src/substrate_framework/streamed_material_state.py',
        ROOT/'src/substrate_framework/streamed_material_evolution.py',
        HERE/'streamed-evolution-registration.json')]
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()