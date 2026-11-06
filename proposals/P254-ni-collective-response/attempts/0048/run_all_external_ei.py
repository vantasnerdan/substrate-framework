"""All external original moving E-I kernels, with bounded nonduplicating output."""
from pathlib import Path
import argparse
import json
import os
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / 'src')] + [str(HERE.parent / name) for name in
    ('0041', '0029', '0028', '0027', '0025', '0024')]
import numpy as np
from run_material_electron_ion_kernel import _compute_q, MATCHING
from material_electron_ion_inputs import (DEFAULT_MANIFEST, SavedTransitions,
    load_material, reference_diagnostics)
from substrate_framework.streamed_material_state import authenticate, file_binding
from streamed_ei_output import CHANNELS, EIOutput


def source_bindings():
    return [file_binding(HERE / name) for name in
            ('run_all_external_ei.py', 'streamed_ei_output.py')] + [
        file_binding(HERE.parent / '0041' / name) for name in
        ('run_material_electron_ion_kernel.py', 'material_electron_ion_inputs.py',
         'electron_ion_functional.py')]


def own_banner():
    membership = Path('/proc/self/cgroup').read_text()
    own = next(line.split(':', 2)[2] for line in membership.splitlines()
               if line.startswith('0::'))
    root = Path('/sys/fs/cgroup') / own.lstrip('/')
    limits = {}
    for name in ('memory.max', 'memory.swap.max', 'cpuset.cpus.effective'):
        try:
            limits[name] = dict(status='available', value=(root / name).read_text().strip())
        except OSError as error:
            limits[name] = dict(status='unavailable', reason=str(error))
    return dict(stage='actual-own-native-start', pid=os.getpid(),
                affinity=sorted(os.sched_getaffinity(0)), cgroup=membership,
                limits=limits, thread_environment={key: os.environ.get(key) for key in
                    ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')})


def compute_streamed(c, reference, store, progress):
    completed = store.consume()
    if completed is not None:
        return completed, dict(new_q=0, reused_q=store.nq, new_edge_pairs=0)
    sigma, done, last = store.restore()
    original_done = done
    electrons = {(a, b): reference.electron(a, b)
                 for a in range(len(reference.times)) for b in range(len(reference.times))}
    actual_pairs, new_q = 0, 0
    defects = dict(vertex_paired_adjoint_max=0., boson_ordered_ccr_max=0.,
                   boson_paired_symplectic_max=0.)
    for q in range(done, store.nq):
        started = perf_counter()
        have_pi = store.pi.has(store.key(q))
        arrays, count, q_defects = _compute_q(c, reference, q, list(range(c.ns)),
                                             not have_pi, electrons)
        for key in CHANNELS:
            sigma[key] += arrays['sigma_' + key]
        pi = (store.pi_values(q) if have_pi else
              {key: arrays['pi_' + key] for key in CHANNELS})
        for key in defects:
            defects[key] = max(defects[key], q_defects[key])
        last = store.commit(sigma, q + 1, pi,
                            dict(q_index=q, edge_pairs=count, diagnostics=q_defects))
        actual_pairs += count
        new_q += 1
        progress(dict(stage='full-external-EI-q-complete', q=q, completed_q=q+1,
                      all_q=store.nq, actual_edge_pairs=count, reused_pi=have_pi,
                      elapsed_seconds=perf_counter()-started,
                      maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    sigma_norms = {key: float(np.linalg.norm(value)) for key, value in sigma.items()}
    pi_norms = [{key: float(np.linalg.norm(value)) for key, value in store.pi_values(q).items()}
                for q in range(store.nq)]
    store.publish(last, dict(sigma_norms=sigma_norms, pi_norms=pi_norms,
        analytic_edge_pair_evaluations=len(reference.times)**2*c.ns*store.nq,
        new_edge_pair_evaluations=actual_pairs, new_q=new_q, reused_q=original_done,
        diagnostics=defects, matching=MATCHING,
        recovery_scope='committed original-order Sigma prefix and immutable Pi-q; current q only may repeat',
        selfconsistent=False, lifetime_or_rate_computed=False, full_material_goal_complete=False))
    return store.consume(), dict(new_q=new_q, reused_q=original_done,
                                new_edge_pairs=actual_pairs)


def make_store(output, specification):
    return EIOutput(output, specification, specification['sigma_shape'],
                    specification['pi_shape'], specification['nq'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path(DEFAULT_MANIFEST))
    parser.add_argument('--saved-transitions', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError('numerical output must stay outside repository')
    print(json.dumps(own_banner()), flush=True)
    started = perf_counter()
    sources = source_bindings()
    if args.resume and (output / 'result.json').exists():
        specification = json.loads((output / 'specification.json').read_bytes())
        if (specification['coefficient_manifest'] != file_binding(args.manifest)
                or specification['transition_bundle'] != file_binding(args.saved_transitions)
                or specification['source_bindings'] != sources):
            raise ValueError('complete output original inputs/implementation differ')
        for binding in specification['transition_bindings']:
            authenticate(binding)
        result = make_store(output, specification).consume()
        if result is None:
            raise ValueError('complete saved output missing')
        print(json.dumps(dict(stage='complete-saved-output-consumed', new_edge_pairs=0,
              result=file_binding(output / 'result.json'), elapsed_seconds=perf_counter()-started)), flush=True)
        return
    output.mkdir(parents=True, exist_ok=args.resume)
    c, manifest, binding = load_material(args.manifest)
    times = [0., 1., 10., 30.]
    reference = SavedTransitions(c, times, args.saved_transitions, binding)
    specification = dict(kind='full-external-original-moving-electron-ion',
        coefficient_manifest=binding, transition_bundle=file_binding(args.saved_transitions),
        transition_bindings=reference.bindings, source_bindings=sources,
        times_inverse_ry=times, nq=len(c.q_indices), sigma_external_slots=list(range(c.ns)),
        pi_external_q=list(range(len(c.q_indices))),
        sigma_shape=[len(times), len(times), c.ns, c.m, c.m],
        pi_shape=[len(times), len(times), c.d, c.d], original_boundary=manifest['specification']['original_boundary'])
    try:
        store = make_store(output, specification)
        result, counts = compute_streamed(c, reference, store,
                                         lambda value: print(json.dumps(value), flush=True))
        print(json.dumps(dict(stage='full-external-EI-complete', counts=counts,
            result=file_binding(output / 'result.json'),
            elapsed_seconds=perf_counter()-started,
            maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            original_reference_diagnostics=reference_diagnostics(reference))), flush=True)
    finally:
        reference.close()


if __name__ == '__main__':
    main()
