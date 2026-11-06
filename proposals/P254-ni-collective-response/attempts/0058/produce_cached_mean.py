"""Common mean/U/S update with paid non-Gamma force and Phi000 reuse."""
from pathlib import Path
import argparse
import json
import faulthandler
import os
import resource
from time import perf_counter
import numpy as np
import feedback_inputs
from mean_feedback_functional import SOURCE, CompactCubic, CommonMeanFunctional, binding, save, authenticate
from paid_mean_feedback_cache import PaidMeanFeedbackCache
from produce_mean_feedback import integrate
from mean_feedback_reference import MeanFeedbackReference
from material_electron_ion_inputs import load_material


def produce_cached_mean(dataset, baseline_catalog, cache_result, template_result, output,
                        *, rtol=2e-8, atol=2e-10):
    started = perf_counter()
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    def stage(name):
        entry = dict(stage=name, elapsed_seconds=perf_counter()-started)
        save(out/'progress.json', entry)
        print(json.dumps(entry), flush=True)
    faulthandler.dump_traceback_later(20, repeat=True)
    print(json.dumps(dict(stage='actual-cached-common-mean-start',
        pid=os.getpid(), affinity=sorted(os.sched_getaffinity(0)),
        cgroup=Path('/proc/self/cgroup').read_text(),
        threads={key: os.environ.get(key) for key in (
            'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})), flush=True)
    c, manifest, coefficients = load_material()
    stage('material-loaded')
    template = json.loads(Path(template_result).read_bytes())
    source = CompactCubic.from_source(SOURCE)
    stage('compact-source-loaded')
    cache = PaidMeanFeedbackCache(cache_result, coefficients)
    stage('paid-cache-loaded')
    try:
        nodes = cache.forcing_for(dataset, source)
        checkpoint = out/'forcing-checkpoint.npz'
        np.savez(checkpoint, times=nodes['times'], rho=nodes['rho'],
            ordered_gamma_covariance=nodes['ordered_gamma_covariance'],
            cubic_covariance_force=nodes['cubic_covariance_force'],
            gamma_force=nodes['gamma_force'])
        save(out/'forcing-checkpoint.json', dict(archive=binding(checkpoint),
            input_dataset=binding(dataset), cache=binding(cache_result),
            coefficient_manifest=coefficients,
            recovery='Complete node forcing only; no DOP853 iteration/history checkpoint.'))
        stage('updated-node-forcing-persisted')
        times, rho = nodes['times'], nodes['rho']
        if not np.array_equal(rho[0], c.rho0):
            raise ValueError('original prepared density changed')
        model = CommonMeanFunctional(c, source, phi=cache.phi000)
        stage('common-mean-model-ready')
        gradients = np.asarray([model.electronic_gradient(r) for r in rho])
        stage('electronic-gradients-ready')
        force = nodes['cubic_covariance_force']
        mean, u, s, work = integrate(model, times, rho, gradients, force,
                                     rtol=rtol, atol=atol)
        stage('common-mean-U-S-integration-complete')
        archive = out/'mean-feedback.npz'
        np.savez(archive, times=times, mean_gamma=mean, U=u, S_gamma=s,
            rho0=c.rho0, W0_gamma=c.mechanical.w0[0], node_rho=rho,
            node_ordered_gamma_covariance=nodes['ordered_gamma_covariance'],
            node_cubic_covariance_force=force, node_electronic_gradient=gradients)
        j = c.mechanical.j[0]
        defects = dict(
            initial_mean_error=float(np.max(np.abs(mean[0]))),
            initial_U_error=float(np.max(np.abs(u[0]-np.eye(c.m)))),
            initial_S_error=float(np.max(np.abs(s[0]-np.eye(model.nx)))),
            U_unitarity=[float(np.linalg.norm(v@v.conj().swapaxes(-1, -2)-np.eye(c.m))) for v in u],
            Gamma_S_symplectic=[float(np.linalg.norm(v@j@v.T-j)) for v in s],
            raw_rho_particle_trace=[float(np.trace(r, axis1=-2, axis2=-1).real.sum()) for r in rho])
        spec = dict(schema='p254-common-mean-feedback-reference-v1', status='complete',
            archive=binding(archive), coefficient_manifest=coefficients,
            input_dataset=binding(dataset), baseline_catalog=binding(baseline_catalog),
            physical_history_bindings=[*template['physical_history_bindings'], *cache.bindings],
            source_bindings=[*template['source_bindings'], binding(__file__),
                binding(Path(__file__).parents[1]/'0057/paid_mean_feedback_cache.py')],
            times_inverse_ry=times.tolist(),
            method=dict(solver='DOP853', rtol=rtol, atol=atol, segments=work),
            carried_boundary_oracle=template['oracle'], oracle_replayed=False,
            paid_force_cache=binding(cache_result),
            Gamma_node_contractions=nodes['Gamma_contractions'],
            nonGamma_node_contractions=nodes['nonGamma_contractions'],
            reference_source_kind='bare-mechanical-with-contacts',
            reference_source_detail=template['reference_source_detail'],
            functional=template['functional'], approximation=template['approximation'],
            bath_scope=template['bath_scope'], defects=defects,
            elapsed_seconds=perf_counter()-started,
            maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            limits=['One later mean update, not a converged coupled state.',
                    'Paid full force plus new-minus-paid Gamma uses disclosed binary64 arithmetic.',
                    'No other-q bath update, temporal enclosure, calibrated loss, heat or fullT30 claim.'])
        result_path = out/'reference-result.json'
        save(result_path, spec)
        reference = MeanFeedbackReference(c, result_path, coefficients)
        try:
            eg, el = reference.electron(len(times)-1, 0)
            dg, dl = reference.ordered(0, len(times)-1, 0, 'positions')
            consumer = dict(electron_shape=list(el.shape), Gamma_shape=list(dl.shape),
                electron_lesser_norm=float(np.linalg.norm(el)),
                Gamma_lesser_norm=float(np.linalg.norm(dl)),
                U_integration_calls=reference.u_integration_calls,
                Gamma_integration_calls=reference.gamma_integration_calls)
            save(out/'saved-consumer.json', consumer)
        finally:
            reference.close()
        result = dict(result=binding(result_path), archive=binding(archive),
            saved_consumer=binding(out/'saved-consumer.json'), defects=defects,
            Gamma_node_contractions=nodes['Gamma_contractions'],
            nonGamma_node_contractions=nodes['nonGamma_contractions'],
            elapsed_seconds=perf_counter()-started)
        print(json.dumps(dict(stage='actual-cached-common-mean-complete', **result)), flush=True)
        return result
    finally:
        cache.close()
        faulthandler.cancel_dump_traceback_later()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('dataset', 'baseline-catalog', 'cache-result', 'template-result', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    produce_cached_mean(a.dataset, a.baseline_catalog, a.cache_result, a.template_result, a.output)


if __name__ == '__main__':
    main()
