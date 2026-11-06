"""Actual G1 full31 common-mean feedback and its saved public consumption.

The immutable all64 force bank supplies the unchanged non-Gamma bath. Only its
Gamma contribution is replaced using the authenticated raw G1 equal-time data.
The unchanged integrator has no persistent in-flight checkpoint: interruption
before it returns loses this mean integration, not the paid G1/cache work.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)] + [str(HERE.parent / n) for n in (
    '0074', '0063', '0062', '0060', '0057', '0055', '0052', '0051',
    '0050', '0047', '0046', '0044', '0045', '0041', '0030', '0029')]
from dressed_mean import DressedMean, SCHEMA, REFERENCE
from fullwindow_force_cache import FullWindowForceCache
from saved_mean import SavedMean, POSITIONS
from history_storage import authenticate, binding, atomic_json
from mean_feedback_functional import CommonMeanFunctional, CompactCubic, SOURCE
from produce_mean_feedback import integrate
from material_electron_ion_inputs import load_material
from history_readout import gamma_readout, json_value


def consume(result_path, model):
    """Read the new saved reference into physical observables and the functional."""
    started = time.monotonic()
    with DressedMean(result_path) as saved:
        nt = len(saved.times)
        trace = np.zeros(nt, dtype=complex)
        retarded_maxabs = 0.
        for slot in range(saved.rho0.shape[0]):
            ret, less = saved.electron_slot(slot)
            if not np.isfinite(ret).all() or not np.isfinite(less).all():
                raise ValueError('nonfinite saved electron reference')
            trace += np.trace(-1j*less[np.arange(nt), np.arange(nt)],
                              axis1=-2, axis2=-1)
            retarded_maxabs = max(retarded_maxabs, float(np.abs(ret).max()))
        gr, gl = saved.gamma_slot()
        if not np.isfinite(gr).all() or not np.isfinite(gl).all():
            raise ValueError('nonfinite saved Gamma reference')
        gamma = gamma_readout(gr, gl, int(saved.n))
        scalars, rhs_norms, curvature_norms = [], [], []
        for a in range(nt):
            scalar = model.scalar(saved.node_rho[a], saved.mean[a],
                                  saved.node_cubic_covariance_force[a])
            rhs = model.mean_rhs(saved.mean[a], saved.node_electronic_gradient[a],
                                 saved.node_cubic_covariance_force[a])
            curvature = model.curvature(saved.mean[a])
            if not all(np.isfinite(value).all() for value in (scalar, rhs, curvature)):
                raise ValueError('nonfinite saved common-functional output')
            scalars.append(dict(real=float(scalar.real), imag=float(scalar.imag)))
            rhs_norms.append(float(np.linalg.norm(rhs)))
            curvature_norms.append(float(np.linalg.norm(curvature)))
        findings = dict(result=saved.binding, archive=saved.spec['archive'],
            reader=binding(HERE/'dressed_mean.py'),
            consumer_implementation=binding(Path(__file__)),
            input_identity=saved.input_identity, preparation_origin=saved.preparation_origin,
            electron_sectors_consumed=int(saved.rho0.shape[0]),
            reference_fermion_trace=[dict(real=float(x.real), imag=float(x.imag)) for x in trace],
            reference_electron_retarded_maxabs=retarded_maxabs,
            Gamma_shape=list(gl.shape),
            reference_Gamma=json.loads(json.dumps(gamma, default=json_value)),
            final_mean_norm=float(np.linalg.norm(saved.mean[-1])),
            common_scalar_at_saved_nodes=scalars, saved_mean_rhs_norms=rhs_norms,
            saved_curvature_norms=curvature_norms, times=saved.times.tolist(),
            elapsed_seconds=time.monotonic()-started,
            scope='Actual saved Gaussian mean reference and common functional consumed. Held forcing remains raw dressed input; reference histories start from original preparation. This does not repair that input, advance the next Dyson iterate, or establish temporal accuracy/conservation.')
    consumer_path = Path(result_path).parent/'saved-consumer.json'
    atomic_json(consumer_path, findings)
    return dict(consumer=binding(consumer_path), findings=findings)


def produce(history, force_cache, output, *, historical_source_snapshots,
            rtol=2e-8, atol=2e-10):
    started = time.monotonic()
    history_binding, force_binding = binding(history), binding(force_cache)
    provider = json.loads(authenticate(history_binding).read_bytes())
    force_record = json.loads(authenticate(force_binding).read_bytes())
    reference_binding = provider['input_identity']['numerical_reference']
    preparation_binding = force_record['original_input_dataset']
    authenticate(reference_binding)
    authenticate(preparation_binding)
    snapshots = json.loads(Path(historical_source_snapshots).read_bytes())['records']
    for item in (provider['input_mean'], provider['input_kernels']):
        authenticate(item)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    c, _, coeff = load_material()
    original = json.loads(authenticate(reference_binding).read_bytes())
    source = CompactCubic.from_source(SOURCE)
    reader = binding(HERE/'dressed_mean.py')
    reader_snapshot = output/'dressed_mean_producer.py.snapshot'
    reader_snapshot.write_bytes(authenticate(reader).read_bytes())
    with SavedMean(authenticate(preparation_binding)) as bare:
        if (not bare.is_bare or bare.reference_identity != reference_binding or
                not np.array_equal(bare.rho0, c.rho0) or
                not np.array_equal(bare.w0, c.mechanical.w0[0]) or
                int(bare.archive['n']) != c.n or
                (c.n, c.nq, c.ns, c.m, c.d, c.f) != (64, 64, 128, 12, 111, 3)):
            raise ValueError('original preparation/material/reference differs')
    with FullWindowForceCache(force_cache, coeff,
            historical_source_snapshots=snapshots) as cache:
        if cache.spec['original_input_dataset'] != preparation_binding:
            raise ValueError('paid cache original preparation identity differs')
        actual = cache.forcing_for(history, source)
        times = np.asarray(actual['times']).copy()
        rho = actual['rho']
        w = actual['ordered_gamma_covariance']
        force = actual['cubic_covariance_force']
        if (actual['input_dataset'] != history_binding or
                not np.array_equal(times, np.arange(31, dtype=float)) or
                not np.array_equal(rho[0], c.rho0) or
                not np.array_equal(w[0], c.mechanical.w0[0][np.ix_(POSITIONS, POSITIONS)]) or
                original['coefficient_manifest'] != coeff or
                actual['nonGamma_contractions'] != 0 or actual['Gamma_contractions'] != 31):
            raise ValueError('actual G1 all31 forcing/time/preparation differs')
        model = CommonMeanFunctional(c, source, phi=cache.phi000)
        gradients = np.asarray([model.electronic_gradient(r) for r in rho])
        print('A86_ACTUAL_G1_INTEGRATION_BEGIN', flush=True)
        mean, u, s, work = integrate(model, times, rho, gradients, force, rtol=rtol, atol=atol)
        arrays = dict(times=times, mean_gamma=mean, U=u, S_gamma=s, rho0=c.rho0,
            W0_gamma=c.mechanical.w0[0], n=np.asarray(c.n, dtype=np.int64),
            node_rho=rho, node_ordered_gamma_covariance=w,
            node_cubic_covariance_force=force, node_electronic_gradient=gradients)
        archive = output/'dressed-mean.npz'
        np.savez(archive, **arrays)
        j = c.mechanical.j[0]
        defects = dict(
            initial_mean_max=float(np.max(np.abs(mean[0]))),
            initial_U_max=float(np.max(np.abs(u[0]-np.eye(12)))),
            initial_S_max=float(np.max(np.abs(s[0]-np.eye(228)))),
            initial_rho_max=float(np.max(np.abs(rho[0]-c.rho0))),
            initial_ordered_covariance_max=float(np.max(np.abs(
                w[0]-c.mechanical.w0[0][np.ix_(POSITIONS, POSITIONS)]))),
            U_unitarity=[float(np.linalg.norm(v@v.conj().swapaxes(-1, -2)-np.eye(12))) for v in u],
            Gamma_S_symplectic=[float(np.linalg.norm(v@j@v.T-j)) for v in s],
            Gamma_hessian_transpose=[float(np.linalg.norm(model.curvature(x)-model.curvature(x).T)) for x in mean],
            raw_rho_hermiticity=[float(np.linalg.norm(r-r.conj().swapaxes(-1, -2))) for r in rho],
            raw_rho_particle_trace=[float(np.trace(r, axis1=-2, axis2=-1).real.sum()) for r in rho],
            mean_norm=[float(np.linalg.norm(x)) for x in mean],
            covariance_force_norm=[float(np.linalg.norm(f)) for f in force])
        paths = [Path(__file__), reader_snapshot]
        paths += [HERE.parent/d/n for d, n in (
            ('0057', 'mean_feedback_functional.py'), ('0057', 'produce_mean_feedback.py'),
            ('0029', 'reduced_electronic_dynamics.py'), ('0047', 'cubic_feedback.py'),
            ('0051', 'cubic_contour.py'), ('0074', 'fullwindow_force_cache.py'),
            ('0063', 'saved_mean.py'), ('0063', 'history_provider.py'),
            ('0063', 'history_storage.py'), ('0041', 'material_electron_ion_inputs.py'))]
        spec = dict(schema=SCHEMA, status='complete', reference_source_kind=REFERENCE,
            input_identity=dict(actual_dataset=actual['input_dataset'],
                numerical_reference=reference_binding),
            preparation_origin=preparation_binding, archive=binding(archive),
            coefficient_manifest=coeff, force_cache=force_binding,
            compact_source=source.source_binding,
            compact_source_arrays=source.source_record['arrays'],
            physical_history_bindings=[*original['physical_history_bindings'],
                reference_binding, preparation_binding,
                provider['input_mean'], provider['input_kernels'],
                *actual['input_bindings'], *actual['cache_bindings'],
                source.source_binding, source.source_record['arrays']],
            source_bindings=[binding(p) for p in paths],
            producer_reader=dict(original=reader, preserved=binding(reader_snapshot),
                interpretation='Actual DressedMean implementation used for this saved public consumer; exact production bytes preserved.'),
            times_inverse_ry=times.tolist(),
            method=dict(solver='DOP853', rtol=rtol, atol=atol, segments=work),
            force_work=dict(nonGamma_contractions=actual['nonGamma_contractions'],
                Gamma_contractions=actual['Gamma_contractions'],
                completed_cache_internal_q=list(range(64)),
                Phi000='Original authenticated stored cache tensor; no recontraction',
                historical_source_snapshots=snapshots,
                cache_source_authentication=cache.source_authentication,
                cache_reader=cache.current_reader_binding,
                retained_paid_force_cache=force_binding,
                retained_paid_cache_elapsed_seconds=cache.spec['elapsed_seconds'],
                retained_dressed_history=history_binding),
            initial_preparation='Literal original rho0/W0/n and unchanged initial_baseline mean; U0=S0=I. Actual forcing input is G1, not this preparation.',
            defects=defects, elapsed_seconds=time.monotonic()-started,
            approximation='Existing first-order interpolation of actual G1 forcing nodes only (raw rho, electronic gradient, all64 contracted covariance force). No U/S or physical-kernel interpolation. Forcing-time error unresolved; not temporal accuracy.',
            functional='Unchanged CommonMeanFunctional/electronic_gradient/integrate with original CompactCubic/Phi000, all128/12/fullphase228 and original scalar/contacts/passive complement; original nonGamma bath retained.',
            recovery=dict(checkpoint_granularity='Whole mean case only: integrate returns mean/U/S and segment work only after all30 intervals finish; no persistent in-flight checkpoint or resume.',
                immutable_reusable_inputs='Completed G1 histories and paid all64/full31 force bank remain reusable without replay.',
                manager_cap_source='Actual Root allocation and manager receipt; no producer-internal reset.'),
            limits=[
                'One actual changed common mean/U/Gamma-S feedback step from completed G1 nodes0..30; not a collision kernel, full state or fixed point.',
                'No fullmesh, physical calibration, rate, conservation, temporal accuracy, ranking, loss, heat or closed all-q feedback claim.',
                'Raw G1 number drift and all numerical defects retained; no clipping, refill, mode reduction or invariant repair.',
                'Historical MeanFeedbackReference requires a distinct schema and baseline catalog; none is invented. DressedMean is the supported saved all128/fullGamma consumer.'])
        result_path = output/'result.json'
        atomic_json(result_path, spec)
        consumed = consume(result_path, model)
    result = dict(result=binding(result_path), archive=binding(archive),
        consumer=consumed['consumer'], findings=defects,
        consumer_findings=consumed['findings'], elapsed_seconds=time.monotonic()-started)
    atomic_json(output/'handoff.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', required=True, type=Path)
    parser.add_argument('--force-cache', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--historical-source-snapshots', required=True, type=Path)
    args = parser.parse_args()
    result = produce(args.history, args.force_cache, args.output,
                     historical_source_snapshots=args.historical_source_snapshots)
    print('A86_DRESSED_MEAN_COMPLETE '+json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
