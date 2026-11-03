"""Run saved-coefficient Ni/H dynamics; retain requested states, not dense delta history.

Each completed requested time is a resumable finite-state boundary. Resuming
recomputes the small baseline interpolation; it does not restore DOP853 history.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0028')]
import numpy as np
from clamped_scientific_readout import encode
from matched_reference_io import save_arrays, runtime
from substrate_framework.streamed_material_state import file_binding, authenticate
from reduced_projection_cache import atomic_json, load_projected_coefficients
from reduced_electronic_dynamics import ActiveBandDynamics, DeltaState
from reduced_electronic_correlations import ActiveBandCorrelations
from reduced_material_readout import original_geometry, spatial_pulse, finite_kinetic_exchange


def store_arrays(path, arrays):
    return save_arrays(path, arrays, sum(np.asarray(v).nbytes for v in arrays.values())+(1 << 20))


def delta_arrays(state):
    return dict(rho=state.rho, **{'x_%02d' % q: x for q, x in enumerate(state.x)})


def restore_delta(binding, nq):
    with np.load(authenticate(binding), allow_pickle=False) as data:
        return DeltaState(data['rho'], tuple(data['x_%02d' % q] for q in range(nq)))


def interrupted(signum, frame):
    del frame
    raise InterruptedError('simulation interrupted by signal '+str(signum))


def configuration(args, manifest):
    return dict(schema='P254-activeband-material-run-configuration-v1',
                coefficient_specification=manifest['specification'],
                geometry_archive=file_binding(args.geometry_archive),
                bands=args.bands, features=args.features, times=args.times,
                atom=args.atom, axis=args.axis, amplitude=args.amplitude,
                baseline_only=args.baseline_only, correlations=args.correlations,
                correlation_reference=args.correlation_reference,
                rtol=args.rtol, atol=args.atol)


def run(args):
    started = perf_counter()
    runtime(args.output)
    manifest_bytes = args.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    config = configuration(args, manifest)
    identity = sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    result_path = args.output/'result.json'
    if result_path.exists():
        if not args.resume:
            raise ValueError('existing run requires explicit --resume')
        result = json.loads(result_path.read_text())
        if result['identity'] != identity:
            raise ValueError('resumed model, dose, geometry, times or tolerances differ')
        if result['status'] == 'complete':
            print('[dft-done] reduced material already complete '+str(result_path), flush=True)
            return
    else:
        result = dict(schema='P254-activeband-material-run-v1', identity=identity,
            configuration=config, status='running', completed=[],
            approximation='dynamical active saved bands plus explicit original passive background',
            full_PW_material_accuracy=False, whole_campaign_complete=False,
            resume_scope='finite requested-time state only; baseline and adaptive history are not reused')
        snapshot = args.output/'coefficient-manifest.json'
        with snapshot.open('xb') as stream:
            stream.write(manifest_bytes)
        result['original_coefficient_manifest_snapshot'] = file_binding(snapshot)
        result['original_input_manifest'] = file_binding(args.manifest)
    result['status'] = 'running'
    result['consumer_sources'] = [file_binding(HERE/name) for name in (
        'run_reduced_material_dynamics.py', 'reduced_material_readout.py',
        'reduced_electronic_dynamics.py', 'reduced_coefficients.py',
        'reduced_electronic_correlations.py', 'reduced_correlation_transpose.py',
        'reduced_nested_q_views.py',
        'reduced_projection_cache.py', 'reduced_projection_storage.py')]
    atomic_json(result_path, result)
    try:
        c = load_projected_coefficients(args.manifest, bands=args.bands,
            features=args.features, allow_partial=args.baseline_only)
        with np.load(authenticate(config['geometry_archive']), allow_pickle=False) as archive:
            geometry_values = {key: archive[key] for key in (
                'atom_indices', 'atom_names', 'cell_bohr', 'valence_charges')}
        geometry = original_geometry(geometry_values, c.mesh)
        atom = np.flatnonzero(geometry['ids'] == args.atom)
        if len(atom) != 1:
            raise ValueError('selected original atom is absent or duplicated')
        coordinate = 3*int(atom[0])+'xyz'.index(args.axis)
        mechanical_item = manifest['items']['mechanical']['archive']
        with np.load(authenticate(mechanical_item), allow_pickle=False) as data:
            field_mass = float(data['field_mass'])
        model = ActiveBandDynamics(c)
        baseline_started = perf_counter()
        baseline = model.integrate_baseline((0., args.times[-1]),
            rtol=args.rtol, atol=args.atol)
        result['baseline_recomputed_seconds'] = perf_counter()-baseline_started
        corr = ActiveBandCorrelations(c, baseline, rtol=args.rtol, atol=args.atol)
        completed = len(result['completed'])
        delta = model.initial_delta(coordinate=coordinate, amplitude=args.amplitude)
        previous = 0.
        if completed and not args.baseline_only:
            last = result['completed'][-1]
            delta = restore_delta(last['finite_state'], c.nq)
            previous = last['time_inverse_ry']
        for index in range(completed, len(args.times)):
            time = args.times[index]
            b = baseline.state(time)
            record = dict(time_inverse_ry=time)
            record['baseline_state'] = store_arrays(args.output/('baseline-%03d.npz' % index),
                                                    dict(rho=b.rho, x=b.x))
            if not args.baseline_only:
                if time > previous:
                    solved = model.integrate_finite((previous, time), delta,
                        baseline=baseline, t_eval=[time], dense_output=False,
                        rtol=args.rtol, atol=args.atol)
                    delta = solved.state(time)
                record['finite_state'] = store_arrays(args.output/('finite-%03d.npz' % index),
                                                       delta_arrays(delta))
                record['energy_change_ry'] = encode(model.energy_difference(b, delta).as_dict())
                record['spatial_pulse'] = encode(spatial_pulse(delta.x, geometry))
                record['finite_kinetic_exchange'] = encode(finite_kinetic_exchange(
                    b.x, delta.x, geometry, c.mechanical.masses,
                    geometry_values['valence_charges'], field_mass))
                if args.correlations:
                    equal = corr.spatial(time, time, geometry['cells'], [coordinate])
                    reference = min(time, args.correlation_reference)
                    two = corr.spatial(time, reference, geometry['cells'], [coordinate])
                    record['connected_correlation'] = encode(dict(
                        reference_time_inverse_ry=reference, noise_scope=equal.noise_scope,
                        equal_time=vars(equal), two_time=vars(two)))
                previous = time
            result['completed'].append(record)
            result['elapsed_seconds'] = perf_counter()-started
            atomic_json(result_path, result)
            print(json.dumps(dict(event='physical-time-complete', time_inverse_ry=time,
                result=str(result_path), mode='baseline' if args.baseline_only else 'finite-dose')), flush=True)
        result['status'] = 'complete'
    except BaseException as exc:
        result['status'] = 'failed'
        result['failure'] = dict(kind=type(exc).__name__, message=str(exc))
        raise
    finally:
        result['elapsed_seconds'] = perf_counter()-started
        atomic_json(result_path, result)
        print(('[dft-done]' if result['status'] == 'complete' else '[dft-fail]')+
              ' reduced material '+str(result_path), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--geometry-archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bands', type=int, choices=(8, 12), default=8)
    parser.add_argument('--features', type=int, choices=(7, 27), default=7)
    parser.add_argument('--times', default='0,1,10,30')
    parser.add_argument('--atom', type=int, default=1)
    parser.add_argument('--axis', choices=('x', 'y', 'z'), default='x')
    parser.add_argument('--amplitude', type=float, default=1.)
    parser.add_argument('--rtol', type=float, default=2e-8)
    parser.add_argument('--atol', type=float, default=2e-10)
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--correlations', action='store_true')
    parser.add_argument('--correlation-reference', type=float, default=0.)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    args.times = sorted(set(float(t) for t in args.times.split(',')))
    if not args.times or args.times[0] < 0 or not np.isfinite(args.times).all() or args.times[-1] <= 0:
        raise ValueError('finite nonnegative requested times with positive endpoint required')
    if args.baseline_only and args.correlations:
        raise ValueError('full spatial correlations require a complete all-q model')
    if args.correlation_reference < 0 or not np.isfinite(args.correlation_reference):
        raise ValueError('finite nonnegative correlation reference required')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    run(args)


if __name__ == '__main__':
    main()
