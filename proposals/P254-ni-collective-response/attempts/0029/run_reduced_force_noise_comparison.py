"""Compare the retained full-PW H37/x Gamma force pair with active-band noise.

Uses the existing nonstationary baseline, true tangent transpose and arbitrary
terminal-probe pullback. Total canonical force is the momentum-derivative row,
L(t).T e_p, matching joint_response_observables.total_force_probe. No original
PW trajectory or electronic initialization is rerun. This is an approximation
assessment, not a gate on the independent finite spatial simulation.
"""
import argparse
import json
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0028')]

import numpy as np
from clamped_scientific_readout import encode
from substrate_framework.streamed_material_state import authenticate, file_binding
from reduced_comparison_inputs import bound_json, physical_context, require, same_payload
from reduced_comparison_early import difference
from reduced_correlation_transpose import wick_edge_pair
from reduced_electronic_correlations import ActiveBandCorrelations, NOISE_SCOPE
from reduced_electronic_dynamics import ActiveBandDynamics
from reduced_projection_cache import atomic_json, load_projected_coefficients


def matched_inputs(args):
    manifest = bound_json(file_binding(args.manifest))
    intake = bound_json(file_binding(args.full_pw_intake))
    require(intake['schema'] == 'P254-original-material-force-covariance-complete-intake-v1',
            'retained original full-PW material force-pair intake required')
    common = intake['actual_consumer']['common']
    source = common['initialization']
    initialization = bound_json(dict(path=source[0], bytes=source[1], sha256=source[2]))
    saved = initialization['saved']
    specification = manifest['specification']
    require(saved['model_identity'] == common['model_identity']
            and common['case'] == specification['case'] == 'loaded-ni31h6'
            and common['observable'] == 'original-H37x-total-force'
            and common['atom'] == 37 and common['axis'] == 'x'
            and common['dynamics'] == 'screened-joint', 'original observable/model differs')
    require(physical_context(saved['specification']) == physical_context(specification['context']),
            'original and reduced physical contexts differ')
    require(same_payload(saved['original_boundary'], specification['original_boundary']),
            'original electronic preparation payload differs')
    authenticate(saved['original_boundary'])
    original = bound_json(intake['covariance_publication_original'])
    require(original['complete'] is True, 'original force-pair publication incomplete')
    reference = original['covariance']
    require(reference['time'] == common['times'][-1]
            and reference['second_time'] == common['second_times'][-1], 'original sample times differ')
    with np.load(authenticate(saved['preparation']), allow_pickle=False) as old, np.load(
            authenticate(manifest['items']['mechanical']['archive']), allow_pickle=False) as new:
        for source_key, target_key in (('masses_ry', 'masses'), ('local_ordered', 'ionic_ordered'),
                                      ('field_ordered_covariance', 'field_ordered'), ('omega_ry', 'trap_omega')):
            require(np.array_equal(old[source_key], new[target_key]), 'original mechanical preparation differs: '+source_key)
    with np.load(authenticate(file_binding(args.geometry_archive)), allow_pickle=False) as geometry:
        ids, names = geometry['atom_indices'], geometry['atom_names']
        selected = np.flatnonzero(ids == common['atom'])
        require(len(selected) == 1 and str(names[selected[0]]) == 'H', 'original H37 atom ordering differs')
        coordinate = 3*int(selected[0])
    return intake, reference, coordinate


def origins(correlation, left, right):
    c = correlation.c
    lx, ls, _ = correlation._split(0, left)
    rx, rs, _ = correlation._split(0, right)
    w = c.mechanical.w0[0]
    cut = 2*c.d
    ionic = (lx[:cut].T @ w[:cut, :cut] @ rx[:cut])[0, 0]
    field = (lx[cut:].T @ w[cut:, cut:] @ rx[cut:])[0, 0]
    cross = (lx[:cut].T @ w[:cut, cut:] @ rx[cut:]
             + lx[cut:].T @ w[cut:, :cut] @ rx[:cut])[0, 0]
    target = c.target[0]
    electronic = wick_edge_pair(ls, rs[target], c.rho0, c.rho0[target])[0, 0]
    return dict(electronic=electronic, ionic=ionic, field=field,
                initial_ionic_field_cross=cross, total=ionic+field+cross+electronic)


def reduced_pair(args, coordinate, time, second_time, bands, features):
    c = load_projected_coefficients(args.manifest, bands=bands, features=features, allow_partial=True)
    model = ActiveBandDynamics(c)
    baseline = model.integrate_baseline((0., max(time, second_time)), rtol=args.rtol, atol=args.atol)
    correlation = ActiveBandCorrelations(c, baseline, rtol=args.rtol, atol=args.atol)
    momentum = np.zeros((correlation.state_size(0), 1), complex)
    momentum[c.d+coordinate, 0] = 1
    left = correlation.backward_linear_probes(0, time, correlation.transpose(0, time, momentum))
    right = correlation.backward_linear_probes(0, second_time, correlation.transpose(0, second_time, momentum))
    ab, ba = origins(correlation, left, right), origins(correlation, right, left)
    return dict(bands=bands, features=features, AB=ab, BA=ba,
                symmetrized=(ab['total']+ba['total'])/2,
                commutator=ab['total']-ba['total'], noise_scope=NOISE_SCOPE)


def run(args):
    started = perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(dict(status='initializing'), stream)
    result = dict(schema='P254-activeband-original-force-noise-comparison-v1', status='running', rows=[],
        inputs=dict(manifest=file_binding(args.manifest), original_full_PW_intake=file_binding(args.full_pw_intake),
                    geometry_archive=file_binding(args.geometry_archive)),
        consumer_sources=[file_binding(HERE/name) for name in (
            'run_reduced_force_noise_comparison.py', 'reduced_comparison_inputs.py', 'reduced_comparison_early.py',
            'reduced_electronic_correlations.py', 'reduced_correlation_transpose.py', 'reduced_electronic_dynamics.py', 'reduced_ode_initial_step.py',
            'reduced_coefficients.py', 'reduced_projection_cache.py', 'reduced_projection_storage.py', 'reduced_nested_q_views.py')],
        limits=['Selected external Gamma H37/x force pair, not a local force or full spatial noise spectrum.',
                NOISE_SCOPE,
                'Raw differences assess this approximation at the existing early sample pair only; not error bounds, late-time calibration, entanglement or a lifetime.',
                'Both propagations retain their own nonstationary baseline; original full-PW discrete RK4 and reduced DOP853 errors are not separated by this comparison.'],
        units='(Ry/bohr)^2', comparison_sign='reduced minus retained original full-PW')
    try:
        intake, reference, coordinate = matched_inputs(args)
        time, second_time = reference['time'], reference['second_time']
        result.update(original_publication=intake['covariance_publication_original'],
                      original_preparation_identity=reference['preparation_identity'],
                      q_fractional=[0., 0., 0.], atom=37, axis='x', time_inverse_ry=time,
                      second_time_inverse_ry=second_time, original_full_PW=reference,
                      rtol=args.rtol, atol=args.atol)
        atomic_json(args.output, encode(result))
        for bands, features in ((8, 7), (8, 27), (12, 7), (12, 27)):
            row = reduced_pair(args, coordinate, time, second_time, bands, features)
            row['difference_by_origin'] = {order: {key: difference(
                np.asarray(complex(*reference[order][key])), np.asarray(row[order][key]))
                for key in ('electronic', 'ionic', 'field', 'total')} for order in ('AB', 'BA')}
            result['rows'].append(row)
            atomic_json(args.output, encode(result))
            print(json.dumps(dict(event='actual-force-pair-complete', bands=bands, features=features,
                                  result=str(args.output))), flush=True)
        result['status'] = 'complete'
    except BaseException as exc:
        result.update(status='failed', failure=dict(kind=type(exc).__name__, message=str(exc)))
        raise
    finally:
        result.update(elapsed_seconds=perf_counter()-started,
                      maxrss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        atomic_json(args.output, encode(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--full-pw-intake', type=Path, required=True)
    parser.add_argument('--geometry-archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rtol', type=float, default=2e-8)
    parser.add_argument('--atol', type=float, default=2e-10)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
