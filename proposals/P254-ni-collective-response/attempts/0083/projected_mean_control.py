"""Actual coupled Gaussian endpoint control at the archived A78 mean nodes."""
from pathlib import Path
import argparse
import json
import os
import resource
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/name) for name in
    ('0063', '0074', '0057', '0069', '0050', '0047', '0052', '0041', '0030', '0029')]
import numpy as np
from history_storage import authenticate, binding, atomic_json
from saved_mean import SavedMean
from fullwindow_force_cache import FullWindowForceCache
from material_electron_ion_inputs import load_material
from mean_feedback_functional import CommonMeanFunctional, CompactCubic
from projected_gamma_control import jsonable, save_json
from reduced_coherence_readout import normalized_ordered
from mean_matched_tangent import CoupledMeanControl

SCHEMA = 'p254-node-mean-coupled-gamma-control-v1'


def consume(path):
    """Read the persisted numerical covariance through the original normalizer."""
    path = Path(path).resolve()
    record = json.loads(path.read_bytes())
    if record['schema'] != SCHEMA or record['status'] != 'complete-at-stated-scope':
        raise ValueError('actual completed node-mean coupled control required')
    with np.load(authenticate(record['raw_archive']), allow_pickle=False) as raw:
        covariance, times = raw['covariance'], raw['times']
        if (covariance.shape != (2, 2) or covariance.dtype != np.dtype('<c16') or
                not np.array_equal(times, [0., 30.]) or not np.isfinite(covariance).all()):
            raise ValueError('literal finite complex128 endpoint covariance required')
        normalized = normalized_ordered(covariance[1, 0], covariance[1, 1], covariance[0, 0])
    if jsonable(normalized) != record['normalized_ordered']:
        raise ValueError('saved numerical normalization differs from produced control')
    return dict(schema='p254-a83-saved-control-consumer-v1', control=binding(path),
        probe=record['probe'], times=times.tolist(), raw_complex_covariance=covariance,
        normalized_ordered=normalized, input_mean=record['input_mean'],
        scope=record['approximation'])


def produce(case):
    started = time.monotonic()
    output = Path(case['output']).resolve()
    output.mkdir(parents=True, exist_ok=False)
    c, _manifest, coefficient_binding = load_material()
    restoration = json.loads(authenticate(case['cache_restoration']).read_bytes())['records']
    with SavedMean(authenticate(case['mean'])) as mean:
        if (not mean.is_full_mean or mean.spec['coefficient_manifest'] != coefficient_binding or
                mean.reference_identity != case['numerical_reference'] or
                not np.array_equal(mean.rho0, c.rho0) or
                not np.array_equal(mean.w0, c.mechanical.w0[0])):
            raise ValueError('literal fullhistory mean, original basis/reference/preparation required')
        rho = mean.archive['node_rho']
        if rho.shape != (31, 128, 12, 12) or not np.array_equal(rho[0], c.rho0):
            raise ValueError('actual original held density nodes required')
        source = CompactCubic.from_source(authenticate(mean.spec['compact_source']))
        with FullWindowForceCache(authenticate(mean.spec['force_cache']), coefficient_binding,
                historical_source_snapshots=restoration) as cache:
            functional = CommonMeanFunctional(c, source, phi=cache.phi000)
            model = CoupledMeanControl(c, functional, mean.times, rho, mean.mean,
                rtol=2e-8, atol=2e-10, method='DOP853')
            probe = model.ionic_probes(0, (0,))/c.sqrt_n
            initial = model.backward_linear_probes(0, 0., probe)
            print('A83_ACTUAL_COUPLED_BACKWARD_BEGIN', flush=True)
            pulled = model.backward_linear_probes(0, 30., probe)
            duals = (initial, pulled)
            covariance = np.asarray([[model.initial_ordered(0, left, right).item()
                for right in duals] for left in duals], dtype=np.complex128)
            expected_initial = c.mechanical.w0[0][0, 0]/c.n
            if not np.isclose(covariance[0, 0], expected_initial,
                    rtol=4*np.finfo(float).eps, atol=0.):
                raise ValueError('original canonical initial covariance/normalization differs')
            with (output/'covariance.npz').open('xb') as stream:
                np.savez(stream, covariance=covariance, times=np.array([0., 30.]))
            with (output/'pullbacks.npy').open('xb') as stream:
                np.save(stream, np.stack(duals), allow_pickle=False)
            record = dict(schema=SCHEMA, status='complete-at-stated-scope',
                kind='density-K-coupled-Gaussian-at-archived-common-mean-nodes',
                input_mean=mean.binding, input_identity=mean.input_identity,
                numerical_reference=mean.reference_identity, coefficient_manifest=coefficient_binding,
                cache=mean.spec['force_cache'], compact_source=mean.spec['compact_source'],
                probe=dict(name='canonicalGammaNi1/x', index=0, normalization='1sqrtN',
                    coefficient=float(1/c.sqrt_n), q=0), n=str(int(c.n)), times=[0., 30.],
                dimensions=dict(slots=c.ns, bands=c.m, ionic=c.d, field=c.f,
                    joint_positions=c.d+c.f, full_phase=model.nx,
                    coupled_state=model.state_size(0), propagated_probe_columns=1),
                raw_archive=binding(output/'covariance.npz'), pullbacks=binding(output/'pullbacks.npy'),
                normalized_ordered=normalized_ordered(covariance[1, 0], covariance[1, 1], covariance[0, 0]),
                initial_preparation=dict(mean_added=False, original_rho0_W0_equal=True),
                settings=dict(method='DOP853', rtol=2e-8, atol=2e-10, precision='complex128'),
                approximation='Prescribed literal A78mean/rho nodes with A57linear node interpolation, including mean-dependent Gamma curvature and original coupled density/K feedback. No U/S or kernel interpolation, no mean reintegration. Node-matched, not the unarchived continuous A78DOP853mean; interpolation/temporal error unresolved. Not intrinsic-only loss or a full114 saved control history.',
                new_block_oracle=case['new_block_oracle'], source_bindings=case['source_bindings'],
                preregistration=case['preregistration'], elapsed_seconds=time.monotonic()-started,
                nonclaims=['Fullintrinsictrajectory', 'Continuousmeanmatch', 'Temporalaccuracy',
                    'Resolvedsmallmeandifference', 'Fixedpoint', 'Physicalcalibration/rate/heat'])
            save_json(output/'control.json', record)
    consumed = consume(output/'control.json')
    save_json(output/'saved-consumer.json', consumed)
    return dict(control=binding(output/'control.json'), consumer=binding(output/'saved-consumer.json'),
        elapsed_seconds=time.monotonic()-started, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    args = parser.parse_args()
    case = json.loads(args.case.read_bytes())
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    if (sorted(os.sched_getaffinity(0)) != [case['cpu']] or sys.executable != contract['python'] or
            any(os.environ.get(k) != v for k, v in contract['threads'].items())):
        raise ValueError('actual owned CPU/pinned runtime/thread contract required')
    cgroup = Path('/proc/self/cgroup').read_text()
    unified = [line.split('::', 1)[1] for line in cgroup.splitlines() if line.startswith('0::')]
    if len(unified) != 1:
        raise ValueError('actual declared unified user-unit resource context required')
    root = Path('/sys/fs/cgroup')/unified[0].lstrip('/')
    memory, swap = (root/'memory.max').read_text().strip(), (root/'memory.swap.max').read_text().strip()
    if memory != '3221225472' or swap != '0':
        raise ValueError('actual declared3GiB/noSwap resource limits required')
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)), executable=sys.executable,
        threads={k: os.environ.get(k) for k in contract['threads']}, cgroup=cgroup,
        memory_max=memory, swap_max=swap, pythonpath=os.environ.get('PYTHONPATH'))
    receipts = Path(case['receipts'])
    atomic_json(receipts/'runtime.json', runtime)
    for item in [case['mean'], case['preregistration'], case['new_block_oracle'], *case['source_bindings']]:
        authenticate(item)
    print('A83_ACTUAL_OWNED_CONTROL_BEGIN', flush=True)
    result = produce(case)
    atomic_json(receipts/'result.json', result)
    print('A83_ACTUAL_CONTROL_COMPLETE '+json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
