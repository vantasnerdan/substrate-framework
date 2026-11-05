"""Consume all actual producer arrays, all128 slots/fullGamma and original cache."""
from pathlib import Path
import json
import numpy as np
from produce_fullhistory_mean import HERE, OUT, SOURCE, CompactCubic
from history_storage import authenticate, binding, atomic_json, StoredNPZ
from saved_mean import SavedMean, POSITIONS
from fullwindow_force_cache import FullWindowForceCache


def consume(result_path):
    spec = json.loads(Path(result_path).read_bytes())
    snapshots = json.loads((HERE/'historical-source-snapshots.json').read_bytes())['records']
    with StoredNPZ(authenticate(spec['archive'])) as producer, SavedMean(result_path) as saved:
        if (saved.binding == saved.reference_identity or saved.reference_identity != spec['input_identity']['numerical_reference'] or
                saved.input_identity != spec['input_identity'] or saved.is_bare):
            raise AssertionError('actual artifact and numerical reference conflated')
        nt = len(saved.times)
        electron_errors = []
        for slot in range(128):
            ret,less = saved.electron_slot(slot)
            maximum = 0.
            for a in range(nt):
                for b in range(nt):
                    end = producer['U'][b,slot].conj().T
                    expected_less = 1j*(producer['U'][a,slot]@producer['rho0'][slot]@end)
                    expected_ret = (-1j*producer['U'][a,slot]@end) if a >= b else np.zeros((12,12),complex)
                    maximum = max(maximum,float(np.max(np.abs(less[a,b]-expected_less))),
                        float(np.max(np.abs(ret[a,b]-expected_ret))))
            electron_errors.append(maximum)
        gr,gl = saved.gamma_slot()
        gamma_error = 0.
        for a in range(nt):
            for b in range(nt):
                forward = producer['S_gamma'][a]@producer['W0_gamma']@producer['S_gamma'][b].T
                reverse = producer['S_gamma'][b]@producer['W0_gamma']@producer['S_gamma'][a].T
                expected_less = -1j*reverse[np.ix_(POSITIONS,POSITIONS)].T
                expected_ret = -1j*forward[np.ix_(POSITIONS,POSITIONS)]-expected_less if a >= b else np.zeros((114,114),complex)
                gamma_error = max(gamma_error,float(np.max(np.abs(gl[a,b]-expected_less))),
                    float(np.max(np.abs(gr[a,b]-expected_ret))))
        source = CompactCubic.from_source(SOURCE)
        with FullWindowForceCache(authenticate(spec['force_cache']),spec['coefficient_manifest'],
                historical_source_snapshots=snapshots) as cache:
            actual = cache.forcing_for(authenticate(spec['input_identity']['actual_dataset']),source)
            force_checks = {key:np.array_equal(actual[key],producer[member]) for key,member in
                (('rho','node_rho'),('ordered_gamma_covariance','node_ordered_gamma_covariance'),
                 ('cubic_covariance_force','node_cubic_covariance_force'))}
            restoration = cache.source_authentication
        findings = dict(result=saved.binding,archive=spec['archive'],input_identity=saved.input_identity,
            reader=binding(HERE.parent/'0063/saved_mean.py'),
            force_reader=binding(HERE.parent/'0074/fullwindow_force_cache.py'),
            consumer_implementation=binding(Path(__file__)),
            historical_source_snapshots=snapshots,cache_source_authentication=restoration,
            electron_slots=list(range(128)),electron_ret_lesser_max_errors=electron_errors,
            fullGamma_ret_lesser_max_error=gamma_error,Gamma_shape=list(gl.shape),
            saved_mean_byte_identical=np.array_equal(saved.mean,producer['mean_gamma']),
            saved_U_byte_identical=np.array_equal(saved.u,producer['U']),
            saved_S_byte_identical=np.array_equal(saved.s,producer['S_gamma']),
            original_cache_reopen_equal=force_checks,
            nonGamma_contractions=actual['nonGamma_contractions'],Gamma_contractions=actual['Gamma_contractions'],
            times=saved.times.tolist(),initial_preparation=spec['defects'],
            final_electron_lesser_norm=float(np.linalg.norm(less[-1,-1])),
            final_fullGamma_lesser_norm=float(np.linalg.norm(gl[-1,-1])),
            final_mean_norm=float(np.linalg.norm(saved.mean[-1])),
            MeanFeedbackReference_compatibility='Not compatible: legacy reader requires its own historical schema and genuine baseline_catalog. Not executed and no catalog fabricated.',
            scope='Actual persisted producer arrays consumed; no independent ODE/physics oracle or reintegration.')
        if (max(electron_errors)>1e-11 or gamma_error != 0 or not all(force_checks.values()) or
                not all(findings[k] for k in ('saved_mean_byte_identical','saved_U_byte_identical','saved_S_byte_identical'))):
            raise AssertionError('actual saved consumer differs from actual producer')
        atomic_json(OUT/'saved-consumer.json',findings)
        return dict(consumer=binding(OUT/'saved-consumer.json'),findings=findings)
