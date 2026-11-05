"""One bounded A51 comparison: Gamma, internal q1 and its momentum partner."""
from pathlib import Path
import hashlib
import json
import os
import resource
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in ('0051','0052','0041')]
from cubic_contour import CubicContour, ordered_matrix
from material_electron_ion_inputs import load_material
from general_moving_reference import SavedMovingReference
from gaussian_cubic_time import GaussianCubicTimeSeparated, stationarity_defect


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1048576), b''):
            digest.update(data)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest.hexdigest())


def save(path, value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')


def run(output):
    output = Path(output)
    started = time.monotonic()
    tick = time.monotonic()
    c, manifest, coefficient_binding = load_material()
    source_path = '/var/tmp/p254-a47-full37-source-001/source.json'
    provider = CubicContour.from_source(source_path)
    reference = SavedMovingReference(c,'/var/tmp/p254-a52-requested-reference-001/producer/bundle.json',coefficient_binding,times=[0.,30.])
    load_seconds = time.monotonic()-tick
    tick = time.monotonic()
    partner,left,right = provider.pair(0,1)
    vertex_seconds = time.monotonic()-tick
    evaluator = GaussianCubicTimeSeparated(reference,chunk=16)
    qs = sorted({1,partner,int(c.minus[1]),int(c.minus[partner])})
    raw, stationarity = {}, {}
    for q in qs:
        defect = stationarity_defect(reference,q)
        a = c.mechanical.generators[q]
        w = c.mechanical.w0[q]
        am = c.mechanical.generators[int(c.minus[q])]
        raw['generator_'+str(q)] = a
        raw['w0_'+str(q)] = w
        raw['stationarity_defect_'+str(q)] = defect
        for index,t in enumerate(reference.times):
            raw['S_'+str(q)+'_'+str(index)] = reference.S(q,t)
        scale = np.linalg.norm(a@w)+np.linalg.norm(w@am.T)
        stationarity[str(q)] = dict(raw_defect_frobenius=float(np.linalg.norm(defect)),
            raw_defect_max=float(np.max(np.abs(defect))),
            derivative_term_norm_sum=float(scale),relative_defect=float(np.linalg.norm(defect)/max(scale,np.finfo(float).tiny)),
            exact_entrywise_zero=bool(np.all(defect==0)),time_translation_used=False)
    raw.update(left=left,right=right,times=reference.times,external_q=np.array(0),internal_q=np.array(1),partner_q=np.array(partner))
    np.savez(output/'actual-inputs.npz',**raw)
    baseline = evaluator.zeros(111)
    tick = time.monotonic()
    leg_seconds = contraction_seconds = 0.
    for a in range(2):
        for b in range(2):
            lap = time.monotonic()
            dg1,dl1 = reference.boson(1,a,b)
            dg2,dl2 = reference.boson(partner,a,b)
            leg_seconds += time.monotonic()-lap
            for key,first,second in (('greater',dg1,dg2),('lesser',dl1,dl2)):
                lap = time.monotonic()
                baseline[key][a,b] = ordered_matrix(left,right,first,second)
                contraction_seconds += time.monotonic()-lap
    baseline_seconds = time.monotonic()-tick
    baseline['retarded'] = evaluator.retarded(baseline)
    candidate = evaluator.zeros(111)
    costs = evaluator.accumulate_q(left,right,1,partner,candidate)
    candidate['retarded'] = evaluator.retarded(candidate)
    archive = {'old_'+key:value for key,value in baseline.items()}
    archive.update({'new_'+key:value for key,value in candidate.items()})
    findings = []
    all_passed = True
    for key in ('greater','lesser','retarded'):
        delta = candidate[key]-baseline[key]
        scale = max(float(np.linalg.norm(baseline[key])),np.finfo(float).tiny)
        absolute = float(np.linalg.norm(delta))
        passed = absolute <= 2e-10*max(1.,scale)
        findings.append(dict(component=key,old_frobenius=scale,new_frobenius=float(np.linalg.norm(candidate[key])),
            difference_frobenius=absolute,difference_max=float(np.max(np.abs(delta))),relative_difference=absolute/scale,
            acceptance_threshold=2e-10*max(1.,scale),passed=passed))
        all_passed = all_passed and passed and np.isfinite(candidate[key]).all()
        archive['difference_'+key] = delta
    np.savez(output/'comparison-arrays.npz',**archive,times=reference.times,
        external_coordinates=np.arange(111),external_q=np.array(0),internal_q=np.array(1),partner_q=np.array(partner))
    with np.load(output/'comparison-arrays.npz',allow_pickle=False) as stored:
        persisted = stored['new_retarded'][1,0]
        probe = np.arange(1,112)/111 + 1j*np.arange(111)/113
        consumer = dict(saved_retarded_action_norm=float(np.linalg.norm(persisted@probe)),
                        saved_retarded_trace=[float(np.trace(persisted).real),float(np.trace(persisted).imag)])
    source = provider.source.source_record
    source_files = sorted({Path(module.__file__).resolve() for module in list(sys.modules.values())
                           if getattr(module,'__file__',None) and str(module.__file__).endswith('.py')
                           and ('substrate-ni-collective' in str(module.__file__) or 'substrate-framework/src' in str(module.__file__))})
    raw_archives = {entry['archive']['path']:entry['archive'] for entry in manifest['items'].values()}
    record = dict(schema='p254-a65-cubic-time-factor-material-pilot-v1',status='passed' if all_passed else 'failed',
        scope=dict(external_q=0,internal_q=1,partner_q=partner,ordered_internal_routes=[1],
            full_internal_mesh_evaluated=False,full_coordinates=111,source_electron_slots=128,bands=12,spins=2,
            times_inverse_ry=[0.,30.],time_pairs=4,bare_Gaussian_first_map_only=True,stationarity_assumed=False,
            field_context='original canonical generator retained; cubic Pi is ionic111 only',dressed_supported=False),
        findings=findings,stationarity=stationarity,
        timings=dict(load_seconds=load_seconds,vertex_seconds=vertex_seconds,baseline_seconds=baseline_seconds,
                     baseline_leg_seconds=leg_seconds,baseline_contraction_seconds=contraction_seconds,
                     candidate=costs,measured_contraction_speedup=contraction_seconds/costs['total_seconds'],
                     measured_with_legs_speedup=baseline_seconds/costs['total_seconds'],
                     material_wall_seconds=time.monotonic()-started),
        residency=dict(maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                       peak_feature_bank_bytes=costs['peak_feature_bank_bytes']),
        input_bindings=dict(material_manifest=coefficient_binding,material_archives=list(raw_archives.values()),
                            cubic_source=binding(source_path),cubic_arrays=source['arrays'],
                            cubic_source_bindings=source['source_bindings'],reference=reference.bindings,
                            raw_actual_factors_vertices=binding(output/'actual-inputs.npz')),
        implementation_bindings=[binding(p) for p in source_files],preregistration=binding(HERE/'preregistration.json'),
        outputs=dict(arrays=binding(output/'comparison-arrays.npz')),saved_consumer=consumer,
        nonclaims=['full31 or full64 production','dressed nonGaussian factorization','selfconsistency','conservation','rank truncation','stationary physical theorem','rigorous floating point enclosure'])
    save(output/'result.json',record)
    print(json.dumps(record),flush=True)
    reference.close()
    if not all_passed:
        raise AssertionError('actual A51 comparison exceeds frozen finite algebra tolerance')
    return record
