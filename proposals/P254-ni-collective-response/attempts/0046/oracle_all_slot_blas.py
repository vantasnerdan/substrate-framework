"""Only the new exact BLAS regrouping is exercised; old suites not replayed."""
import json
from pathlib import Path
from types import SimpleNamespace
from time import perf_counter
import numpy as np
from all_slot_coulomb_checkpoint import (AllSlotPrefix,binding,atomic_json,publish_result)
from all_slot_coulomb_blas import evaluate_time_grid,local_exchange,bounds
from run_all_slot_coulomb import BareCoulombFactors
from second_born_functional import contract
HERE = Path(__file__).resolve().parent


def run():
    start = perf_counter()
    out = HERE/'optimization-oracle-001'
    out.mkdir(exist_ok=False)
    fixture = HERE/'finite-oracle-output'/'actual-independent-inputs.npz'
    expected_path = HERE/'finite-oracle-output'/'all-slots-all-pairs.npz'
    with np.load(fixture,allow_pickle=False) as a:
        density = tuple(a[f'density_{q}'] for q in range(4))
        electrons = {(i,j):tuple(a[f'leg_{i}_{j}_{o}'] for o in range(2)) for i in range(4) for j in range(4)}
    gs = (np.array([[0,0,0],[1,0,0],[-1,0,0]]),np.array([[0,0,0],[1,0,0]]),
        np.array([[0,0,0],[0,0,1],[1,0,0],[-1,0,1]]),np.array([[0,0,0],[-1,0,0]]))
    p = (np.array([0,2,1]),np.array([0,1]),np.array([1,0,3,2]),np.array([0,1]))
    mesh = np.array([[0,0,k] for k in range(4)])
    c = SimpleNamespace(n=4,ns=8,nq=4,m=2,density=density,slot_spin=np.repeat([1,2],4),
        slot_k=np.tile(np.arange(4),2),k_indices=mesh,q_indices=mesh,
        target=np.array([[s//4*4+(s%4+q)%4 for s in range(8)] for q in range(4)]),
        minus=np.array([0,3,2,1]),reciprocal=np.eye(3),tpiba=1.,feature_g=gs,feature_partner=p)
    spec = dict(nq=4,shape=[8,4,4,2,2],external_batches=[[0,1,2],[3,4,5],[6,7]],
        authenticated_sources_and_inputs=[binding(fixture),binding(expected_path),
            *[binding(HERE/n) for n in ('all_slot_coulomb_blas.py','oracle_all_slot_blas.py',
            'all_slot_coulomb_checkpoint.py','all_slot_coulomb_contractions.py')]])
    factors = BareCoulombFactors(c); times = [0.,.1,.4,.9]
    def interrupted(q):
        if q == 1:
            raise InterruptedError('new BLAS prefix interruption afterq1')
    cache = AllSlotPrefix(out/'prefix',spec,spec['shape'],4)
    try:
        evaluate_time_grid(factors,electrons,times,cache,external_batch=3,after_commit=interrupted)
        raise AssertionError('new BLAS partial branch not exercised')
    except InterruptedError:
        pass
    actual,counts = evaluate_time_grid(factors,electrons,times,
        AllSlotPrefix(out/'prefix',spec,spec['shape'],4),external_batch=3)
    with np.load(expected_path,allow_pickle=False) as a:
        errors = {k:float(np.linalg.norm(actual[k]-a[k])/max(1e-30,np.linalg.norm(a[k]))) for k in a.files}
    assert max(errors.values())<2e-11,errors
    saved = publish_result(out,spec,actual,counts,[])
    assert all(np.array_equal(saved['arrays'][k],actual[k]) for k in actual)
    completed,reuse = evaluate_time_grid(factors,electrons,times,
        AllSlotPrefix(out/'prefix',spec,spec['shape'],4),external_batch=3)
    assert reuse['new_q']==0 and all(np.array_equal(completed[k],actual[k]) for k in actual)
    # Full12-band/unequal46x27-feature matrix kernel versus original7-operand
    # contraction. This exposes every local band ordering and true complex transpose.
    rng = np.random.default_rng(4612)
    rnd = lambda s:(rng.normal(size=s)+1j*rng.normal(size=s))*.05
    x,y,z,w = rnd((46,12,12)),rnd((46,12,12)),rnd((27,12,12)),rnd((27,12,12))
    fb,fc,rev = rnd((12,12)),rnd((12,12)),rnd((12,12))
    left = (x.transpose(1,2,0).reshape(144,46)@y.reshape(46,144)).reshape(1,12,12,12,12)
    right = (z.transpose(1,2,0).reshape(144,27)@w.reshape(27,144)).reshape(1,12,12,12,12)
    expected = -contract(('rib','rac','sdf','sej','bd','ce','fa'),(x,y,z,w,fb,fc,rev),'ij')
    t = perf_counter()
    full = local_exchange(left,right,fb[None],fc[None],rev[None])[0]
    local_seconds = perf_counter()-t
    local_error = float(np.linalg.norm(full-expected)/np.linalg.norm(expected))
    wrong = local_exchange(left,right,fb[None],fc[None],rev.conj()[None])[0]
    mutation = float(np.linalg.norm(wrong-expected)/np.linalg.norm(expected))
    assert local_error<2e-11 and mutation>1e-4
    # Representative actual-width batch256 through the exact new local kernel;
    # no material inputs consumed by this synthetic source performance probe.
    large_left = np.broadcast_to(left,(256,12,12,12,12)).copy()
    large_right = np.broadcast_to(right,(256,12,12,12,12)).copy()
    t = perf_counter()
    large = local_exchange(large_left,large_right,np.broadcast_to(fb,(256,12,12)),
        np.broadcast_to(fc,(256,12,12)),np.broadcast_to(rev,(256,12,12)))
    batch_seconds = perf_counter()-t
    assert np.allclose(large,np.broadcast_to(full,large.shape),rtol=2e-11,atol=2e-11)
    report = dict(schema='p254-a46-new-blas-regrouping-oracle-v1',status='passed',
        seconds=perf_counter()-start,finite_all_slot_relative_errors=errors,
        full12_band_46x27_feature_relative_error=local_error,conjugating_reverse_mutation=mutation,
        single_local_kernel_seconds=local_seconds,synthetic_batch256_local_seconds=batch_seconds,
        synthetic_only_material_speed_unmeasured=True,partial_resume=counts,completed_resume=reuse,
        source_input_specification=spec,result=saved['record']['archive'],bounds=bounds(),
        findings=[dict(claim='new BLAS complete all-slot/all16pair finite result matches original saved unreduced physical definition',status='passed'),
                  dict(claim='full12band unequal46x27-feature local regrouping preserves complex/bilinear convention',status='passed')],blockers=[])
    atomic_json(out/'report.json',report)
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    run()
