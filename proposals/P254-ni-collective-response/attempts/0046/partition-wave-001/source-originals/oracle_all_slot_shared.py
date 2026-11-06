"""One owning comparison for distinct shared-Green144x144 regrouping only."""
import json
from pathlib import Path
from types import SimpleNamespace
from time import perf_counter
import numpy as np
from all_slot_coulomb_checkpoint import AllSlotPrefix,binding,atomic_json,publish_result
from all_slot_coulomb_shared import evaluate_time_grid,exchange_superoperator,bounds
from run_all_slot_coulomb import BareCoulombFactors
from second_born_functional import contract
HERE = Path(__file__).resolve().parent


def run():
    start = perf_counter()
    out = HERE/'shared-oracle-001'
    out.mkdir(exist_ok=False)
    fixture = HERE/'finite-oracle-output'/'actual-independent-inputs.npz'
    expected_path = HERE/'finite-oracle-output'/'all-slots-all-pairs.npz'
    with np.load(fixture,allow_pickle=False) as a:
        density = tuple(a[f'density_{q}'] for q in range(4))
        electrons = {(i,j):tuple(a[f'leg_{i}_{j}_{o}'] for o in range(2)) for i in range(4) for j in range(4)}
    gs = (np.array([[0,0,0],[1,0,0],[-1,0,0]]),np.array([[0,0,0],[1,0,0]]),
        np.array([[0,0,0],[0,0,1],[1,0,0],[-1,0,1]]),np.array([[0,0,0],[-1,0,0]]))
    partners = (np.array([0,2,1]),np.array([0,1]),np.array([1,0,3,2]),np.array([0,1]))
    mesh = np.array([[0,0,k] for k in range(4)])
    c = SimpleNamespace(n=4,ns=8,nq=4,m=2,density=density,slot_spin=np.repeat([1,2],4),
        slot_k=np.tile(np.arange(4),2),k_indices=mesh,q_indices=mesh,
        target=np.array([[s//4*4+(s%4+q)%4 for s in range(8)] for q in range(4)]),
        minus=np.array([0,3,2,1]),reciprocal=np.eye(3),tpiba=1.,feature_g=gs,feature_partner=partners)
    spec = dict(nq=4,shape=[8,4,4,2,2],external_batches=[[0,1,2],[3,4,5],[6,7]],
        authenticated_sources_and_inputs=[binding(fixture),binding(expected_path),
            *[binding(HERE/n) for n in ('all_slot_coulomb_shared.py','oracle_all_slot_shared.py',
            'all_slot_coulomb_blas.py','all_slot_coulomb_checkpoint.py','all_slot_coulomb_contractions.py')]])
    factors = BareCoulombFactors(c); times = [0.,.1,.4,.9]
    def interrupt(q):
        if q==1:
            raise InterruptedError('shared-Green committedq1 interruption')
    try:
        evaluate_time_grid(factors,electrons,times,AllSlotPrefix(out/'prefix',spec,spec['shape'],4),
            external_batch=3,time_batch=7,after_commit=interrupt)
        raise AssertionError('shared partial checkpoint path not exercised')
    except InterruptedError:
        pass
    actual,counts = evaluate_time_grid(factors,electrons,times,
        AllSlotPrefix(out/'prefix',spec,spec['shape'],4),external_batch=3,time_batch=7)
    with np.load(expected_path,allow_pickle=False) as a:
        errors = {k:float(np.linalg.norm(actual[k]-a[k])/max(1e-30,np.linalg.norm(a[k]))) for k in a.files}
    assert max(errors.values())<2e-11,errors
    saved = publish_result(out,spec,actual,counts,[])
    assert all(np.array_equal(saved['arrays'][k],actual[k]) for k in actual)
    completed,reuse = evaluate_time_grid(factors,electrons,times,
        AllSlotPrefix(out/'prefix',spec,spec['shape'],4),external_batch=3,time_batch=7)
    assert reuse['new_q']==0 and all(np.array_equal(completed[k],actual[k]) for k in actual)
    rng = np.random.default_rng(4627)
    rnd = lambda s:(rng.normal(size=s)+1j*rng.normal(size=s))*.05
    x,y,z,w = rnd((46,12,12)),rnd((46,12,12)),rnd((27,12,12)),rnd((27,12,12))
    fb,fc,rev = rnd((12,12)),rnd((12,12)),rnd((12,12))
    raw = (z.transpose(1,2,0).reshape(144,27)@w.reshape(27,144)).reshape(12,12,12,12)
    matrix = raw.transpose(0,3,1,2).reshape(1,144,144)
    xx = (x@fb)[None,None]
    yy = (rev[None]@y@fc[None])[None,None]
    expected = -contract(('rib','rac','sdf','sej','bd','ce','fa'),(x,y,z,w,fb,fc,rev),'ij')
    full = exchange_superoperator(matrix,xx,yy)[0,0]
    error12 = float(np.linalg.norm(full-expected)/np.linalg.norm(expected))
    wrong = exchange_superoperator(matrix,xx,(rev.conj()[None]@y@fc[None])[None,None])[0,0]
    mutation = float(np.linalg.norm(wrong-expected)/np.linalg.norm(expected))
    assert error12<2e-11 and mutation>1e-4
    # Width12,46feature,256localquartets,8timeorders: actual source routine.
    # Independent random perturbations avoid identical-operand benchmark tricks.
    mm = np.broadcast_to(matrix,(256,144,144)).copy()
    xb = np.broadcast_to(xx,(256,8,46,12,12)).copy()
    yb = np.broadcast_to(yy,(256,8,46,12,12)).copy()
    xb[:,1:] *= np.linspace(.9,1.1,7)[None,:,None,None,None]
    t = perf_counter()
    large = exchange_superoperator(mm,xb,yb)
    batch_seconds = perf_counter()-t
    assert np.allclose(large[:,0],np.broadcast_to(full,large[:,0].shape),rtol=2e-11,atol=2e-11)
    projection = batch_seconds*(32/8)*(128/4)*64
    report = dict(schema='p254-a46-shared-Green-superoperator-oracle-v1',status='passed',
        seconds=perf_counter()-start,all_slot16pair_relative_errors=errors,
        full12band_46x27feature_relative_error=error12,conjugating_reverse_mutation=mutation,
        representative_kernel=dict(local_quartets=256,time_orders=8,features=46,bands=12,
            seconds=batch_seconds,projected_full_exchange_only_seconds=projection,
            forecast_scope='synthetic full-width source-kernel timing; NOT actual material throughput; excludes static/shared/direct/checkpoint/input costs'),
        partial_resume=counts,completed_resume=reuse,source_input_specification=spec,
        saved_result=saved['record']['archive'],bounds=bounds(),
        findings=[dict(claim='new shared(q,c,timepair)superoperator matches complete original finite allslot/spin/nonzeroq/Nyquist/unequalfeature/16pair outputs',status='passed'),
                  dict(claim='actual12band unequal46x27feature independent complex/bilinear local equality and reversed-leg sensitivity',status='passed')],blockers=[])
    atomic_json(out/'report.json',report)
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    run()
