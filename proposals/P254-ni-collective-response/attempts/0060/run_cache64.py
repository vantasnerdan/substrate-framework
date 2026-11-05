"""Second measured layout: bounded64 static slots, actual saved numeric oracle."""
from pathlib import Path
import os,sys,json,time,hashlib,resource,traceback
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE)]+[str(HERE.parent/n) for n in ('0056','0046','0044','0041')]
import numpy as np
from streamed_ee import stream_ee,read_blocks
from all_slot_coulomb_contractions import KEYS
from bare_coulomb_factors import BareCoulombFactors
from coulomb_material_inputs import load_coulomb_inputs
from dressed_reference import greater
START=time.monotonic()
OUT=HERE/'cache64-002';OUT.mkdir(exist_ok=False)
report=dict(status='running',findings=[],failures=[],runtime=dict(affinity=sorted(os.sched_getaffinity(0)),env={k:os.environ.get(k) for k in ('HERDR_ENV','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},no_swap='systemd MemorySwapMax=0'),prior_pilot_exit=json.loads((HERE/'receipts-001/exit.json').read_bytes()))
def binding(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(1048576),b''):h.update(x)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
report['sources']=[binding(HERE/'streamed_ee.py'),binding(HERE/'run_cache64.py')]
save('preregistered-contract.json',dict(claim='Observed remaining memory permits64 static slots; two groups rather than eight, no static128 over3GiB',comparison='Actual dressed3 q1 all128/all64inner/all12bands/allfeatures exact complex saved numeric reference from paid first pilot, all9channels blockwise',budget_seconds=180,combined_native_ceiling=600,source_inventory=report['sources']))
try:
    t=time.monotonic();c=load_coulomb_inputs();report['load_seconds']=time.monotonic()-t
    factors=BareCoulombFactors(c)
    dataset=Path('/var/tmp/p254-a58-common-feedback-step-002/actual-material-step.npz')
    db=binding(dataset)
    assert db['sha256']=='e98de95ce506e9db54514c72a5dcd6ef410aed64d7b5e4a5310d1488b3b3a33b'
    with np.load(dataset,allow_pickle=False) as z:times=z['times'];er=z['electron_ret'];el=z['electron_less']
    def electron(a,b):return greater(er,el,a,b),el[a,b]
    # Original saved result is explicitly admitted as reference. The following
    # small3-time reference is not evaluator working state or31-history storage.
    reference={k:np.zeros((128,3,3,12,12),complex) for k in KEYS}
    for block in read_blocks(HERE/'pilot-001/dressed3.npyblocks'):
        for j,(a,b) in enumerate(block.pairs):
            for k in KEYS:reference[k][block.external,a,b]=block.values[k][:,j]
    coverage=np.zeros((128,3,3),np.int16)
    diff={k:0. for k in KEYS};maxabs={k:0. for k in KEYS};norm={k:0. for k in KEYS}
    def consume(block):
        for j,(a,b) in enumerate(block.pairs):
            coverage[block.external,a,b]+=1
            for k in KEYS:
                value=block.values[k][:,j];d=value-reference[k][block.external,a,b]
                diff[k]+=float(np.vdot(d,d).real);norm[k]+=float(np.vdot(value,value).real)
                maxabs[k]=max(maxabs[k],float(np.max(np.abs(d))))
                if not np.allclose(value,reference[k][block.external,a,b],rtol=3e-12,atol=3e-12):raise AssertionError((k,a,b,maxabs[k]))
    t=time.monotonic()
    manifest=stream_ee(factors,electron,times,consume,owned_q=(1,),pair_batch=2,external_batch=2,static_slots=64,source_bindings=(c.bindings['manifest'],db))
    report['evaluator_seconds']=time.monotonic()-t
    assert np.all(coverage==1)
    report['findings'].append(dict(label='actual-dressed3-static64-allchannels',coverage_exact=True,manifest=manifest,norm={k:float(np.sqrt(v)) for k,v in norm.items()},difference_norm={k:float(np.sqrt(v)) for k,v in diff.items()},maxabs=maxabs))
    report['status']='complete'
except BaseException as e:
    report['status']='failed';report['failures'].append(dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc()));raise
finally:
    report['elapsed_seconds']=time.monotonic()-START
    report['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    save('report.json',report)
    print('A60_CACHE64_NUMERIC_TERMINAL '+report['status'],flush=True)
