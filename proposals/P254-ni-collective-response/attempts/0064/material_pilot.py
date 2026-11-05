"""One actual q1 full-slot comparison at genuine saved times0,30."""
from pathlib import Path
import os,sys,json,time,resource,traceback,hashlib
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE)]+[str(HERE.parent/n) for n in ('0060','0052','0046','0044','0041','0045')]
import numpy as np
from gaussian_ee import stream_gaussian_ee
from streamed_ee import stream_ee,NumericSink,read_blocks
from all_slot_coulomb_contractions import KEYS
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors
from general_moving_reference import SavedMovingReference
OUT=HERE/'pilot-001';OUT.mkdir(exist_ok=False)
START=time.monotonic()
report=dict(schema='p254-a64-one-actual-comparison-v1',status='running',findings=[],failures=[],timings=[],scope='actual q1, all128 externals, all64 inner h, both spins, all12 bands, full original features; bare Gaussian first map only',tolerance=dict(rtol=3e-12,atol=3e-12))
def save():
    (OUT/'result.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
def binding(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def measure(label,fn):
    start=time.monotonic();value=fn()
    report['timings'].append(dict(label=label,seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
    print(json.dumps(report['timings'][-1]),flush=True);save();return value
try:
    report['runtime']=dict(affinity=sorted(os.sched_getaffinity(0)),environment={k:v for k,v in os.environ.items() if k.startswith('HERDR_') or k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    assert report['runtime']['affinity']==[6]
    cgroup=Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    report['runtime']['cgroup']={key:(Path('/sys/fs/cgroup')/cgroup.lstrip('/')/key).read_text().strip() for key in ('memory.max','memory.swap.max')}
    assert report['runtime']['cgroup']=={'memory.max':'3221225472','memory.swap.max':'0'}
    report['analytic']=json.loads((HERE/'symbolic-findings.json').read_bytes())
    report['frozen_preregistration']=binding(HERE/'preregistration.json')
    c=measure('load-authenticated-original-inputs',load_coulomb_inputs)
    factors=BareCoulombFactors(c)
    bundle=Path('/var/tmp/p254-a52-requested-reference-001/producer/bundle.json')
    reference=measure('open-genuine-saved-U0-U30',lambda:SavedMovingReference(c,bundle,c.bindings['manifest'],times=[0.,30.]))
    report['inputs']=dict(material=c.bindings,reference=reference.bindings,dimensions=dict(ns=c.ns,nq=c.nq,m=c.m,n=c.n,spins=2,features=[len(w) for w in factors.weights]))
    times=reference.times
    def run(label,fn):
        sink=NumericSink(OUT/(label+'.npyblocks'))
        manifest=fn(sink)
        return sink.close(manifest)
    new=measure('new-time-separated-q1-full-2times',lambda:run('new',lambda sink:stream_gaussian_ee(factors,reference,times,sink,owned_q=(1,),inner_batch=8,internal_batch=8,pair_batch=4,source_bindings=reference.bindings)))
    report['new']=new
    old=measure('old-A60-q1-full-2times',lambda:run('old',lambda sink:stream_ee(factors,reference.electron,times,sink,owned_q=(1,),static_slots=16,external_batch=2,pair_batch=2,source_bindings=reference.bindings)))
    report['old']=old
    def load(label):
        arrays={k:np.zeros((128,2,2,12,12),complex) for k in KEYS};coverage=np.zeros((128,2,2),np.int8)
        for block in read_blocks(OUT/(label+'.npyblocks')):
            assert block.q==1
            for j,(a,b) in enumerate(block.pairs):
                coverage[block.external,a,b]+=1
                for k in KEYS:arrays[k][block.external,a,b]+=block.values[k][:,j]
        assert np.all(coverage==1)
        np.savez(OUT/(label+'-all9.npz'),**arrays)
        return arrays
    def compare():
        a,b=load('new'),load('old')
        for k in KEYS:
            diff=a[k]-b[k];norm=float(np.linalg.norm(b[k]));dn=float(np.linalg.norm(diff))
            finding=dict(channel=k,maxabs=float(np.max(np.abs(diff))),difference_norm=dn,reference_norm=norm,relative=dn/max(norm,1e-300),allclose=bool(np.allclose(a[k],b[k],rtol=3e-12,atol=3e-12)))
            report['findings'].append(finding)
        assert all(v['allclose'] for v in report['findings'])
    measure('persisted-output-all9-consumer',compare)
    elapsed={v['label']:v['seconds'] for v in report['timings']}
    report['observed_speedup_old_over_new']=elapsed['old-A60-q1-full-2times']/elapsed['new-time-separated-q1-full-2times']
    report['long_history_verdict']='Temporal factors are reused over pairs but31-time speed is unmeasured; no full31 or dressed/selfconsistent physics claim.'
    report['outputs']=[binding(OUT/(label+'-all9.npz')) for label in ('new','old')]
    report['status']='complete'
except BaseException as exc:
    report['status']='failed';report['failures'].append(dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise
finally:
    report['native_elapsed_seconds']=time.monotonic()-START
    report['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    report['source_bindings']=[binding(m.__file__) for m in list(sys.modules.values()) if getattr(m,'__file__',None) and str(HERE.parent) in str(m.__file__) and str(m.__file__).endswith('.py')]
    report['source_bindings'] += [binding(HERE/'material_pilot.py'),binding(HERE/'symbolic_blocks.py'),binding(HERE/'symbolic-findings.json')]
    save();print('A64_NATIVE_TERMINAL '+report['status'],flush=True)
