"""Owned, preregistered A60 boundary and material pilot; never full all-q replay."""
from pathlib import Path
import os,sys,json,time,hashlib,resource,traceback
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE)]+[str(HERE.parent/n) for n in ('0056','0046','0044','0041')]
import numpy as np
from streamed_ee import stream_ee, NumericSink, read_blocks, _operands
from all_slot_coulomb_shared import all_time_q, static_superoperators, exchange_superoperator
from all_slot_coulomb_blas import all_time_q as independent_time_q
from all_slot_coulomb_contractions import KEYS
from bare_coulomb_factors import BareCoulombFactors
from coulomb_material_inputs import load_coulomb_inputs
from dressed_reference import greater
from material_electron_ion_inputs import FrozenReference

START=time.monotonic()
OUT=HERE/'pilot-001'
OUT.mkdir(exist_ok=False)
REPORT=dict(schema='p254-a60-pilot-v1',findings=[],timings=[],failures=[],prior_ee_seconds=15722.9054025717266)

def binding(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())

def save(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def measure(label,fn):
    t=time.monotonic();value=fn();elapsed=time.monotonic()-t
    REPORT['timings'].append(dict(label=label,seconds=elapsed,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
    print(json.dumps(REPORT['timings'][-1]),flush=True)
    save('progress.json',REPORT)
    return value

def comparison(label,new,old):
    findings=[]
    for k in KEYS:
        d=new[k]-old[k]
        entry=dict(label=label,channel=k,maxabs=float(np.max(np.abs(d))),difference_norm=float(np.linalg.norm(d)),reference_norm=float(np.linalg.norm(old[k])),relative=float(np.linalg.norm(d)/max(np.linalg.norm(old[k]),1e-300)))
        findings.append(entry)
        if not np.allclose(new[k],old[k],rtol=3e-12,atol=3e-12):raise AssertionError(entry)
    REPORT['findings'].extend(findings)
    return findings

def collect(factors,electron,times,**kwargs):
    # Deliberately small finite oracle only. Production uses NumericSink.
    c=factors.c;nt=len(times)
    arrays={k:np.zeros((c.ns,nt,nt,c.m,c.m),complex) for k in KEYS}
    coverage=np.zeros((c.nq,c.ns,nt,nt),np.int16)
    def consume(block):
        for j,(a,b) in enumerate(block.pairs):
            coverage[block.q,block.external,a,b]+=1
            for k in KEYS:arrays[k][block.external,a,b]+=block.values[k][:,j]
    manifest=stream_ee(factors,electron,times,consume,**kwargs)
    return arrays,coverage,manifest

def boundary():
    rng=np.random.default_rng(605031)
    def complex_array(shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
    c=SimpleNamespace(ns=4,nq=2,n=2,m=2,minus=np.arange(2),
        target=np.asarray([[0,1,2,3],[1,0,3,2]]),slot_spin=np.array([0,0,1,1]),
        feature_partner=(np.arange(3),np.arange(3)),
        density=tuple(complex_array((4,3,2,2)) for _ in range(2)))
    factors=SimpleNamespace(c=c,weights=(np.array([.7,1.2,.3]),np.array([1.1,.8,.4])))
    # Asymmetric independent legs expose a,b and conjugation mistakes. No Green
    # symmetry is required by the exact contraction identity at this boundary.
    coeff=complex_array((2,4,2,2))
    def electron(a,b):
        return tuple(coeff[j]*(1+.013*a+.007j*b)+(.002*a*b+.003j*a*a)/(j+1) for j in range(2))
    times=np.arange(31,dtype=float)[::-1] # not index-ordered retarded masks
    old=independent_time_q(factors,1,{(a,b):electron(a,b) for a in range(31) for b in range(31)},times,external_batch=3)
    # Original low-level q evaluator returns primitives, not assembled channels.
    for order in ('greater','lesser'):old[order+'_total']=old[order+'_direct']+old[order+'_exchange']
    mask=times[:,None]>=times[None,:]
    for ch in ('direct','exchange','total'):old['retarded_'+ch]=(old['greater_'+ch]-old['lesser_'+ch])*mask[None,:,:,None,None]
    new,coverage,manifest=collect(factors,electron,times,owned_q=(1,),pair_batch=7,external_batch=3,static_slots=3)
    assert np.all(coverage[1]==1) and np.all(coverage[0]==0)
    comparison('finite31-all961-tail-independent-local-band-chain',new,old)
    # Mutation magnitudes are numerical, not boolean source-text checks.
    swapped=np.linalg.norm(new['greater_exchange']-new['greater_exchange'].swapaxes(1,2))
    def wrong(a,b):return tuple(x.conj() for x in electron(a,b))
    bad,_,_=collect(factors,wrong,times,owned_q=(1,),owned_pairs=((30,29),),pair_batch=2,external_batch=3,static_slots=3)
    conjugate=np.linalg.norm(new['greater_exchange'][:,30,29]-bad['greater_exchange'][:,30,29])
    tail=np.linalg.norm(new['greater_total'][:,30,30])
    assert min(swapped,conjugate,tail)>1e-5
    REPORT['findings'].append(dict(label='finite31-sensitive-mutations',swapped_time_axes_norm=float(swapped),conjugated_provider_norm=float(conjugate),dropped_final_pair_norm=float(tail),coverage_exact=True,manifest=manifest))
    save('finite31-boundary.json',dict(manifest=manifest,findings=REPORT['findings']))

def assembled(old,times):
    for order in ('greater','lesser'):old[order+'_total']=old[order+'_direct']+old[order+'_exchange']
    mask=np.asarray(times)[:,None]>=np.asarray(times)[None,:]
    for ch in ('direct','exchange','total'):old['retarded_'+ch]=(old['greater_'+ch]-old['lesser_'+ch])*mask[None,:,:,None,None]
    return old

def persisted_run(label,factors,electron,times,bindings,pair_batch):
    sink=NumericSink(OUT/(label+'.npyblocks'))
    manifest=stream_ee(factors,electron,times,sink,owned_q=(1,),pair_batch=pair_batch,external_batch=2,static_slots=16,source_bindings=bindings)
    record=sink.close(manifest)
    return record

def read_consumer(path,nt,reference=None,selected=None):
    # Real persisted consumer reduces norms/checks reference blockwise, never
    # allocates nt² kernel. Coverage is tiny metadata, not physics payload.
    coverage=np.zeros((128,nt,nt),np.int16)
    norms={k:0. for k in KEYS};diff={k:0. for k in KEYS};maxabs={k:0. for k in KEYS}
    for block in read_blocks(path):
        assert block.q==1
        for j,(a,b) in enumerate(block.pairs):
            coverage[block.external,a,b]+=1
            for k in KEYS:
                value=block.values[k][:,j];norms[k]+=float(np.vdot(value,value).real)
                if reference is not None and (selected is None or (a,b) in selected):
                    aa,bb=(a,b) if selected is None else selected[(a,b)]
                    d=value-reference[k][block.external,aa,bb]
                    diff[k]+=float(np.vdot(d,d).real);maxabs[k]=max(maxabs[k],float(np.max(np.abs(d))))
                    if not np.allclose(value,reference[k][block.external,aa,bb],rtol=3e-12,atol=3e-12):raise AssertionError((k,a,b,maxabs[k]))
    assert np.all(coverage==1)
    return dict(complete_owned_q_coverage=True,norms={k:float(np.sqrt(v)) for k,v in norms.items()},difference_norm={k:float(np.sqrt(v)) for k,v in diff.items()},maxabs=maxabs)

def material():
    c=measure('load-original-all128-all64-density',load_coulomb_inputs)
    factors=BareCoulombFactors(c)
    dataset=Path('/var/tmp/p254-a58-common-feedback-step-002/actual-material-step.npz')
    db=binding(dataset)
    assert db['sha256']=='e98de95ce506e9db54514c72a5dcd6ef410aed64d7b5e4a5310d1488b3b3a33b'
    with np.load(dataset,allow_pickle=False) as z:times=z['times'];er=z['electron_ret'];el=z['electron_less']
    def electron(a,b):return greater(er,el,a,b),el[a,b]
    bindings=[c.bindings['manifest'],db]
    REPORT['material_inputs']=dict(density_bindings=c.bindings,dataset=db,physical_dimensions=dict(external=128,outer_global=64,owned_outer=[1],inner=64,bands=12,features=[len(w) for w in factors.weights],spins=2))
    # Benchmark real useful layout alternatives at unchanged physical size.
    def layouts():
        matrix,b,cc=static_superoperators(factors,1,np.arange(2))
        x,y,_=_operands(factors,1,electron,((0,1),(2,0)))
        xb=x[:,b].transpose(1,0,2,3,4);yb=y[:,cc].transpose(1,0,2,3,4)
        t=time.monotonic();batched=exchange_superoperator(matrix,xb,yb);tb=time.monotonic()-t
        t=time.monotonic();single=np.concatenate([exchange_superoperator(matrix,xb[:,i:i+1],yb[:,i:i+1]) for i in range(4)],axis=1);ts=time.monotonic()-t
        defect=float(np.max(np.abs(single-batched)))
        assert np.allclose(single,batched,rtol=3e-12,atol=3e-12)
        result=dict(label='actual-layout-batched4-versus-single1',batched_seconds=tb,single_seconds=ts,maxabs=defect,selected_pair_batch=2 if tb<=ts else 1,external=2,inner=64,orders=4)
        REPORT['findings'].append(result)
        return result['selected_pair_batch']
    pair_batch=measure('actual-layout-alternatives',layouts)
    old=measure('actual-dressed-original-A46-q1-all3times',lambda:assembled(all_time_q(factors,1,{(a,b):electron(a,b) for a in range(3) for b in range(3)},times,external_batch=2,time_batch=2),times))
    record=measure('actual-dressed-streamed-q1-all3times',lambda:persisted_run('dressed3',factors,electron,times,bindings,pair_batch))
    consumer=measure('actual-dressed-persisted-consumer',lambda:read_consumer(OUT/'dressed3.npyblocks',3,old))
    REPORT['findings'].append(dict(label='actual-NiH-dressed-old-new-all9channels',record=record,consumer=consumer))
    del old,er,el
    # Genuine five saved nodes, not U/rho reconstruction of dressed Green legs.
    bundle=Path('/var/tmp/p254-a52-requested-reference-001/producer/bundle.json')
    spec=json.loads(bundle.read_bytes());ub=binding(spec['U']['path'])
    assert ub==spec['U']
    assert binding(bundle)['sha256']=='c38f087b57e818517cb268d30c3be53a9305ead7eee107ebe85dd569c0fa99f0'
    u31=np.load(spec['U']['path'],mmap_mode='r',allow_pickle=False)
    ix=[0,1,15,29,30]
    reference=SimpleNamespace(c=c,u=u31[ix])
    five=np.asarray(spec['times_inverse_ry'])[ix]
    def free(a,b):return FrozenReference.electron(reference,a,b)
    ubindings=[c.bindings['manifest'],binding(bundle),ub]
    record5=measure('actual-U31-genuine5-streamed-q1-all25pairs',lambda:persisted_run('genuine5',factors,free,five,ubindings,pair_batch))
    tailtimes=five[-2:]
    oldtail=measure('actual-U31-independent-old-A46-tail2times',lambda:assembled(all_time_q(factors,1,{(a,b):free(a+3,b+3) for a in range(2) for b in range(2)},tailtimes,external_batch=2,time_batch=2),tailtimes))
    selected={(a,b):(a-3,b-3) for a in range(3,5) for b in range(3,5)}
    consumer5=measure('actual-genuine5-persisted-consumer',lambda:read_consumer(OUT/'genuine5.npyblocks',5,oldtail,selected))
    REPORT['findings'].append(dict(label='actual-NiH-genuine5-U31-all25pairs-nonzero-tail',times=five.tolist(),U_indices=ix,record=record5,consumer=consumer5,scope='bare moving-reference actual electron Green; not dressed31/full31 collisions'))

sources=[HERE/'streamed_ee.py',HERE/'run_pilot.py']+[HERE.parent/n for n in ('0046/all_slot_coulomb_shared.py','0046/all_slot_coulomb_blas.py','0046/all_slot_coulomb_contractions.py','0044/bare_coulomb_factors.py','0044/coulomb_material_inputs.py','0041/material_electron_ion_inputs.py','0056/dressed_reference.py','0052/general_moving_reference.py')]
REPORT['sources']=[binding(p) for p in sources]
save('preregistered-contract.json',dict(sources=REPORT['sources'],claim='Exact bounded arbitrary-time EE implementation; no full31 physical closure claimed',boundary='independent old BLAS local-band chain at all961 finite31 pairs, odd7 tails, external3 tails, unsorted times; numeric orientation/tail mutations',actual='q1 full128/all64inner/all12bands physical features/spins/weights, actual complex dressed3 all9 pairs old/new and persisted numeric consumer; genuine5 nodes0,1,15,29,30 from actualU31 all25pairs plus original-A46 final2node comparison',limits=dict(cpu=8,threads=1,memory=3221225472,swap=0,native_seconds=600,new_bytes=67108864),accepted='A46 static factors/math/physical pair identity admitted; no Fock proof replay',units='unchanged A46 physical EE units and weights',tolerance='maxabs/relative allclose rtol3e-12 atol3e-12, all channel findings retained',plan='one boundary then actual feature/layout profile; source inventory before dependent work'))
try:
    measure('finite31-narrow-boundary-once',boundary)
    material()
    REPORT['status']='complete'
except BaseException as e:
    REPORT['status']='failed';REPORT['failures'].append(dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc()))
    raise
finally:
    REPORT['native_elapsed_seconds']=time.monotonic()-START
    REPORT['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    REPORT['ee_ledger_upper_seconds']=REPORT['prior_ee_seconds']+REPORT['native_elapsed_seconds']
    REPORT['new_payload_bytes']=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
    save('report.json',REPORT)
    print('A60_NUMERIC_TERMINAL '+REPORT['status'],flush=True)
