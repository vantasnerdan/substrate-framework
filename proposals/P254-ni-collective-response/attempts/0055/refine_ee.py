"""A46 physical backend on only five genuine midpoint pairs and disjoint q."""
from refinement_common import *
import time,resource
from all_slot_coulomb_shared import static_superoperators,exchange_superoperator
from all_slot_coulomb_blas import direct_q
from run_all_slot_coulomb import BareCoulombFactors,load_coulomb_inputs
CHANNELS=('direct','exchange','total')
KEYS=tuple(o+'_'+c for o in ('greater','lesser','retarded') for c in CHANNELS)

def sparse_q(factors,q,electrons):
    c=factors.c;qm=int(c.minus[q]);p=c.feature_partner[q][:len(factors.weights[q])]
    aa=c.target[qm,np.arange(c.ns)];raw=c.density[qm][:,p].conj();lq=c.density[q].conj()
    shape=(c.ns,len(NEW_PAIRS),c.m,c.m)
    out={key:np.zeros(shape,complex) for key in KEYS}
    x=np.empty((10,c.ns,len(p),c.m,c.m),complex);y=np.empty_like(x);labels=[]
    for k,(a,b) in enumerate(NEW_PAIRS):
        f,r=electrons[a,b],electrons[b,a]
        for order,ff,rr in (('greater',f[0],r[1]),('lesser',f[1],r[0])):
            j=len(labels);x[j]=(lq@ff[:,None])*(factors.weights[q]/c.n)[None,:,None,None]
            y[j]=rr[aa,None]@raw@ff[:,None]
            labels.append((order,k));out[order+'_direct'][:,k]=direct_q(factors,q,ff,rr)
    for start in range(0,c.ns,2):
        i=np.arange(start,min(start+2,c.ns));matrix,b,cc=static_superoperators(factors,q,i)
        for begin in range(0,len(labels),4):
            end=min(begin+4,len(labels))
            xb=x[begin:end,b].transpose(1,0,2,3,4);yb=y[begin:end,cc].transpose(1,0,2,3,4)
            value=exchange_superoperator(matrix,xb,yb).reshape(len(i),c.nq,end-begin,c.m,c.m).sum(axis=1)
            for j,(order,k) in enumerate(labels[begin:end]):out[order+'_exchange'][i,k]=value[:,j]
    theta=np.array([TIMES[a]>=TIMES[b] for a,b in NEW_PAIRS])[None,:,None,None]
    for order in ('greater','lesser'):out[order+'_total']=out[order+'_direct']+out[order+'_exchange']
    for channel in CHANNELS:out['retarded_'+channel]=theta*(out['greater_'+channel]-out['lesser_'+channel])
    return out

def affected_oracle(output):
    """One finite mask/backend/oldmerge comparison; no old proof replay."""
    from types import SimpleNamespace
    from all_slot_coulomb_shared import all_time_q
    fixture=HERE.parent/'0046/finite-oracle-output/actual-independent-inputs.npz'
    with np.load(fixture,allow_pickle=False) as z:
        density=tuple(z[f'density_{q}'] for q in range(4))
        electrons={(a,b):tuple(z[f'leg_{a}_{b}_{o}'] for o in range(2)) for a in range(3) for b in range(3)}
    gs=(np.array([[0,0,0],[1,0,0],[-1,0,0]]),np.array([[0,0,0],[1,0,0]]),np.array([[0,0,0],[0,0,1],[1,0,0],[-1,0,1]]),np.array([[0,0,0],[-1,0,0]]))
    mesh=np.array([[0,0,k] for k in range(4)])
    c=SimpleNamespace(n=4,ns=8,nq=4,m=2,density=density,slot_spin=np.repeat([1,2],4),slot_k=np.tile(np.arange(4),2),k_indices=mesh,q_indices=mesh,target=np.array([[s//4*4+(s%4+q)%4 for s in range(8)] for q in range(4)]),minus=np.array([0,3,2,1]),reciprocal=np.eye(3),tpiba=1.,feature_g=gs,feature_partner=(np.array([0,2,1]),np.array([0,1]),np.array([1,0,3,2]),np.array([0,1])))
    f=BareCoulombFactors(c);new=sparse_q(f,0,electrons);full=all_time_q(f,0,electrons,TIMES,external_batch=2,time_batch=4)
    error=max(float(np.max(np.abs(new[key][:,k]-full[key][:,a,b]))) for key in ('greater_direct','lesser_direct','greater_exchange','lesser_exchange') for k,(a,b) in enumerate(NEW_PAIRS))
    old=np.arange(8*4*4*2*2).reshape(8,4,4,2,2).astype(complex)*(1+.37j)
    merged=old_merge(old,new['greater_total'],True)
    literal=all(np.array_equal(merged[:,a,b],old[:,oa,ob]) for a,oa in OLD_MAP.items() for b,ob in OLD_MAP.items())
    mutation=float(np.max(np.abs(new['greater_total'][:,0]-(full['greater_direct'][:,1,0]+full['greater_exchange'][:,1,0]))))
    if error>2e-10 or not literal or mutation<1e-9:raise AssertionError((error,literal,mutation))
    receipt=dict(status='passed',fixture=binding(fixture),new_five_pair_backend_maxabs=error,old_four_literal_equal=literal,swapped_orientation_mutation_maxabs=mutation,scope='one affected sparse mask/oldmerge/orientation finite comparison only')
    save(output/'affected-boundary-oracle.json',receipt);return receipt

def run(name,q_start,q_stop,cap):
    start=time.monotonic();out=HERE/'material'/name;out.mkdir(parents=True,exist_ok=False)
    if q_start==0:affected_oracle(out)
    c=load_coulomb_inputs(density_cap_bytes=900*1024**2);ref=reference(c,c.bindings['manifest']);f=BareCoulombFactors(c)
    electrons={(a,b):ref.electron(a,b) for a,b in set(NEW_PAIRS)}
    arrays={key:np.zeros((c.ns,5,c.m,c.m),complex) for key in KEYS};timings=[];owned=list(range(q_start,q_stop+1))
    save(out/'source-input.json',dict(reference=ref.bindings,inputs=c.bindings,sources=[binding(__file__),binding(HERE/'refinement_common.py')],owned_q=owned,new_pairs=NEW_PAIRS,scope='partial q contribution only on five midpoint pairs'))
    completed=[]
    for q in owned:
        began=time.monotonic();delta=sparse_q(f,q,electrons)
        for key in KEYS:arrays[key]+=delta[key]
        completed.append(q);seconds=time.monotonic()-began;timings.append(dict(q=q,seconds=seconds))
        path=out/'five-pairs.npz';np.savez(path,**arrays,new_pairs=np.array(NEW_PAIRS),completed_outer_q=completed,times_inverse_ry=TIMES)
        forecast=(64-len(completed))*seconds+time.monotonic()-start
        scope=dict(status='complete' if completed==owned else 'partial',contribution_scope='PARTIAL_Q_CONTRIBUTION',archive=binding(path),owned_q=owned,completed_outer_q=completed,external_slots=128,inner_h=64,bands=12,new_ordered_pairs=NEW_PAIRS,timings=timings,elapsed_seconds=time.monotonic()-start,remaining_global_serial_forecast_seconds=forecast,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        save(out/'scope.json',scope)
        print(json.dumps(dict(name=name,q=q,seconds=seconds,elapsed=scope['elapsed_seconds'],global_forecast=forecast)),flush=True)
        if len(completed)==1:
            decision=dict(first_q=timings[0],all64_serial_forecast_seconds=64*seconds,per_worker_forecast_seconds=len(owned)*seconds+time.monotonic()-start,EE_grant_seconds=1800,total_new_native_cap_seconds=6000,continue_within_cap=len(owned)*seconds+time.monotonic()-start<cap-12)
            save(out/'first-q-decision.json',decision)
            if not decision['continue_within_cap']:break
        if time.monotonic()-start+1.5*seconds>cap-12:break
    ref.close();return scope
