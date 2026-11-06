"""Complete one-map dressed EE/EI/common-Gamma/cubic functional evaluation."""
from dressed_reference import *
import time, resource
from all_slot_coulomb_shared import all_time_q
from all_slot_coulomb_contractions import zeros
from run_all_slot_coulomb import BareCoulombFactors,load_coulomb_inputs
from electron_ion_functional import EdgePair,evaluate_components
from material_electron_ion_inputs import edge_pair
from joint_gamma_adapter import vertices,common_pair as bare_common_pair
from cubic_contour import CubicContour,ordered_matrix
SOURCE=Path('/var/tmp/p254-a47-full37-source-001/source.json')

def common_pair(c,ref,a,b,ionic_only=False):
    forward,reverse=vertices(c);dg,dl=ref.ordered(0,a,b)
    if ionic_only:forward,reverse=forward[:,:c.d],reverse[:,:c.d];dg,dl=dg[:c.d,:c.d],dl[:c.d,:c.d]
    gg,gl=ref.electron(a,b);bg,bl=ref.electron(b,a);total={}
    for source in range(c.ns):
        target=c.target[0,source];left=forward[source];right=reverse[source].transpose(0,2,1)
        pair=EdgePair(left,right,left,right,gg[source],gl[source],gg[target],gl[target],bg[source],bl[source],dg,dl)
        for key,value in evaluate_components(pair,n=1,t=ref.times[a],s=ref.times[b]).items():
            if key.startswith('sigma'):
                if key not in total:total[key]=np.zeros((c.ns,c.m,c.m),complex)
                total[key][target]+=value
            else:total[key]=total.get(key,0)+value
    return total

def oracle(dataset,catalog,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    c,ref,bindings=inputs(dataset,catalog);bare=ref.bare;nt=len(ref.times)
    er=np.zeros_like(ref.er);el=np.empty_like(er);dr=np.zeros_like(ref.dr);dl=np.empty_like(dr)
    genuine={}
    for a,b in PAIRS:
        eg,ell=bare.electron(a,b);dg,dll=bare.ordered(0,a,b)
        el[a,b]=ell;dl[a,b]=dll;genuine[a,b]=(eg,dg)
        if a>=b:er[a,b]=eg-ell;dr[a,b]=dg-dll
    errors={'electron_greater_maxabs':max(float(np.max(np.abs(greater(er,el,a,b)-genuine[a,b][0]))) for a,b in PAIRS),'gamma_greater_maxabs':max(float(np.max(np.abs(greater(dr,dl,a,b)-genuine[a,b][1]))) for a,b in PAIRS)}
    # Actual callback with independently obtained genuine bare input, not rho/U reconstruction.
    dressed=(ref.er,ref.el,ref.dr,ref.dl);ref.er,ref.el,ref.dr,ref.dl=er,el,dr,dl
    actual=common_pair(c,ref,1,0);expected=bare_common_pair(c,bare,1,0)
    errors['common_callback_maxabs']=max(float(np.max(np.abs(actual[k]-expected[k]))) for k in actual)
    mutations=dict(wrong_equal_time_R_minus_advanced=float(np.linalg.norm((er[1,1]-er[1,1].swapaxes(-1,-2).conj())-er[1,1])),wrong_future_reverse_sign=float(np.linalg.norm((el[0,1]+er[1,0].swapaxes(-1,-2).conj())-genuine[0,1][0])),wrong_boson_lesser_adjoint=float(np.linalg.norm(dl[0,1].conj().T+spectral(dr,0,1)-genuine[0,1][1])))
    provider=CubicContour.from_source(SOURCE);indices=provider.source.indices
    routed=[provider.source.qindex(indices[0]-indices[k]) for k in range(64)]
    changed=[k for k,l in enumerate(routed) if k==0 or l==0]
    if changed!=[0] or routed[0]!=0:raise AssertionError('q0 replacement is not disjoint')
    ref.er,ref.el,ref.dr,ref.dl=dressed
    record=dict(schema='p254-a56-new-callback-boundary-oracle-v1',status='passed',input_dataset=ref.dataset,baseline_catalog=ref.catalog_binding,genuine_bare_bindings=bare.bindings,errors=errors,mutations=mutations,theta_at_zero=1,raw_diagonal_preserved=True,cubic_changed_internal_q=changed,cubic_routed_partner_q=routed,cubic_identity='sum_all64_new = paid_sum_all64_old - old_internal_q0 + new_internal_q0; for external Gamma l=-k, only k=0 uses changed Gamma D',no_old_functional_proofs_replayed=True)
    passed=max(errors.values())<=2e-10 and min(mutations.values())>=1e-9
    record['status']='passed' if passed else 'failed'
    save(output/'oracle.json',record);ref.close()
    if not passed:raise AssertionError(record)
    return record

def ee(dataset,catalog,output,q_start,q_stop,wall_cap):
    started=time.monotonic();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    c=load_coulomb_inputs(density_cap_bytes=900*1024**2)
    # Coulomb coefficient object has only electronic fields; use the same saved full G without a mechanical loader.
    with np.load(dataset,allow_pickle=False) as z:times=z['times'];er=z['electron_ret'];el=z['electron_less']
    electrons={(a,b):(greater(er,el,a,b),el[a,b]) for a,b in PAIRS}
    f=BareCoulombFactors(c);arrays=zeros(c,len(times));timings=[];completed=[]
    checkpoint=output/'ee.npz'
    if checkpoint.exists():
        with np.load(checkpoint,allow_pickle=False) as z:
            arrays={k:z[k] for k in arrays};completed=z['completed_outer_q'].tolist()
    for q in range(q_start,q_stop):
        if q in completed:continue
        began=time.monotonic();delta=all_time_q(f,q,electrons,times,external_batch=2,time_batch=4)
        for order in ('greater','lesser'):
            for channel in ('direct','exchange'):arrays[order+'_'+channel]+=delta[order+'_'+channel]
            arrays[order+'_total'][:]=arrays[order+'_direct']+arrays[order+'_exchange']
        theta=(times[:,None]>=times[None,:])[None,:,:,None,None]
        for channel in ('direct','exchange','total'):arrays['retarded_'+channel][:]=theta*(arrays['greater_'+channel]-arrays['lesser_'+channel])
        completed.append(q);seconds=time.monotonic()-began;timings.append(dict(q=q,seconds=seconds));elapsed=time.monotonic()-started
        np.savez(checkpoint,**arrays,completed_outer_q=completed,times_inverse_ry=times)
        receipt=dict(status='complete' if len(completed)==q_stop-q_start else 'partial',completed_outer_q=completed,owned_outer_q=list(range(q_start,q_stop)),inner_h=64,slots=128,bands=12,pairs=9,timings=timings,elapsed_seconds=elapsed,forecast_seconds=elapsed+(q_stop-q-len([x for x in completed if x>q])-1)*seconds,archive=binding(checkpoint),input_dataset=binding(dataset),coefficient_bindings=c.bindings)
        save(output/'scope.json',receipt);print(json.dumps(dict(stage='EE-q',q=q,seconds=seconds,elapsed=elapsed,forecast=receipt['forecast_seconds'])),flush=True)
        if elapsed+1.5*seconds>wall_cap-10 and q+1<q_stop:raise RuntimeError('preserved q checkpoint before cap')
    return receipt

def mixed(dataset,catalog,output):
    started=time.monotonic();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    c,ref,bindings=inputs(dataset,catalog);shape=(3,3,c.ns,c.m,c.m)
    sigma={o:np.zeros(shape,complex) for o in ORDERS};pi={o:np.zeros((3,3,114,114),complex) for o in ORDERS};electrons={p:ref.electron(*p) for p in PAIRS};timings=[]
    from substrate_framework.partitioned_q_storage import release_mapped_pages
    for q in range(64):
        began=time.monotonic();data=c.q(q)
        for a,b in PAIRS:
            boson=ref.boson(q,a,b)
            for source in range(c.ns):
                target=int(c.target[q,source]);pair=edge_pair(data,source,target,electrons[a,b],electrons[b,a],boson)
                values=evaluate_components(pair,n=c.n,t=ref.times[a],s=ref.times[b],components=('sigma',))
                for o in ORDERS:sigma[o][a,b,target]+=values['sigma_'+o]
        timings.append(dict(q=q,seconds=time.monotonic()-began))
        np.savez(output/'ei.npz',**sigma,completed_internal_q=np.arange(q+1))
        save(output/'ei-scope.json',dict(completed_internal_q=list(range(q+1)),timings=timings,elapsed_seconds=time.monotonic()-started))
        print(json.dumps(dict(stage='EI-q',q=q,seconds=timings[-1]['seconds'],elapsed=time.monotonic()-started)),flush=True)
        for name in ('ionic','ionic_reverse','density','potential','kernel'):release_mapped_pages(getattr(data,name))
    joint={};ionic={}
    for a,b in PAIRS:
        full=common_pair(c,ref,a,b);ion=common_pair(c,ref,a,b,True)
        for key,value in full.items():
            if key not in joint:joint[key]=np.empty((3,3,*value.shape),complex)
            joint[key][a,b]=value
        for key,value in ion.items():
            if key not in ionic:ionic[key]=np.empty((3,3,*value.shape),complex)
            ionic[key][a,b]=value
        for o in ORDERS:sigma[o][a,b]+=full['sigma_'+o]-ion['sigma_'+o];pi[o][a,b]=full['pi_'+o]
    np.savez(output/'joint.npz',**{'joint_'+k:v for k,v in joint.items()},**{'ionic_'+k:v for k,v in ionic.items()})
    np.savez(output/'mixed.npz',**{'sigma_'+o:sigma[o] for o in ORDERS},**{'pi_'+o:pi[o] for o in ORDERS})
    result=dict(status='complete',archive=binding(output/'mixed.npz'),joint_archive=binding(output/'joint.npz'),ei_archive=binding(output/'ei.npz'),elapsed_seconds=time.monotonic()-started,timings=timings,coefficient_bindings=bindings,reference_bindings=ref.bindings,combination='raw all64 EI Sigma + joint Gamma Sigma - ionic Gamma Sigma; common Pi not stacked')
    save(output/'result.json',result);ref.close();return result

def cubic(dataset,catalog,output):
    started=time.monotonic();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    c,ref,bindings=inputs(dataset,catalog);cat=ref.catalog;provider=CubicContour.from_source(SOURCE)
    l,left,right=provider.pair(0,0)
    if l!=0:raise ValueError('Gamma q0 routing differs')
    authenticate(cat['archives']['cubic'])
    with np.load(cat['archives']['cubic']['path'],allow_pickle=False) as z:
        if not np.array_equal(z['completed_internal_q'],np.arange(64)):raise ValueError('paid cubic sum incomplete')
        arrays={o:z[o].copy() for o in ORDERS}
    old0={o:np.empty((3,3,111,111),complex) for o in ORDERS};new0={o:np.empty_like(old0[o]) for o in ORDERS}
    # Independently paid, genuine five-pair q0 checkpoint fixes orientation.
    paid_q0=HERE.parent/'0055/material/cubic/Gamma-Pi.npz'
    paid_q0_binding=binding(paid_q0);orientation_errors=[]
    with np.load(paid_q0,allow_pickle=False) as paid:
        if not np.array_equal(paid['completed_internal_q'],[0]):raise ValueError('orientation anchor is not pure q0')
        paid_pairs=set(map(tuple,paid['new_time_pairs'].tolist()))
        orientation_anchor={o:paid[o] for o in ORDERS}
    for a,b in PAIRS:
        old=ref.bare.boson(0,a,b);new=ref.boson(0,a,b)
        for i,o in enumerate(('greater','lesser')):
            old0[o][a,b]=ordered_matrix(left,right,old[i],old[i]);new0[o][a,b]=ordered_matrix(left,right,new[i],new[i])
        for bank in (old0,new0):bank['retarded'][a,b]=(bank['greater'][a,b]-bank['lesser'][a,b]) if a>=b else 0
        if (a,b) in paid_pairs:
            error=max(float(np.max(np.abs(old0[o][a,b]-orientation_anchor[o][a,b]))) for o in ORDERS)
            orientation_errors.append(dict(pair=[a,b],maxabs=error))
            if error>2e-10:raise ValueError('paid q0 original orientation identity failed')
        for o in ORDERS:arrays[o][a,b]+=new0[o][a,b]-old0[o][a,b]
        np.savez(output/'Gamma-Pi.npz',**arrays,updated_pairs=np.array([p for p in PAIRS if p<=(a,b)]),completed_internal_q=np.arange(64),times_inverse_ry=ref.times)
        print(json.dumps(dict(stage='cubic-q0-pair',a=a,b=b,elapsed=time.monotonic()-started)),flush=True)
    np.savez(output/'q0-replacement.npz',**{'old_'+o:old0[o] for o in ORDERS},**{'new_'+o:new0[o] for o in ORDERS})
    result=dict(status='complete',archive=binding(output/'Gamma-Pi.npz'),replacement=binding(output/'q0-replacement.npz'),elapsed_seconds=time.monotonic()-started,paid_complete_sum=cat['archives']['cubic'],source=binding(SOURCE),source_arrays=provider.source.source_record['arrays'],scope='exact full64 paid sum minus original q0 plus dressed q0; original otherq raw opposite-q bath unchanged',routed_partner_q=0)
    result['paid_original_q0_orientation_anchor']=paid_q0_binding
    result['paid_original_q0_orientation_errors']=orientation_errors
    save(output/'result.json',result);ref.close();return result
