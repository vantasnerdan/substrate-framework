"""Reusable complete dressed-functional map; no represented-control recomputation."""
from dressed_reference import *
import subprocess, time, concurrent.futures

def consume(catalog_path):
    cat=json.loads(Path(catalog_path).read_bytes());path=authenticate(cat['archive']);result={}
    with np.load(path,allow_pickle=False) as z:
        for kind,shape in (('sigma',(3,3,128,12,12)),('pi',(3,3,114,114))):
            for o in ORDERS:
                value=z[kind+'_'+o]
                if value.shape!=shape or not np.isfinite(value).all():raise ValueError('complete finite '+kind+' required')
                vector=np.linspace(.1,.9,shape[-1])*(1+.37j);action=value@vector
                result[kind+'_'+o]=dict(shape=list(value.shape),norm=float(np.linalg.norm(value)),action_norm=float(np.linalg.norm(action)),action_sum=[float(action.sum().real),float(action.sum().imag)],sample=[float(value[1,0].flat[0].real),float(value[1,0].flat[0].imag)])
        times=np.array(cat['times_inverse_ry']);theta=(times[:,None]>=times[None,:])
        result['raw_contour_defects']={kind:float(np.max(np.abs(z[kind+'_retarded']-theta.reshape(3,3,*([1]*(z[kind+'_retarded'].ndim-2)))*(z[kind+'_greater']-z[kind+'_lesser'])))) for kind in ('sigma','pi')}
    return result

def assemble(dataset,baseline_catalog,output):
    output=Path(output);baseline=json.loads(Path(baseline_catalog).read_bytes());times=baseline['selected_times_inverse_ry']
    ee={};coverage=[]
    for name in ('ee-a','ee-b','ee-c'):
        scope=json.loads((output/name/'scope.json').read_bytes())
        if scope['status']!='complete':raise ValueError('partial EE bank')
        coverage+=scope['completed_outer_q']
        with np.load(authenticate(scope['archive']),allow_pickle=False) as z:
            for o in ORDERS:
                for channel in ('direct','exchange','total'):
                    key=o+'_'+channel;ee[key]=ee.get(key,0)+z[key].transpose(1,2,0,3,4)
    if sorted(coverage)!=list(range(64)):raise ValueError('all64 disjoint q coverage required')
    np.savez(output/'ee-complete.npz',**ee,completed_outer_q=np.arange(64),times_inverse_ry=times)
    with np.load(output/'mixed/mixed.npz',allow_pickle=False) as mixed,np.load(output/'cubic/Gamma-Pi.npz',allow_pickle=False) as cubic:
        arrays={}
        for o in ORDERS:
            arrays['sigma_'+o]=ee[o+'_total']+mixed['sigma_'+o]
            arrays['pi_'+o]=mixed['pi_'+o].copy();arrays['pi_'+o][:,:,:111,:111]+=cubic[o]
    # Bare mechanical source keeps contacts, not an already represented bubble.
    # Common Pi is intrinsic; represented-only Pi remains a separate control.
    represented=baseline['archives']['represented'];authenticate(represented)
    archive=output/'kernels.npz';np.savez(archive,**arrays)
    # Actual old/new numeric changes use original A53 consumer, not invented baseline combinations.
    from material_step_inputs import load_step
    c,reference,oldcat,old_sigma,old_pi,control=load_step(baseline_catalog)
    changes={};deltas={}
    for kind,old in (('sigma',old_sigma),('pi',old_pi)):
        for o in ORDERS:
            delta=arrays[kind+'_'+o]-old[o]
            deltas[kind+'_'+o]=delta
            changes[kind+'_'+o]=dict(old_norm=float(np.linalg.norm(old[o])),new_norm=float(np.linalg.norm(arrays[kind+'_'+o])),difference_norm=float(np.linalg.norm(delta)),maxabs=float(np.max(np.abs(delta))))
    reference.close()
    np.savez(output/'kernel-changes.npz',**deltas)
    costs=[]
    for p in sorted((HERE/'invocations').glob('*/exit.json')):
        receipt=json.loads(p.read_bytes());native=p.parent/'native-exit.json'
        if native.exists():receipt['native']=json.loads(native.read_bytes())
        costs.append(receipt)
    elapsed=sum(r.get('native',{}).get('elapsed_seconds',r['elapsed_seconds']) for r in costs)
    ee_seconds=sum(r.get('native',{}).get('elapsed_seconds',r['elapsed_seconds']) for r in costs if r['mode']=='ee')
    cubic_seconds=sum(r.get('native',{}).get('elapsed_seconds',r['elapsed_seconds']) for r in costs if r['mode']=='cubic')
    source_bindings=[binding(p) for p in HERE.glob('*.py')]
    cat=dict(schema='p254-dressed-functional-components-v1',status='complete',
        times_inverse_ry=times,input_dataset=binding(dataset),
        baseline_catalog=binding(baseline_catalog),archive=binding(archive),
        sigma={o:dict(array='sigma_'+o,shape=[3,3,128,12,12]) for o in ORDERS},
        pi={o:dict(array='pi_'+o,shape=[3,3,114,114]) for o in ORDERS},
        source_bindings=source_bindings,represented_control=represented,
        represented_control_unchanged=True,actual_old_new_kernel_changes=changes,
        reference_source_kind='bare-mechanical-with-contacts',
        pi_convention='full-common-joint-plus-cubic; represented-only-control-separate',
        exact_physics_scope=dict(
            electron='all128 all12 bands full two-time dressed G, all64 outer and internal routings and complete features',
            boson='full114 dressed Gamma D; original explicit raw opposite-q otherq bath',
            sigma='all64 EI + common joint Gamma - ionic Gamma + complete EE direct/exchange',
            pi='full common114 Gamma plus full111 cubic once; represented-only control separate',
            unchanged=['original preparation','scalar mean','active K','contacts','passive complement'],
            map_iterations=1,selfconsistent=False,conservation_established=False,accuracy_established=False),
        costs=dict(invocations=costs,new_native_seconds=elapsed,new_native_cap_seconds=5000,
            ee_append_seconds=ee_seconds,ee_cumulative_seconds=8733.680+ee_seconds,ee_ceiling_seconds=24000,
            cubic_append_seconds=cubic_seconds,cubic_cumulative_seconds=9126.572+cubic_seconds,cubic_ceiling_seconds=9500),
        findings=[])
    cat['complete_ee_direct_exchange_total']=binding(output/'ee-complete.npz')
    cat['actual_old_new_kernel_change_arrays']=binding(output/'kernel-changes.npz')
    cat['new_callback_boundary_oracle']=binding(output/'boundary/oracle.json')
    cat['preregistered_source_inventory']=binding(HERE/'preregistration.json')
    cat['execution_amendments']=[binding(p) for p in sorted(HERE.glob('execution-amendment-*.json'))]
    cat['component_bindings']={name:binding(output/name/'result.json') for name in ('mixed','cubic')}
    save(output/'catalog.json',cat);consumer=consume(output/'catalog.json');save(output/'saved-consumer.json',consumer)
    cat['actual_saved_numeric_consumer']=binding(output/'saved-consumer.json');cat['findings']=[dict(status='observed',claim='complete one-map dressed functional evaluation',evidence=cat['archive']),dict(status='observed',claim='authenticated saved numeric consumer applies all six kernels',evidence=cat['actual_saved_numeric_consumer']),dict(status='scope',claim='one map is not selfconsistency/conservation/accuracy; nonGamma bath unchanged'),dict(status='observed',claim='raw contour defects retained without state repairs',values=consumer['raw_contour_defects'])]
    save(output/'catalog.json',cat);return cat

def evaluate_dressed(dataset,baseline_catalog,output,*,cpus=(6,7,8),native_cap_seconds=5000,run_boundary_oracle=True,boundary_oracle=None):
    """Execute complete reusable map via bounded single-thread systemd units.

    dataset is the full A53 saved two-time NPZ; baseline_catalog is the paid
    A53 physical component catalog. output is a fresh owned directory. The
    original paid represented-control bank remains unchanged and separate.
    """
    if tuple(cpus)!=(6,7,8):raise ValueError('this allocation owns CPU6/7/8 only')
    if not 4830<=native_cap_seconds<=5000:raise ValueError('complete bounded grant sum4830s must fit new native ceiling5000s')
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    def launch(mode,name,cpu,cap,*extra):
        name=name+'-'+hashlib.sha256(str(output).encode()).hexdigest()[:8]
        argv=[sys.executable,str(HERE/'launch.py'),mode,name,str(cpu),str(cap),str(dataset),str(baseline_catalog),str(output),*map(str,extra)]
        code=subprocess.call(argv)
        if code:raise RuntimeError('native component failed: '+name)
    if run_boundary_oracle:launch('oracle','boundary-001',6,120)
    else:
        if boundary_oracle is None:raise ValueError('explicit authenticated carried boundary oracle required when not rerunning it')
        entry=boundary_oracle if isinstance(boundary_oracle,dict) else binding(boundary_oracle)
        carried=json.loads(authenticate(entry).read_bytes())
        if carried['status']!='passed':raise ValueError('cannot carry a failed boundary')
        save(output/'boundary/oracle.json',dict(schema='p254-a56-carried-boundary-v1',status='passed',original_oracle=entry,scope='unchanged supplied callback boundary; no old proof replay'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        work=[pool.submit(launch,'ee','ee-'+chr(97+i)+'-001',cpus[i],1300,lo,hi) for i,(lo,hi) in enumerate(((0,22),(22,43),(43,64)))]
        for f in work:f.result()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        work=[pool.submit(launch,'mixed','mixed-001',6,450),pool.submit(launch,'cubic','cubic-001',8,240)]
        for f in work:f.result()
    launch('assemble','assemble-001',6,120)
    from finalize_handoff import run as finalize
    finalize(output)
    return json.loads((output/'catalog.json').read_bytes())
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--dataset',type=Path,default=DEFAULT_DATASET);p.add_argument('--catalog',type=Path,default=DEFAULT_CATALOG);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();evaluate_dressed(a.dataset,a.catalog,a.output)
