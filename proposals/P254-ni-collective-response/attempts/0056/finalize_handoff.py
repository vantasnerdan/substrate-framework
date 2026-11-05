"""Finalize exact execution costs/findings after all native terminal receipts."""
from dressed_reference import *

def run(output,closure=None):
    output=Path(output);cat=json.loads((output/'catalog.json').read_bytes());costs=[]
    for p in sorted((HERE/'invocations').glob('*/exit.json')):
        receipt=json.loads(p.read_bytes());receipt['launch']=binding(p.parent/'launch.json')
        receipt['manager_exit']=binding(p)
        for name in ('native-start','native-exit'):
            path=p.parent/(name+'.json')
            if path.exists():receipt[name.replace('-','_')]=dict(binding=binding(path),record=json.loads(path.read_bytes()))
        progress=[]
        for line in (p.parent/'stdout').read_text().splitlines():
            if line.startswith('{'):
                event=json.loads(line)
                if 'stage' in event:progress.append(event)
        receipt['actual_numeric_progress']=progress
        costs.append(receipt)
    # Manager elapsed conservatively includes imports/startup and every failure.
    def seconds(r):return r['elapsed_seconds']
    native=sum(seconds(r) for r in costs);ee=sum(seconds(r) for r in costs if r['mode']=='ee');cubic=sum(seconds(r) for r in costs if r['mode']=='cubic')
    ledger=dict(invocations=costs,new_native_seconds=native,new_native_cap_seconds=5000,ee_append_seconds=ee,ee_prior_seconds=8733.680,ee_cumulative_seconds=8733.680+ee,ee_ceiling_seconds=24000,cubic_append_seconds=cubic,cubic_prior_seconds=9126.572,cubic_cumulative_seconds=9126.572+cubic,cubic_new_ceiling_seconds=300,cubic_ceiling_seconds=9500,all_new_exitcodes=[r['exitcode'] for r in costs],failed_attempts=[r for r in costs if r['exitcode']!=0])
    if native>5000 or cubic>300 or 9126.572+cubic>9500:raise RuntimeError('actual cumulative runtime cap exceeded')
    cat['costs']=ledger
    cat['source_bindings']=[binding(p) for p in sorted(HERE.glob('*.py'))]
    cat['execution_source_revisions']=[binding(HERE/'preregistration.json'),*[binding(p) for p in sorted(HERE.glob('execution-amendment-*.json'))]]
    original=json.loads((HERE.parent/'0055/actual-genuine-midpoint-handoff-001.json').read_bytes())
    cat['physical_backend_source_bindings']=original['physical_backend_source_bindings']
    for entry in cat['physical_backend_source_bindings']:authenticate(entry)
    raw={}
    with np.load(authenticate(cat['input_dataset']),allow_pickle=False) as z:
        er=z['electron_ret'];el=z['electron_less'];dr=z['gamma_ret'];dl=z['gamma_less'];rho=-1j*np.array([el[a,a] for a in range(3)])
        raw['fermion_trace']=np.trace(rho,axis1=-2,axis2=-1).sum(axis=-1).real.tolist()
        raw['fermion_density_hermiticity_maxabs']=float(np.max(np.abs(rho-rho.conj().swapaxes(-1,-2))))
        raw['electron_lesser_two_time_antihermiticity_maxabs']=float(np.max(np.abs(el+el.swapaxes(0,1).conj().swapaxes(-1,-2))))
        raw['gamma_lesser_two_time_antihermiticity_maxabs']=float(np.max(np.abs(dl+dl.swapaxes(0,1).conj().swapaxes(-1,-2))))
        raw['electron_future_retarded_maxabs']=max(float(np.max(np.abs(er[a,b]))) for a,b in PAIRS if a<b)
        raw['gamma_future_retarded_maxabs']=max(float(np.max(np.abs(dr[a,b]))) for a,b in PAIRS if a<b)
        raw['electron_retarded_diagonal_norms']=[float(np.linalg.norm(er[a,a])) for a in range(3)]
        raw['gamma_retarded_diagonal_norms']=[float(np.linalg.norm(dr[a,a])) for a in range(3)]
    raw['repairs_applied']=[];cat['raw_dressed_input_defects']=raw
    cat['payload_bytes']=sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file())
    if cat['payload_bytes']>256*1024**2:raise RuntimeError('owned payload cap exceeded')
    cat['findings'] += [dict(status='resolved-interface',claim='A58 starts bare mechanical with contacts: full joint+cubic Pi intrinsic, represented Pi only separate control. No response/noise stacked or subtracted twice.'),dict(status='observed',claim='original paid represented control unchanged and separately authenticated',evidence=cat['represented_control']),dict(status='observed',claim='complete old/new numeric kernel differences retained',evidence=cat['actual_old_new_kernel_change_arrays']),dict(status='observed',claim='raw input defects without state repair',values=raw),dict(status='execution-failure-preserved',claim='startup cpuset path unavailable; observer source fixed with actual affinity enforcement',evidence=binding(HERE/'execution-amendment-002.json')),dict(status='oracle-failure-preserved',claim='bare Gamma retarded is real so transpose/adjoint mutation was insensitive; corrected exposing lesser-adjoint mutation, no physical formula/tolerance change',evidence=binding(HERE/'execution-amendment-003.json'))]
    cat['findings']=list({finding['claim']:finding for finding in cat['findings']}.values())
    if closure:cat['terminal_view_closure']=binding(closure)
    save(output/'catalog.json',cat)
    handoff=dict(schema='p254-a56-complete-dressed-map-handoff-v1',status='complete',catalog=binding(output/'catalog.json'),input_dataset=cat['input_dataset'],archive=cat['archive'],consumer=cat['actual_saved_numeric_consumer'],changes=cat['actual_old_new_kernel_change_arrays'],complete_ee=cat['complete_ee_direct_exchange_total'],api=dict(source=binding(HERE/'evaluate_dressed.py'),signature='evaluate_dressed(dataset,baseline_catalog,output,*,cpus=(6,7,8),native_cap_seconds=5000,run_boundary_oracle=True)'),costs=ledger,findings=cat['findings'],raw_dressed_input_defects=raw,source_bindings=cat['source_bindings'],physical_backend_source_bindings=cat['physical_backend_source_bindings'],payload_bytes=cat['payload_bytes'],scope=cat['exact_physics_scope'])
    handoff['api']['signature']='evaluate_dressed(dataset,baseline_catalog,output,*,cpus=(6,7,8),native_cap_seconds=5000,run_boundary_oracle=True,boundary_oracle=None)'
    handoff['execution_source_revisions']=cat['execution_source_revisions']
    if closure:handoff['terminal_view_closure']=binding(closure)
    save(output/'handoff.json',handoff);return handoff
if __name__=='__main__':
    run(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else None)
