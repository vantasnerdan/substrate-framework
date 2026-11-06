"""One revised build, focused caller-limit boundary and original material smoke."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import numpy as np
from production_transport import (plan_layout, validate, FinalBank, ResultView,
    solve_supplied, load_solution, binding, array_digest, atomic_json, run_recorded)
from exercise_transport import free_control
HERE=Path(__file__).resolve().parent
ROOT=Path('/var/tmp/p254-a61-lean-002')
OLDROOT=Path('/var/tmp/p254-a61-lean-001')
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
PILOT=64*1024*1024
PRIOR=47.54140062816441


def rejection(label,fn,findings):
    try:fn()
    except (ValueError,RuntimeError) as exc:findings.append(dict(label=label,rejected=True,error=str(exc)))
    else:raise AssertionError('accepted invalid boundary: '+label)


def boundaries(binary):
    findings=[]
    plans={}
    for label,args in (('full31x114Gamma',(31,114,1)),('full31x128x12fermions',(31,12,128))):
        p=plan_layout(*args)
        assert p['bank_bytes']>PILOT
        assert p==plan_layout(*args,max_input_bytes=p['input_bytes'],max_result_bytes=p['result_bytes'],max_bank_bytes=p['bank_bytes'])
        plans[label]=p
        for name,size in (('max_input_bytes',p['input_bytes']),('max_result_bytes',p['result_bytes']),('max_bank_bytes',p['bank_bytes'])):
            rejection(label+'-'+name,lambda args=args,name=name,size=size:plan_layout(*args,**{name:size-1}),findings)
        findings.append(dict(label=label,accepted_without_pilot_cap=True,accepted_at_exact_caller_limits=True,allocated_physical_history=False))
    for value in (0,-1,True,1.5):
        rejection('invalid-caller-limit-'+str(value),lambda value=value:plan_layout(3,2,max_bank_bytes=value),findings)
    rejection('native-header-time-overflow',lambda:plan_layout(2**31+1,1),findings)
    rejection('native-matrix-overflow',lambda:plan_layout(3,50000),findings)
    times,arrays=free_control(5,3)
    small=plan_layout(5,3,2)
    rejection('bank-cap-before-creation',lambda:FinalBank(ROOT/'rejected-bank',times,3,2,statistics=-1,order=2,max_bank_bytes=small['bank_bytes']-1),findings)
    assert not (ROOT/'rejected-bank').exists()
    with ResultView(OLDROOT/'boundary-minimal.bin',max_result_bytes=plan_layout(5,3)['result_bytes']) as view:
        assert np.shares_memory(view.ret,view._mapped)
    rejection('mapped-result-caller-cap',lambda:ResultView(OLDROOT/'boundary-minimal.bin',max_result_bytes=plan_layout(5,3)['result_bytes']-1),findings)
    rejection('validate-caller-input-cap',lambda:validate(times,arrays,-1,2,1.,max_input_bytes=plan_layout(5,3)['input_bytes']-1),findings)
    # Distinct slots expose accidental slot-offset aliasing; zero selfenergy is
    # an explicitly free finite boundary fixture, not a material substitution.
    second=[np.asarray(a*.75,dtype='<c16') for a in arrays]
    with FinalBank(ROOT/'distinct-slots',times,3,2,statistics=-1,order=2,max_bank_bytes=small['bank_bytes']) as bank:
        for slot,source in enumerate((arrays,second)):
            solve_supplied(times,*source,statistics=-1,order=2,binary=binary,bank=bank,slot=slot,source_bindings=[binding(HERE/'exercise_transport.py'),binding(__file__)],max_input_bytes=None if slot==0 else PILOT,max_result_bytes=PILOT)
        assert np.array_equal(bank.ret[0],arrays[0]) and np.array_equal(bank.less[0],arrays[1])
        assert np.array_equal(bank.ret[1],second[0]) and np.array_equal(bank.less[1],second[1])
        bank.finalize()
    assert sum((ROOT/'distinct-slots'/(n+'.npy')).stat().st_size for n in ('ret','less'))==small['bank_bytes']
    with FinalBank(ROOT/'distinct-slots',resume=True,max_bank_bytes=small['bank_bytes']) as recovered:
        assert set(recovered.completed)=={0,1}
    findings.append(dict(label='distinct-slots-streamed-exact-layout-and-resume',passed=True,plan=small))
    # New native optional limit must reject before allocating the full history.
    rejection('native-explicit-byte-cap',lambda:run_recorded([str(binary),'/var/tmp/p254-a58-common-feedback-step-002/fermions/slot-000/input.bin',str(ROOT/'rejected-native.bin'),'1'],ROOT,'expected-native-byte-cap',10),findings)
    assert json.loads((ROOT/'expected-native-byte-cap.json').read_text())['exit']==2
    assert not (ROOT/'rejected-native.bin').exists()
    return dict(plans=plans,findings=findings,scope='exact planning only for31x114/31x128x12; no large physical arrays or banks created, no material31 solve')


def material(binary):
    old=Path('/var/tmp/p254-a58-common-feedback-step-002/fermions/slot-000')
    manifest=json.loads((old/'manifest.json').read_text())
    assert binding(old/'input.npz')==manifest['input']
    assert binding(old/'solution.npz')==manifest['solution']
    with np.load(old/'input.npz',allow_pickle=False) as data:
        times=data['times'];arrays=[data[n] for n in ('source_ret','source_less','kernel_ret','kernel_less')]
        order,beta=int(data['order']),float(data['beta'])
    plan=plan_layout(len(times),12)
    target=ROOT/'original-material-slot000'
    with FinalBank(target,times,12,1,statistics=-1,order=order,beta=beta,max_bank_bytes=plan['bank_bytes']) as bank:
        receipt=solve_supplied(times,*arrays,statistics=-1,order=order,beta=beta,binary=binary,bank=bank,slot=0,source_bindings=[binding(old/'manifest.json'),manifest['input']],max_input_bytes=plan['input_bytes'],max_result_bytes=plan['result_bytes'])
        bank.finalize()
    saved=load_solution(target)
    with np.load(old/'solution.npz',allow_pickle=False) as prior:
        oldr,oldl=array_digest(prior['ret']),array_digest(prior['less'])
        newr,newl=array_digest(saved['ret'][0]),array_digest(saved['less'][0])
        assert oldr==newr and oldl==newl
        np.savez(ROOT/'material-complete-findings.npz',ret_difference=saved['ret'][0]-prior['ret'],less_difference=saved['less'][0]-prior['less'])
    for key in ('ret','less'):saved[key]._mmap.close()
    return dict(old_source=manifest['input'],old_solution=manifest['solution'],old_ret=oldr,new_ret=newr,old_less=oldl,new_less=newl,canonical_bytes_equal=True,ret_max=0.,less_max=0.,exact_plan=plan,native=binding(target/'.staging'/'slot-000'/'native.json'),profile=binding(target/'.staging'/'slot-000'/'native.time'),findings=binding(ROOT/'material-complete-findings.npz'),completion=binding(target/'complete.json'),released=binding(target/'.staging'/'slot-000'/'released.json'))


def native():
    started=time.monotonic();code=1
    outcome=dict(affinity=sorted(os.sched_getaffinity(0)),herdr_env=os.environ.get('HERDR_ENV'),thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},sources=[binding(p) for p in sorted(HERE.iterdir()) if p.is_file()],original_snapshot=binding(OLDROOT/'pre-limit-revision-001.tar'),original_handoff=binding(HERE/'complete-handoff.json'))
    atomic_json(ROOT/'native-bindings.json',outcome)
    try:
        if outcome['affinity'] != [9] or outcome['herdr_env']!='1' or any(v!='1' for v in outcome['thread_environment'].values()):raise RuntimeError('resource grant mismatch')
        run_recorded(['cmake','-S',str(HERE),'-B',str(ROOT/'build'),'-DCMAKE_BUILD_TYPE=Release'],ROOT,'configure',30)
        run_recorded(['/usr/bin/time','-v','-o',str(ROOT/'build.time'),'cmake','--build',str(ROOT/'build'),'--parallel','1'],ROOT,'build',180)
        binary=ROOT/'build'/'solve_result'
        boundary=boundaries(binary)
        atomic_json(ROOT/'boundary-findings.json',boundary)
        actual=material(binary)
        record=dict(schema='p254-a61-supplemental-handoff-v1',reason='Remove erroneous fixed pilot allocation ceiling from production API/native IO',original_snapshot=outcome['original_snapshot'],original_handoff=outcome['original_handoff'],api=binding(HERE/'production-api.json'),contract=binding(HERE/'supplemental-contract.json'),source_bindings=binding(ROOT/'native-bindings.json'),new_binary=binding(binary),boundary=boundary,material=actual,all_exit_receipts=[binding(p) for p in sorted(ROOT.rglob('native.json'))]+[binding(ROOT/p) for p in ('configure.json','build.json','expected-native-byte-cap.json')],prior_native_manager_wall_seconds=PRIOR,physics='Unchanged direct A54/A49 headers; original full114Gamma3node and freeNt31 evidence carried via original_handoff, no material31 claim',resource_scope='CPU9/one-thread/3GiB/noSwap; full-history plans allocate no physical arrays; actual combined pilot payload remains below64MiB',layout_scope='Exact numeric/file layout support, not a large-history native working-set or performance qualification; existing full-slot A54 scratch remains')
        atomic_json(HERE/'supplemental-handoff.json',record)
        code=0
    except BaseException as exc:
        outcome['failure']=repr(exc);raise
    finally:
        own=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
        outcome.update(exit=code,wall_seconds=time.monotonic()-started,cpu_seconds=own.ru_utime+own.ru_stime+children.ru_utime+children.ru_stime,owner_maxrss_kib=own.ru_maxrss,child_maxrss_kib=children.ru_maxrss,payload_bytes=sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file()),combined_payload_bytes=sum(p.stat().st_size for root in (ROOT,OLDROOT) for p in root.rglob('*') if p.is_file()))
        atomic_json(ROOT/'native-outcome.json',outcome)
        print('[dft-done]' if code==0 else '[dft-fail]','A61 supplemental exit='+str(code),flush=True)


def launch():
    ROOT.mkdir(exist_ok=False)
    argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus','systemd-run','--user','--wait','--pipe','--collect','--unit=p254-a61-lean-002.service','--property=RuntimeMaxSec=240','--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=9','env','HERDR_ENV=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','BLIS_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1','taskset','-c','9',PYTHON,'-u',str(Path(__file__).resolve()),'--native']
    code=1
    try:run_recorded(argv,ROOT,'systemd-native',245);code=0
    finally:
        atomic_json(ROOT/'terminal.json',dict(exit=code,marker='A61_SUPPLEMENT_TERMINAL_EXIT_'+str(code)))
        print('A61_SUPPLEMENT_TERMINAL_EXIT_'+str(code),flush=True)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--native':native()
    else:launch()
