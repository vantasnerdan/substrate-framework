"""Serialize unique finite component cases on leased Native cores; retain receipts."""
from pathlib import Path
import concurrent.futures
from datetime import datetime
import hashlib
import json
import os
import subprocess
import queue
import struct
import threading
import time

HERE = Path(__file__).resolve().parent
A80 = HERE.parent / '0080'
CONTRACT = json.loads((HERE/'preregistration.json').read_text())
PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'
CAP = CONTRACT['case_caps_seconds']
CEILING = {'ee':160000.0, 'mixed_ei':60000.0, 'cubic':40000.0}
OLD = {'ee':22950.79131010516390, 'mixed_ei':50.15633947495371, 'cubic':9811.458142627235497}
EXTRA_LEASE = HERE/'root-resource-extension-618.json'
EXTRA_DISPATCH_DEADLINE = datetime.fromisoformat(json.loads(EXTRA_LEASE.read_text())['new_dispatch_due']).timestamp()
EXTRA_CASE_QUOTA = int(json.loads(EXTRA_LEASE.read_text())['additional_native_cases_authorized'])
EXTRA_PANE = ('w1:p70','w1:t7H')
EXTRA_EI_LEASE = HERE/'root-resource-extension-621.json'
EXTRA_EI_DISPATCH_DEADLINE = datetime.fromisoformat(json.loads(EXTRA_EI_LEASE.read_text())['new_dispatch_due']).timestamp()
EXTRA_EI_CASE_QUOTA = int(json.loads(EXTRA_EI_LEASE.read_text())['additional_native_cases_authorized'])
EXTRA_EI_PANE = ('w1:p81','w1:t7F')
EXTRA_CUBIC_LEASE = HERE/'root-resource-extension-623.json'
EXTRA_CUBIC_DISPATCH_DEADLINE = datetime.fromisoformat(json.loads(EXTRA_CUBIC_LEASE.read_text())['new_dispatch_due']).timestamp()
EXTRA_CUBIC_CASE_QUOTA = int(json.loads(EXTRA_CUBIC_LEASE.read_text())['additional_native_cases_authorized'])
EXTRA_CUBIC_PANE = ('w1:p82','w1:t7F')
BANK_SLOTS = {6:('w1:p7N','w1:t77'),7:('w1:p7P','w1:t78'),8:('w1:p7Q','w1:t79'),9:('w1:p7R','w1:t7A'),10:('w1:p7T','w1:t1')}
# Exact three full31 numeric NPYs plus canonical slot bitmap, from paid catalogs.
OUTPUT_NPY_BYTES = {'ee':850354816,'mixed_ei':850354816,'cubic':568344561}
OUTPUT_LIMIT_BYTES = 8*1024**3

def ceilings():
    amendment=HERE/'ceiling-amendment.json'
    return dict(CEILING,**json.loads(amendment.read_text())['ceilings_seconds']) if amendment.exists() else CEILING

def save(path, value):
    path = Path(path)
    content=json.dumps(value,indent=2)+'\n'
    if path.exists() and path.read_text()==content:
        return
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(content)
    tmp.replace(path)

def binding(path):
    path=Path(path)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def output_reservations(cases):
    # Snapshot ownership before scanning files: later file creation consumes this
    # reservation, and concurrent bank retirement can only make it conservative.
    pending=[c for c in cases if c['operation']=='component' and
             not (Path(c['receipts'])/'exit.json').exists()]
    allocated=0
    by_output={}
    for path in Path('/var/tmp').glob('p254-a8[01]-*/*.npy'):
        try: size=path.stat().st_size
        except FileNotFoundError: continue  # Sole bank writer retired this array.
        allocated+=size
        by_output[path.parent]=by_output.get(path.parent,0)+size
    unallocated=sum(max(OUTPUT_NPY_BYTES[c['component']]-by_output.get(Path(c['output']),0),0)
                    for c in pending)
    return allocated,unallocated

def cli(*args):
    p=subprocess.run(['herdr',*args],capture_output=True,text=True)
    if p.returncode:
        raise RuntimeError(p.stderr or p.stdout)
    return json.loads(p.stdout) if p.stdout.strip() else dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)

def shell(pane):
    value=cli('pane','process-info','--pane',pane)['result']['process_info']
    processes=value['foreground_processes']
    if len(processes)==1:
        process=processes[0]
        argv=process.get('argv',[])
        if len(argv)==3 and argv[:2]==[PYTHON,str(A80/'launch_production.py')]:
            previous_path=Path(argv[2]).resolve()
            if previous_path.parent in (HERE/'cases',A80/'cases'):
                previous=json.loads(previous_path.read_text())
                terminal_path=Path(previous['receipts'])/'exit.json'
                if previous['pane_id']==pane and terminal_path.exists():
                    terminal=json.loads(terminal_path.read_text())
                    if terminal['returncode'] is not None:
                        # The frozen launcher writes its manager exit before its final
                        # print/return. Wait for that exact completed transport process,
                        # never for or by retrying a physics case.
                        import select
                        try: fd=os.pidfd_open(process['pid'])
                        except ProcessLookupError: fd=None
                        if fd is not None:
                            try:
                                ready,_,_=select.select([fd],[],[],10)
                                if not ready:
                                    raise TimeoutError('completed launcher did not return within10s: '+str(previous_path))
                            finally:
                                os.close(fd)
                        wait_dir=HERE/'launcher-exit-waits'
                        wait_dir.mkdir(exist_ok=True)
                        value_after=cli('pane','process-info','--pane',pane)['result']['process_info']
                        save(wait_dir/(str(process['pid'])+'.json'),dict(previous_case=binding(previous_path),
                             terminal=binding(terminal_path),before=value,after=value_after))
                        value=value_after
                        processes=value['foreground_processes']
    if len(processes)!=1 or processes[0]['pid']!=value['shell_pid'] or processes[0]['name']!='bash':
        raise RuntimeError('not sole shell foreground: '+json.dumps(value))
    return value

def prepare(component,q,cpu,pane,tab):
    name=f'{component}-q{q}-full-001'
    receipts=HERE/('receipts-'+name)
    receipts.mkdir()
    case=json.loads((A80/'cases/ee-q2-full-001.json').read_text())
    case.update(component=component,q=q,cpu=cpu,cap_seconds=CAP[component],pane_id=pane,tab_id=tab,
        label='P254 '+name+' full31 production',unit='p254-a81-'+name+'.service',
        receipts=str(receipts),output='/var/tmp/p254-a81-'+name,
        allocation=binding(EXTRA_LEASE if cpu==0 else EXTRA_EI_LEASE if cpu==1 else EXTRA_CUBIC_LEASE if cpu==2 else HERE/'preregistration.json'))
    path=HERE/'cases'/(name+'.json')
    path.parent.mkdir(exist_ok=True)
    save(path,case)
    return path,case

def launch(component,q,cpu,pane,tab):
    evidence=shell(pane)
    path,case=prepare(component,q,cpu,pane,tab)
    save(Path(case['receipts'])/'sole-shell-before-launch.json',evidence)
    response=cli('pane','run',pane,f'{PYTHON} {A80}/launch_production.py {path}')
    save(Path(case['receipts'])/'pane-admission.json',response)
    return case

def observe(case):
    # One finite observer is bound to this case's terminal marker and receipt.
    # The generic marker may be historical in a reused pane; exit.json is authoritative.
    receipt=Path(case['receipts'])/'exit.json'
    from datetime import datetime
    launch_receipt=Path(case['receipts'])/'launch.json'
    launched=json.loads(launch_receipt.read_text()) if launch_receipt.exists() else {}
    began=datetime.fromisoformat(launched['actual_start_utc'].replace('Z','+00:00')).timestamp() if 'actual_start_utc' in launched else time.time()
    deadline=time.monotonic()+max(0,began+case['cap_seconds']+120-time.time())
    import select
    import ctypes
    libc=ctypes.CDLL(None,use_errno=True)
    fd=libc.inotify_init1(os.O_NONBLOCK|os.O_CLOEXEC)
    if fd<0:
        raise OSError(ctypes.get_errno(),'inotify_init1')
    try:
        watch=libc.inotify_add_watch(fd,os.fsencode(receipt.parent),0x80|0x8)
        if watch<0:
            raise OSError(ctypes.get_errno(),'inotify_add_watch')
        while True:
            try:
                terminal=json.loads(receipt.read_text())
                break
            except (FileNotFoundError,json.JSONDecodeError):
                pass
            remaining=deadline-time.monotonic()
            if remaining<=0:
                raise TimeoutError('terminal receipt absent beyond finite Native cap: '+case['unit'])
            select.select([fd],[],[],remaining)
            try: os.read(fd,65536)
            except BlockingIOError: pass
        result_path=receipt.parent/'result.json'
        result=json.loads(result_path.read_text()) if result_path.exists() else None
        return dict(case=case,terminal=terminal,result=result,
                    receipt=binding(receipt),result_receipt=binding(result_path) if result_path.exists() else None)
    finally:
        os.close(fd)

def lease_events(wakeup):
    import ctypes
    import select
    libc=ctypes.CDLL(None,use_errno=True)
    fd=libc.inotify_init1(os.O_NONBLOCK|os.O_CLOEXEC)
    if fd<0: raise OSError(ctypes.get_errno(),'lease inotify')
    releases=HERE.parent/'0082/releases'
    watched={}
    def watch(path):
        wd=libc.inotify_add_watch(fd,os.fsencode(path),0x80|0x100|0x8)
        if wd<0: raise OSError(ctypes.get_errno(),'lease watch '+str(path))
        watched[wd]=path
    try:
        watch(HERE)
        watch(HERE.parent/'0082')
        if releases.is_dir(): watch(releases)
        while True:
            select.select([fd],[],[])
            try: events=os.read(fd,65536)
            except BlockingIOError: continue
            offset=0
            while offset<len(events):
                wd,_,_,length=struct.unpack_from('iIII',events,offset)
                name=events[offset+16:offset+16+length].split(b'\0',1)[0]
                offset+=16+length
                if watched.get(wd)==HERE.parent/'0082' and name==b'releases' and releases.is_dir():
                    watch(releases)
                    wakeup.put(None)
                elif watched.get(wd)==releases and name.endswith(b'.json'):
                    wakeup.put(None)
                elif watched.get(wd)==HERE and name in (b'cpu10-lease.json',b'budget-carry.json',b'bank-core-request.json',b'ceiling-amendment.json'):
                    wakeup.put(None)
    finally:
        os.close(fd)

def apply_bank_releases(busy):
    record_path=HERE/'applied-bank-releases.json'
    applied=json.loads(record_path.read_text()) if record_path.exists() else {}
    request_path=HERE/'bank-core-request.json'
    request=json.loads(request_path.read_text()) if request_path.exists() else {}
    for path in sorted((HERE.parent/'0082/releases').glob('*.json')):
        if path.name in applied: continue
        try: release=json.loads(path.read_text())
        except json.JSONDecodeError: continue  # Wait for this release's close-write event.
        if not release['production_authorized']:
            raise RuntimeError('bank release does not authorize production: '+str(path))
        if release['cpu'] in busy:
            raise RuntimeError('bank release conflicts with an owned active Native case: '+str(path))
        acknowledgments=[]
        for name in release['admitted_handoffs']:
            if Path(name).name!=name:
                raise ValueError('release handoff is not a basename')
            ack_path=HERE.parent/'0082/intake'/name
            ack=json.loads(ack_path.read_text())
            if ack['status']!='actual-admitted-consumed-retired-reopened' or ack['handoff']!=binding(HERE/'handoff'/name):
                raise ValueError('bank acknowledgment differs from actual handed contribution: '+name)
            acknowledgments.append(binding(ack_path))
        remaining=[p for p in (HERE/'handoff').glob('*.json') if not (HERE.parent/'0082/intake'/p.name).exists()]
        if request.get('pending') and request.get('cpu')==release['cpu']:
            request.update(pending=bool(remaining),claim_first_free=False,after_release=binding(path),
                reason='Actual durable bank release; follow-on pending admission retained on freed core' if remaining else 'Actual durable bank release and all queued intake acknowledgments complete')
            save(request_path,request)
        lease_path=HERE/'cpu10-lease.json'
        lease=json.loads(lease_path.read_text())
        lease.update(production_authorized=not(request.get('pending') and request.get('cpu')==10),
                     reason='Explicit durable bank release applied; CPU10 is held only when it is the admission core')
        save(lease_path,lease)
        applied[path.name]=dict(release=binding(path),acknowledgments=acknowledgments)
        save(record_path,applied)


def main():
    # CPU10 remains excluded until an explicit cooperative lease appears.
    active={}
    all_events=[]
    costs=dict(OLD)
    # Preserve the exact preregistered aggregate and charge paid q1 EI separately.
    prior_ei=HERE.parent/'0071/handoff.json'
    paid_ei=json.loads(prior_ei.read_text())['cost']['launch_elapsed_seconds']
    costs['mixed_ei']+=paid_ei
    save(HERE/'historical-ei-costs.json',dict(base=OLD['mixed_ei'],paid_q1_receipt=binding(prior_ei),
         paid_q1_seconds=paid_ei,total_numeric_known=costs['mixed_ei'],
         unresolved_prior_history=CONTRACT['known_spent_before_first_wave']['mixed_ei_prior_A76'],
         note='All original receipts retained. Root aggregate correction is additive through budget-carry.json.'))
    inherited=[]
    for name in ('ee-q0-remainder-001','ee-q2-full-001','mixed_ei-q0-remainder-001','cubic-q0-remainder-001','ee-q1-pack-001'):
        case=json.loads((A80/'cases'/(name+'.json')).read_text())
        if case['operation']=='pack-q1':
            case.update(component='ee',q=1)
        inherited.append(case)
    inherited.extend(json.loads(p.read_text()) for p in sorted((HERE/'cases').glob('*.json')))
    queues={'ee':list(range(3,64)),'mixed_ei':list(range(2,64)),'cubic':list(range(2,64))}
    slots={cpu:pane for cpu,pane in BANK_SLOTS.items() if cpu!=10}
    slots[0]=EXTRA_PANE
    slots[1]=EXTRA_EI_PANE
    # Root625 holds new CPU2 admission; inherited CPU2 Native remains observed.
    extra_allocation=binding(EXTRA_LEASE)
    extra_dispatch_count=sum(c['cpu']==0 and c.get('allocation')==extra_allocation for c in inherited)
    extra_dispatch_blocked=False
    extra_ei_allocation=binding(EXTRA_EI_LEASE)
    extra_ei_dispatch_count=sum(c['cpu']==1 and c.get('allocation')==extra_ei_allocation for c in inherited)
    extra_ei_dispatch_blocked=False
    extra_cubic_allocation=binding(EXTRA_CUBIC_LEASE)
    extra_cubic_dispatch_count=sum(c['cpu']==2 and c.get('allocation')==extra_cubic_allocation for c in inherited)
    extra_cubic_dispatch_blocked=False
    # One watcher per possible active physics core, including all bounded loans.
    pool=concurrent.futures.ThreadPoolExecutor(max_workers=len(BANK_SLOTS)+3)
    wakeup=queue.Queue()
    threading.Thread(target=lease_events,args=(wakeup,),daemon=True).start()
    def submit(case):
        future=pool.submit(observe,case)
        active[future]=case
        future.add_done_callback(wakeup.put)
    for case in inherited:
        submit(case)
    existing={(c['component'],c['q']) for c in inherited}
    for component in queues:
        queues[component]=[q for q in queues[component] if (component,q) not in existing]
    for component,q,cpu in [('ee',3,8),('mixed_ei',2,9)]:
        if (component,q) not in existing:
            case=launch(component,q,cpu,*slots[cpu]); queues[component].remove(q)
            submit(case)
    turn=0
    def ledger():
        carry=HERE/'budget-carry.json'
        carried=json.loads(carry.read_text()) if carry.exists() else {}
        reservations={k:sum(c['cap_seconds'] for c in active.values() if c['component']==k) for k in costs}
        save(HERE/'actual-cost-ledger.json',dict(cost_seconds=costs,active_cap_reservations=reservations,
             ceilings=ceilings(),historical_costs_preserved=CONTRACT['known_spent_before_first_wave'],
             additive_historical_carry=carried,
             events=all_events,active=list(active.values()),remaining_queues=queues,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    ledger()
    while active or any(queues.values()):
        if active:
            notification=wakeup.get()
            done=set() if notification is None else {notification}
            for future in done:
                case=active.pop(future)
                try: value=future.result()
                except BaseException as error:
                    value=dict(case=case,observer_failure=repr(error))
                    save(HERE/('blocked-'+case['unit']+'.json'),value)
                    all_events.append(value); ledger()
                    # No successor on a core without authoritative terminal ownership.
                    slots.pop(case['cpu'],None)
                    continue
                result=value['result']
                charge=value['terminal']['elapsed_seconds']
                costs[case['component']]+=charge
                value['actual_charged_seconds']=charge
                all_events.append(value)
                stem=case['unit'].removesuffix('.service')
                save(HERE/('settled-'+stem+'.json'),value)
                if value['terminal']['returncode']==0 and result and result['status']=='complete-owned-production' and case['operation']=='component' and (case['unit'].startswith('p254-a81-') or case['component']=='ee'):
                    catalog=json.loads(Path(result['catalog']['path']).read_text())
                    value['eligible_staging_arrays']=[v for v in catalog.get('sigma',catalog.get('arrays',{})).values() if str(v['path']).startswith(('/var/tmp/p254-a80-','/var/tmp/p254-a81-'))]
                    value['retirement_scope']='Only bound own-new A80/A81 sigma/arrays, after assembler commit/reopen; never original inputs.'
                    ready=HERE/'handoff'; ready.mkdir(exist_ok=True)
                    save(ready/(stem+'.json'),value)
                    lease_path=HERE/'cpu10-lease.json'
                    if lease_path.exists() and not (HERE.parent/'0082/intake'/(stem+'.json')).exists():
                        lease_value=json.loads(lease_path.read_text())
                        lease_value.update(production_authorized=False,reason='Pending settled component admission',pending_handoff=str(ready/(stem+'.json')))
                        save(lease_path,lease_value)
                        request_path=HERE/'bank-core-request.json'
                        request=json.loads(request_path.read_text()) if request_path.exists() else {}
                        if not request.get('pending'):
                            bank_cpu=case['cpu'] if case['cpu'] in BANK_SLOTS else 6
                            request.update(cpu=bank_cpu,pending=True,
                                claim_first_free=case['cpu'] not in BANK_SLOTS,
                                pane_id=BANK_SLOTS[bank_cpu][0],tab_id=BANK_SLOTS[bank_cpu][1],
                                reason='Actual settled component: first free ORIGINAL bank core before successor physics',
                                authorization='Root original6–10 bank admission policy; extra producer CPUs never bank cores')
                            save(request_path,request)
                elif value['terminal']['returncode']!=0 or not result or result['status']!='complete-owned-production':
                    save(HERE/('failed-'+stem+'.json'),value)
                    if case['cpu']==0 and case.get('allocation')==extra_allocation: extra_dispatch_blocked=True
                    if case['cpu']==1 and case.get('allocation')==extra_ei_allocation: extra_ei_dispatch_blocked=True
                    if case['cpu']==2 and case.get('allocation')==extra_cubic_allocation: extra_cubic_dispatch_blocked=True
        busy={c['cpu'] for c in active.values()}
        apply_bank_releases({c['cpu'] for c in active.values() if not (Path(c['receipts'])/'exit.json').exists()})
        save(HERE/'cpu10-state.json',dict(native_busy=10 in busy,active_case=next((c for c in active.values() if c['cpu']==10),None),actual_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        bank_request_path=HERE/'bank-core-request.json'
        bank_request=json.loads(bank_request_path.read_text()) if bank_request_path.exists() else {}
        if bank_request.get('pending') and bank_request.get('claim_first_free'):
            free=[cpu for cpu in sorted(BANK_SLOTS) if cpu not in busy]
            if free:
                cpu=free[0]
                pane,tab=BANK_SLOTS[cpu]
                bank_request.update(cpu=cpu,pane_id=pane,tab_id=tab,claim_first_free=False)
                save(bank_request_path,bank_request)
        held_cpu=bank_request.get('cpu') if bank_request.get('pending') else None
        if held_cpu is not None:
            save(HERE/'bank-core-state.json',dict(cpu=held_cpu,native_busy=held_cpu in busy,
                active_case=next((c for c in active.values() if c['cpu']==held_cpu),None),
                pane_id=bank_request['pane_id'],tab_id=bank_request['tab_id'],
                request_binding=binding(bank_request_path),
                boundary='Native terminal receipt required; sole shell must be freshly observed before admission'))
        lease=HERE/'cpu10-lease.json'
        if lease.exists() and json.loads(lease.read_text()).get('production_authorized'):
            lease_value=json.loads(lease.read_text())
            slots[10]=(lease_value['pane_id'],lease_value['tab_id'])
        elif 10 not in busy:
            slots.pop(10,None)
        for cpu,(pane,tab) in sorted(slots.items()):
            if cpu in busy or cpu==held_cpu: continue
            if cpu==0 and (extra_dispatch_blocked or extra_dispatch_count>=EXTRA_CASE_QUOTA or time.time()>=EXTRA_DISPATCH_DEADLINE): continue
            if cpu==1 and (extra_ei_dispatch_blocked or extra_ei_dispatch_count>=EXTRA_EI_CASE_QUOTA or time.time()>=EXTRA_EI_DISPATCH_DEADLINE): continue
            if cpu==2 and (extra_cubic_dispatch_blocked or extra_cubic_dispatch_count>=EXTRA_CUBIC_CASE_QUOTA or time.time()>=EXTRA_CUBIC_DISPATCH_DEADLINE): continue
            preference=['ee','mixed_ei','cubic'] if turn%3!=2 else ['cubic','ee','mixed_ei']
            if cpu==9: preference=['mixed_ei','cubic','ee'] if turn%2==0 else ['cubic','mixed_ei','ee']
            if cpu==0: preference=['ee']  # Root618: bounded distinct EE cases.
            if cpu==1: preference=['mixed_ei']  # Root621: bounded distinct required EI cases.
            if cpu==2: preference=['cubic']  # Root623: bounded distinct required cubic cases.
            selected=None
            carried=json.loads((HERE/'budget-carry.json').read_text()) if (HERE/'budget-carry.json').exists() else {}
            for component in preference:
                reserved=sum(c['cap_seconds'] for c in active.values() if c['component']==component)
                if queues[component] and costs[component]+carried.get(component,0)+reserved+CAP[component]<=ceilings()[component]:
                    selected=component; break
            if selected is None: continue
            q=queues[selected][0]
            allocated,unallocated=output_reservations(active.values())
            needed=OUTPUT_NPY_BYTES[selected]
            if allocated+unallocated+needed>OUTPUT_LIMIT_BYTES:
                save(HERE/'backpressure.json',dict(unadmitted_array_bytes=allocated,
                     unallocated_active_reservation_bytes=unallocated,
                     requested_case_bytes=needed,limit=OUTPUT_LIMIT_BYTES))
                continue
            try: case=launch(selected,q,cpu,pane,tab)
            except BaseException as error:
                save(HERE/('launch-blocked-'+selected+'-q'+str(q)+'-cpu'+str(cpu)+'.json'),dict(component=selected,q=q,error=repr(error)))
                slots.pop(cpu,None)
                continue
            queues[selected].pop(0)
            submit(case)
            if cpu==0: extra_dispatch_count+=1
            if cpu==1: extra_ei_dispatch_count+=1
            if cpu==2: extra_cubic_dispatch_count+=1
            turn+=1
        ledger()
        if not active and any(queues.values()):
            save(HERE/'dispatch-blocked.json',dict(remaining=queues,cost_seconds=costs,reason='No admissible core/storage/budget; no automatic retry.'))
            break
    ledger()
    print(json.dumps(dict(status='terminal-owned-dispatch',remaining=queues,cost_seconds=costs)),flush=True)

if __name__=='__main__':
    main()
