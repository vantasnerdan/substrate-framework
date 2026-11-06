"""Finite all64 queue owner; orchestration only, unchanged Native numerical workers.

Each case is launched once, from the leased real pane, with existing bounded
launch/admit/retire helpers. Any error stops the controller without retries.
"""
from pathlib import Path
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
HERE = Path(__file__).resolve().parent
ATTEMPTS = HERE.parent
PRODUCER = ATTEMPTS/'0081'
BANK = ATTEMPTS/'0077/exercise-002/bank'
PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'


def binding(path):
    path=Path(path).resolve()
    payload=path.read_bytes()
    return dict(path=str(path),bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest())


def authenticated(descriptor):
    if binding(descriptor['path']) != {k:descriptor[k] for k in ('path','bytes','sha256')}:
        raise ValueError('changed metadata binding: '+descriptor['path'])
    return Path(descriptor['path'])


def load(path):
    return json.loads(Path(path).read_bytes())


def save_new(path,record):
    path=Path(path)
    path.parent.mkdir(exist_ok=True)
    with path.open('x') as stream:
        json.dump(record,stream,indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def command(argv,receipt):
    started=time.monotonic()
    result=subprocess.run(argv,capture_output=True,text=True,check=False)
    save_new(receipt,dict(argv=argv,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,elapsed_seconds=time.monotonic()-started))
    if result.returncode:
        raise RuntimeError('orchestration command failed: '+str(receipt))
    return result


def sole_shell(pane,path):
    result=command(['herdr','pane','process-info','--pane',pane],path)
    info=json.loads(result.stdout)['result']['process_info']
    processes=info['foreground_processes']
    if len(processes)!=1 or processes[0]['pid']!=info['shell_pid'] or processes[0]['name']!='bash' or info['pane_id']!=pane:
        raise ValueError('leased pane is not a sole foreground shell')
    return info


def run_native(case_path,logdir):
    case=load(case_path)
    text=shlex.join([PYTHON,str(HERE/'launch.py'),str(case_path)])
    command(['herdr','pane','run',case['pane_id'],text],logdir/'pane-run.json')
    # Unit-specific terminal marker prevents an old successful case satisfying it.
    command(['herdr','pane','wait-output',case['pane_id'],'--regex','A82_NATIVE_TERMINAL_.*'+case['unit'],'--source','recent','--timeout',str((case['cap_seconds']+20)*1000),'--raw'],logdir/'terminal-observer.json')
    terminal=load(Path(case['receipts'])/'exit.json')
    if terminal['returncode']!=0:
        raise ValueError('actual Native failed; no retry: '+case['unit'])


def prepare(batch,lease,paths):
    handoffs=[]
    for path in paths:
        record=load(path)
        authenticated(record['receipt'])
        authenticated(record['result_receipt'])
        authenticated(record['result']['consumer'])
        if record['terminal']['returncode']!=0 or record['result']['status']!='complete-owned-production':
            raise ValueError('unsettled producer handoff')
        if (HERE/'intake'/Path(path).name).exists():
            raise ValueError('already acknowledged handoff in new batch')
        handoffs.append((Path(path),record))
    receipts=HERE/f'receipts-intake-{batch:03d}'
    receipts.mkdir(exist_ok=False)
    case=dict(schema='p254-a82-owned-admission-case-v1',cpu=lease['cpu'],pane_id=lease['pane_id'],unit=f'p254-a82-bank-intake-{batch:03d}.service',cap_seconds=360,receipts=str(receipts),bank=str(BANK),publication=str(HERE/f'partial-intake-{batch:03d}.json'),canonical_output='/var/tmp/p254-a82-complete-bare-map-001/catalog.json',catalogs=[h['result']['catalog']['path'] for p,h in handoffs],settled_handoffs=[binding(p) for p,h in handoffs],source_bindings=[binding(HERE/'admit.py'),binding(HERE/'launch.py'),binding(ATTEMPTS/'0077/component_accumulator.py')])
    path=HERE/f'case-intake-{batch:03d}.json'
    save_new(path,case)
    return path,handoffs


def prepare_retirement(batch,lease,handoffs):
    admission=HERE/f'receipts-intake-{batch:03d}'
    items=[]
    for index,(path,h) in enumerate(handoffs):
        consumer=load(admission/f'{index}-consumer.json')
        if consumer['status']!='actual-saved-bank-consumer-pass':
            raise ValueError('actual saved-bank consumer missing')
        unit=h['case']['unit'].removesuffix('.service')
        items.append(dict(producer_exit=h['receipt'],producer_result=h['result_receipt'],admission_exit=binding(admission/'exit.json'),bank_consumer=binding(admission/f'{index}-consumer.json'),receipt=str(HERE/f'retired-{unit}.json'),handoff=binding(path)))
    receipts=HERE/f'receipts-retire-intake-{batch:03d}'
    receipts.mkdir(exist_ok=False)
    spec=dict(schema='p254-a82-owned-retirement-case-v1',worker='retire.py',cpu=lease['cpu'],pane_id=lease['pane_id'],unit=f'p254-a82-retire-intake-{batch:03d}.service',cap_seconds=100,receipts=str(receipts),bank=str(BANK),items=items,result=str(receipts/'result.json'))
    path=HERE/f'case-retire-intake-{batch:03d}.json'
    save_new(path,spec)
    return path


def acknowledge(batch,handoffs):
    admission=HERE/f'receipts-intake-{batch:03d}'
    retirement=HERE/f'receipts-retire-intake-{batch:03d}'
    if load(admission/'exit.json')['returncode']!=0 or load(retirement/'exit.json')['returncode']!=0:
        raise ValueError('acknowledgment before actual Native0')
    ack_paths=[]
    for index,(path,h) in enumerate(handoffs):
        unit=h['case']['unit'].removesuffix('.service')
        retired=HERE/f'retired-{unit}.json'
        record=load(retired)
        if record['status']!='retired-and-actual-reopen-pass':
            raise ValueError('retirement/reopen gate missing')
        ack=dict(schema='p254-a82-settled-intake-ack-v1',status='actual-admitted-consumed-retired-reopened',handoff=binding(path),catalog=h['result']['catalog'],producer_terminal=h['receipt'],producer_savedconsumer=h['result']['consumer'],actual_producer_charged_seconds=h['actual_charged_seconds'],admission=binding(admission/f'{index}-admission.json'),saved_bank_consumer=binding(admission/f'{index}-consumer.json'),admission_result=binding(admission/'result.json'),admission_terminal=binding(admission/'exit.json'),retirement=binding(retired),retirement_terminal=binding(retirement/'exit.json'),committed_generation_snapshot=record['committed_generation_snapshot'],recovery_loss=record['recovery_loss'])
        ack_path=HERE/'intake'/path.name
        save_new(ack_path,ack)
        ack_paths.append(ack_path)
    return ack_paths


def main(first_batch,after_release,seconds):
    if os.environ.get('HERDR_ENV')!='1':
        raise ValueError('actual inherited Herdr orchestration context required')
    began=time.monotonic()
    controller=HERE/f'controller-{first_batch:03d}'
    controller.mkdir(exist_ok=False)
    save_new(controller/'startup.json',dict(pid=os.getpid(),first_batch=first_batch,after_release=binding(after_release),deadline_seconds=seconds,context={k:v for k,v in os.environ.items() if k.startswith('HERDR_')},source=binding(__file__),numeric_workers='unchanged existing Native helpers',automatic_retries=0))
    batch=first_batch
    last_release=Path(after_release)
    try:
        while time.monotonic()-began < seconds:
            step=controller/f'{batch:03d}'
            step.mkdir(exist_ok=False)
            remaining=seconds-(time.monotonic()-began)
            result=command([PYTHON,str(HERE/'wait_lease.py'),str(last_release),str(min(2700,remaining))],step/'lease-observer.json')
            lease=json.loads(result.stdout)
            if lease['event']!='actual-new-bank-lease-ready':
                raise ValueError('finite observer did not deliver a lease')
            request=lease['request']
            # Explicitly bound current producer generation, not merely stale CPU values.
            if binding(PRODUCER/'bank-core-request.json')!=lease['state']['request_binding']:
                raise ValueError('lease request changed before admission')
            if request['cpu'] not in range(6,11):
                raise ValueError('protected CPU')
            sole_shell(request['pane_id'],step/'before-sole-shell.json')
            # Six catalogs bound one finite batch; excess stays owned in durable queue.
            paths=lease['pending_handoffs'][:6]
            case,handoffs=prepare(batch,request,paths)
            admission_log=step/'admission'
            admission_log.mkdir()
            run_native(case,admission_log)
            retirement=prepare_retirement(batch,request,handoffs)
            retirement_log=step/'retirement'
            retirement_log.mkdir()
            run_native(retirement,retirement_log)
            acks=acknowledge(batch,handoffs)
            shell=sole_shell(request['pane_id'],step/'after-sole-shell.json')
            release=dict(schema='p254-a82-terminal-core-release-v1',cpu=request['cpu'],pane_id=request['pane_id'],tab_id=request['tab_id'],admitted_handoffs=[p.name for p,h in handoffs],production_authorized=True,sole_shell_evidence=shell,intake_acks=[binding(p) for p in acks])
            last_release=HERE/'releases'/f'intake-{batch:03d}.json'
            save_new(last_release,release)
            actual=load(HERE/f'receipts-intake-{batch:03d}/result.json')
            event=dict(event='actual-batch-admitted-retired-reopened-released',batch=batch,coverage=actual['coverage'],result=binding(HERE/f'receipts-intake-{batch:03d}/result.json'),release=binding(last_release),intake_acks=[binding(p) for p in acks])
            save_new(step/'complete.json',event)
            print(json.dumps(event),flush=True)
            if all(x['complete'] for x in actual['coverage'].values()):
                if not actual.get('complete_A62_consumer'):
                    raise ValueError('complete coverage without actual final A62 consumer')
                save_new(controller/'terminal.json',dict(status='actual-all64-complete-A62-consumed',batch=batch,actual_result=binding(HERE/f'receipts-intake-{batch:03d}/result.json'),canonical=actual['publication'],consumer=actual['complete_A62_consumer'],last_release=binding(last_release),elapsed_seconds=time.monotonic()-began))
                print('A82_ROLLING_ALL64_COMPLETE',flush=True)
                return 0
            batch+=1
        raise TimeoutError('finite rolling owner deadline; no Native renewal')
    except BaseException as error:
        save_new(controller/'terminal.json',dict(status='stopped-no-automatic-retry',batch=batch,last_release=binding(last_release),error_type=type(error).__name__,error=str(error),elapsed_seconds=time.monotonic()-began))
        print('A82_ROLLING_STOPPED '+str(error),flush=True)
        raise


if __name__=='__main__':
    raise SystemExit(main(int(sys.argv[1]),sys.argv[2],float(sys.argv[3])))
