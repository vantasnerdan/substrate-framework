"""One bounded invocation, all-exit receipts, real inherited Herdr context."""
from pathlib import Path
import os,json,subprocess,time
HERE=Path(__file__).resolve().parent
contract=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
RECEIPTS=HERE/'receipts-001';RECEIPTS.mkdir(exist_ok=False)
unit='p254-a64-gaussian-ee-001.service'
assert os.environ.get('HERDR_ENV')=='1'
herdr={k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
assert all(herdr.get(k) for k in ('HERDR_PANE_ID','HERDR_TAB_ID','HERDR_WORKSPACE_ID'))
argv=contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
    '--property=RuntimeMaxSec=590','--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=6',
    'env',contract['python_environment_assignment'],*[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in herdr.items()],
    'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(RECEIPTS/'time.txt'),contract['python'],'-u',str(HERE/'material_pilot.py')]
start=time.monotonic()
def save(name,value):(RECEIPTS/name).write_text(json.dumps(value,indent=2)+'\n')
save('launch.json',dict(argv=argv,unit=unit,start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),symbolic_prior_wall_seconds=1.56,native_cap_seconds=600,herdr_context=herdr))
print('A64_OWNED_NATIVE_BEGIN '+unit,flush=True)
rc=None
try:
    with (RECEIPTS/'stdout.txt').open('xb') as out,(RECEIPTS/'stderr.txt').open('xb') as err:
        rc=subprocess.run(argv,stdout=out,stderr=err,check=False).returncode
finally:
    manager=contract['manager_env']+['systemctl','--user','show',unit,'--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
    result=subprocess.run(manager,capture_output=True,check=False)
    save('manager-terminal.json',dict(argv=manager,returncode=result.returncode,stdout=result.stdout.decode(),stderr=result.stderr.decode()))
    save('exit.json',dict(argv=argv,returncode=rc,elapsed_seconds=time.monotonic()-start,symbolic_prior_wall_seconds=1.56,cumulative_this_attempt_upper_seconds=time.monotonic()-start+1.56))
    print('A64_OWNED_LAUNCH_TERMINAL rc='+str(rc),flush=True)
