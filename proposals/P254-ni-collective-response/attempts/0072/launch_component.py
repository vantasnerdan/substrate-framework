"""One bounded Native invocation, with exact argv/stdout/stderr/all exits."""
from pathlib import Path
import json
import os
import subprocess
import time
HERE = Path(__file__).resolve().parent
OUTPUT = Path('/var/tmp/p254-a72-cubic-full-history-q1-001')
UNIT = 'p254-a72-cubic-full-history-q1-001.service'


def launch():
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    receipts = HERE/'receipts-001'
    receipts.mkdir(exist_ok=False)
    context = {k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
    if context.get('HERDR_ENV') != '1' or not all(context.get(k) for k in ('HERDR_PANE_ID','HERDR_TAB_ID','HERDR_WORKSPACE_ID')):
        raise ValueError('actual owned Herdr context required')
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+UNIT,
        '--property=RuntimeMaxSec=600','--property=MemoryMax=3221225472','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=9','--property=StartLimitBurst=1',
        'env',contract['python_environment_assignment'],
        *[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in context.items()],
        'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),
        contract['python'],'-u',str(HERE/'produce_component.py'),'--output',str(OUTPUT)]
    def save(name,value):
        (receipts/name).write_text(json.dumps(value,indent=2)+'\n')
    started = time.monotonic()
    save('launch.json',dict(argv=argv,unit=UNIT,cap_seconds=600,cpu=9,context=context,
        actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),automatic_retries=0))
    print('A72_NATIVE_START '+UNIT,flush=True)
    code = None
    try:
        with (receipts/'stdout').open('xb') as stdout,(receipts/'stderr').open('xb') as stderr:
            process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line)
                stdout.flush()
                print(line.decode(errors='replace'),end='',flush=True)
            code = process.wait()
    finally:
        manager = contract['manager_env']+['systemctl','--user','show',UNIT,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal = subprocess.run(manager,capture_output=True,check=False)
        save('manager-terminal.json',dict(argv=manager,returncode=terminal.returncode,
            stdout=terminal.stdout.decode(),stderr=terminal.stderr.decode()))
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-started,
            actual_end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        print('A72_TERMINAL_'+str(code),flush=True)
    return code if code is not None else 1


if __name__ == '__main__':
    raise SystemExit(launch())
