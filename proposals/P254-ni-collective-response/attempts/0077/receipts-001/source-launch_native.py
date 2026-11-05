"""Exactly one CPU9/3GiB/240s Native exercise; durable argv/output/exit."""
from pathlib import Path
import json
import os
import subprocess
import time

HERE = Path(__file__).resolve().parent
UNIT = 'p254-a77-serial-accumulator-001.service'


def launch():
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    receipts = HERE/'receipts-001'
    receipts.mkdir(exist_ok=False)
    context = {k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
    if context.get('HERDR_ENV') != '1' or not all(context.get(k) for k in ('HERDR_PANE_ID','HERDR_TAB_ID','HERDR_WORKSPACE_ID')):
        raise ValueError('actual owned Herdr context required')
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+UNIT,
        '--property=RuntimeMaxSec=240','--property=MemoryMax=3221225472','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=9','--property=StartLimitBurst=1',
        'env',contract['python_environment_assignment'],
        *[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in context.items()],
        'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),
        contract['python'],'-u',str(HERE/'native_exercise.py')]
    def save(name,value):
        with (receipts/name).open('x') as handle:
            json.dump(value,handle,indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
    started = time.monotonic()
    save('launch.json',dict(argv=argv,unit=UNIT,cap_seconds=240,cpu=9,memory_max=3221225472,
        context=context,actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),automatic_retries=0))
    print('A77_NATIVE_START '+UNIT,flush=True)
    code = None
    try:
        with (receipts/'stdout').open('xb') as stdout,(receipts/'stderr').open('xb') as stderr:
            process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line)
                stdout.flush()
                print(line.decode(errors='replace'),end='',flush=True)
            code = process.wait()
            stdout.flush()
            os.fsync(stdout.fileno())
            stderr.flush()
            os.fsync(stderr.fileno())
    finally:
        manager = contract['manager_env']+['systemctl','--user','show',UNIT,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal = subprocess.run(manager,capture_output=True,check=False)
        save('manager-terminal.json',dict(argv=manager,returncode=terminal.returncode,
            stdout=terminal.stdout.decode(),stderr=terminal.stderr.decode()))
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-started,
            actual_end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        print('A77_TERMINAL_'+str(code),flush=True)
    return code if code is not None else 1


if __name__ == '__main__':
    raise SystemExit(launch())
