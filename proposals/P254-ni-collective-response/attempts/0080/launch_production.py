"""Launch one finite case from its real Herdr pane; retain actual manager I/O."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def launch(case_path):
    case_path = Path(case_path).resolve()
    case = json.loads(case_path.read_bytes())
    context = {key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    if context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != case['pane_id']:
        raise RuntimeError('launch must originate in the assigned real pane')
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    receipts = Path(case['receipts'])
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+case['unit'],
        '--property=RuntimeMaxSec='+str(case['cap_seconds']),
        '--property=MemoryMax='+str(case['memory_bytes']),'--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity='+str(case['cpu']),'env',
        contract['python_environment_assignment'],*[k+'='+v for k,v in contract['threads'].items()],
        *[k+'='+v for k,v in context.items()],'PYTHONDONTWRITEBYTECODE=1',
        '/usr/bin/time','-v','-o',str(receipts/'time.txt'),contract['python'],'-u',
        str(HERE/'run_production.py'),str(case_path)]
    # Exclusive intent is a single-case duplicate-launch guard, not a retry loop.
    with (receipts/'launch.json').open('x') as stream:
        json.dump(dict(argv=argv,case_path=str(case_path),actual_start_utc=time.strftime(
            '%Y-%m-%dT%H:%M:%SZ',time.gmtime()),context=context),stream,indent=2)
        stream.write('\n')
    began = time.monotonic()
    rc = None
    try:
        with (receipts/'stdout.txt').open('wb') as stdout, (receipts/'stderr.txt').open('wb') as stderr:
            process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in iter(process.stdout.readline,b''):
                stdout.write(line)
                stdout.flush()
                sys.stdout.buffer.write(line)
                sys.stdout.buffer.flush()
            rc = process.wait()
    finally:
        save(receipts/'exit.json',dict(argv=argv,returncode=rc,elapsed_seconds=time.monotonic()-began))
        print('A80_MANAGER_TERMINAL_'+str(rc if rc is not None else -1),flush=True)
    if rc:
        raise SystemExit(rc)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: launch_production.py ACTUAL_BOUND_CASE.json')
    launch(sys.argv[1])
