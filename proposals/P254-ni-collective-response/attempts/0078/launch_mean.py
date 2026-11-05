"""Producer <=780s plus actual saved consumer <=120s: total cap <=900s."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
HERE = Path(__file__).resolve().parent


def main(mode):
    if mode not in ('producer','consumer'):
        raise ValueError('unknown execution mode')
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    context = {k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
    if context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != 'w1:p7K':
        raise ValueError('actual assigned owner view required')
    receipts = HERE/('receipts-'+mode)
    receipts.mkdir(exist_ok=False)
    cap = 780 if mode == 'producer' else 120
    unit = 'p254-a78-fullhistory-mean-'+mode+'-001.service'
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec='+str(cap),'--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=6',
        'env',contract['python_environment_assignment'],
        *[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in context.items()],
        'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),
        contract['python'],'-u',str(HERE/'native_mean.py'),mode]
    def save(name,value):
        (receipts/name).write_text(json.dumps(value,indent=2)+'\n')
    started = time.monotonic()
    save('launch.json',dict(argv=argv,unit=unit,cap_seconds=cap,context=context,
        actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    code = None
    try:
        with (receipts/'stdout').open('xb') as stdout,(receipts/'stderr').open('xb') as stderr:
            process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line);stdout.flush()
                print(line.decode(errors='replace'),end='',flush=True)
            code = process.wait()
    finally:
        manager = contract['manager_env']+['systemctl','--user','show',unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal = subprocess.run(manager,capture_output=True,check=False)
        save('manager-terminal.json',dict(argv=manager,returncode=terminal.returncode,
            stdout=terminal.stdout.decode(),stderr=terminal.stderr.decode()))
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-started,
            actual_end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        print('A78_'+mode.upper()+'_OWNER_TERMINAL_'+str(code),flush=True)
    return code if code is not None else 1


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1]))
