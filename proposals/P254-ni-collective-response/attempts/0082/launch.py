"""One explicitly owned finite admission, retaining raw actual Native receipts."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
HERE = Path(__file__).resolve().parent


def main(case_path):
    case = json.loads(Path(case_path).read_bytes())
    receipts = Path(case['receipts'])
    context = {k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
    if context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != case['pane_id']:
        raise ValueError('launch must originate in owned real pane')
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    worker = case.get('worker','admit.py')
    if worker not in ('admit.py','retire.py','consume_paid_nine.py'):
        raise ValueError('only owned admission/retirement workers allowed')
    cpu_id = case.get('cpu',10)
    if cpu_id not in range(6,11):
        raise ValueError('protected CPU')
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+case['unit'],'--property=RuntimeMaxSec='+str(case['cap_seconds']),'--property=MemoryMax=2147483648','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity='+str(cpu_id),'--property=StartLimitBurst=1','env',contract['python_environment_assignment'],*[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in context.items()],'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),contract['python'],'-u',str(HERE/worker),str(Path(case_path).resolve())]
    def save(name,value):
        with (receipts/name).open('x') as stream:
            json.dump(value,stream,indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
    began = time.monotonic()
    save('launch.json',dict(argv=argv,context=context,actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),automatic_retries=0))
    code = None
    try:
        with (receipts/'stdout.txt').open('xb') as stdout,(receipts/'stderr.txt').open('xb') as stderr:
            process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line)
                stdout.flush()
                print(line.decode(errors='replace'),end='',flush=True)
            code = process.wait()
    finally:
        manager = contract['manager_env']+['systemctl','--user','show',case['unit'],'--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal = subprocess.run(manager,capture_output=True,check=False)
        save('manager-terminal.json',dict(argv=manager,returncode=terminal.returncode,stdout=terminal.stdout.decode(),stderr=terminal.stderr.decode()))
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-began,actual_end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        print('A82_NATIVE_TERMINAL_'+str(code)+' '+case['unit'],flush=True)
    return code if code is not None else 1


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1]))
