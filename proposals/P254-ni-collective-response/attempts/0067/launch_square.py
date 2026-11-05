"""Bounded actual full-history owned block, with live and durable Native output."""
from pathlib import Path
import json
import os
import subprocess
import time

HERE=Path(__file__).resolve().parent


def main():
    contract=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    receipts=HERE/'receipts-001';receipts.mkdir(exist_ok=False)
    output=Path('/var/tmp/p254-a67-ee-full31-001')
    context={key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    if context.get('HERDR_ENV')!='1' or not all(context.get(key) for key in ('HERDR_PANE_ID','HERDR_TAB_ID','HERDR_WORKSPACE_ID')):
        raise ValueError('verified real Herdr owner context required')
    unit='p254-a67-ee-full31-001.service'
    argv=contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec=360','--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=6',
        'env',contract['python_environment_assignment'],*[key+'='+value for key,value in contract['threads'].items()],
        *[key+'='+value for key,value in context.items()],'PYTHONDONTWRITEBYTECODE=1',
        '/usr/bin/time','-v','-o',str(receipts/'time.txt'),contract['python'],'-u',
        str(HERE/'produce_square.py'),'--output',str(output)]
    def save(name,value):(receipts/name).write_text(json.dumps(value,indent=2)+'\n')
    started=time.monotonic();code=None
    save('launch.json',dict(argv=argv,unit=unit,actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    print('A67_NATIVE_START '+unit,flush=True)
    try:
        with (receipts/'stdout').open('xb') as stdout,(receipts/'stderr').open('xb') as stderr:
            process=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line);stdout.flush();print(line.decode(errors='replace'),end='',flush=True)
            code=process.wait()
    finally:
        manager=contract['manager_env']+['systemctl','--user','show',unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal=subprocess.run(manager,capture_output=True,check=False)
        save('manager-terminal.json',dict(argv=manager,returncode=terminal.returncode,stdout=terminal.stdout.decode(),stderr=terminal.stderr.decode()))
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-started))
        print('A67_TERMINAL_'+str(code),flush=True)
    return code if code is not None else 1


if __name__=='__main__':raise SystemExit(main())
