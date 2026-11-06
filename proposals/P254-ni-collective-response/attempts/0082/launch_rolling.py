"""Owned Herdr orchestration lifecycle; no new numerical work or Native policy."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent


def main(first,release,deadline):
    context={k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
    if context.get('HERDR_ENV')!='1' or context.get('HERDR_PANE_ID')!='w1:p7V':
        raise ValueError('actual owned labeled runner pane required')
    receipts=HERE/f'controller-launch-{first:03d}'
    receipts.mkdir(exist_ok=False)
    runtime=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    argv=['env',runtime['python_environment_assignment'],*[k+'='+v for k,v in runtime['threads'].items()],'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),runtime['python'],'-u',str(HERE/'rolling_intake.py'),str(first),str(Path(release).resolve()),str(deadline)]
    def save(name,value):
        with (receipts/name).open('x') as stream:
            json.dump(value,stream,indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
    began=time.monotonic()
    save('launch.json',dict(argv=argv,context=context,actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),orchestration_only=True,numerical_workers='Existing unchanged0082 Native launch/admit/retire; leasedCPU6–10/2GiB/noSwap/one-thread/noRestart and existingcaps',automatic_retries=0))
    code=None
    try:
        with (receipts/'stdout.txt').open('xb') as stdout,(receipts/'stderr.txt').open('xb') as stderr:
            process=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line)
                stdout.flush()
                print(line.decode(errors='replace'),end='',flush=True)
            code=process.wait()
    finally:
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-began,actual_end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        print('A82_ROLLING_CONTROLLER_TERMINAL_'+str(code),flush=True)
    return code if code is not None else 1


if __name__=='__main__':
    raise SystemExit(main(int(sys.argv[1]),sys.argv[2],int(sys.argv[3])))
