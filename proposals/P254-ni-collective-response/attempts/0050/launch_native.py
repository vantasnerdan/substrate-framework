"""Launch one bounded CPU9 manager, retaining original terminal evidence."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
BASE=Path('/var/tmp/p254-a50-joint-gamma')

def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())


def main(generation,mode):
    receipts=Path(str(BASE)+'-'+generation+'-receipts');receipts.mkdir(exist_ok=False)
    output=Path(str(BASE)+'-'+generation)
    previous=0.
    for path in BASE.parent.glob(BASE.name+'-*-receipts/exit.json'):
        previous+=json.loads(path.read_bytes())['elapsed_seconds']
    remaining=1800-previous
    if remaining<60:raise RuntimeError('cumulative Native grant exhausted')
    unit='p254-a50-joint-'+generation+'.service'
    manager=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus']
    contract=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())
    argv=manager+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec='+str(int(remaining-5)), '--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=9',
        'env',contract['common_runtime']['python_environment_assignment']]+[
        k+'=1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')]+[
        'taskset','-c','9',PYTHON,'-u',str(HERE/'native_wave.py'),str(output),mode]
    start=time.monotonic()
    save(receipts/'argv.json',dict(argv=argv,start_unix=time.time(),prior_native_seconds=previous,remaining_seconds=remaining))
    with (receipts/'stdout').open('xb') as out,(receipts/'stderr').open('xb') as err:
        proc=subprocess.Popen(argv,stdout=out,stderr=err)
        time.sleep(1)
        query=manager+['systemctl','--user','show',unit,'--property=Id,InvocationID,MainPID,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,ActiveState,SubState']
        startup=subprocess.run(query,capture_output=True,check=False)
        save(receipts/'manager-startup.json',dict(argv=query,exitcode=startup.returncode,stdout=startup.stdout.decode(),stderr=startup.stderr.decode()))
        code=proc.wait()
    elapsed=time.monotonic()-start
    save(receipts/'exit.json',dict(argv=argv,exitcode=code,elapsed_seconds=elapsed,cumulative_native_seconds=previous+elapsed))
    print((receipts/'stdout').read_text(),end='',flush=True)
    print((receipts/'stderr').read_text(),end='',file=sys.stderr,flush=True)
    print('A50_TERMINAL_'+generation+'_EXIT_'+str(code),flush=True)
    return code


if __name__=='__main__':sys.exit(main(sys.argv[1],sys.argv[2]))
