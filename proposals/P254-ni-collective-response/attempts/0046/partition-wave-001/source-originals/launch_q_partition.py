"""Stream one authorized partition Native through its exclusively owned live view."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
import sys
import threading
import time
HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name',choices=('A','B','C'),required=True)
    p.add_argument('--q-start',type=int,required=True); p.add_argument('--q-stop',type=int,required=True)
    p.add_argument('--cpu',type=int,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    if (out/'actual-native-argv.json').exists():
        raise FileExistsError('original actual invocation already retained; do not overwrite/relaunch')
    contract=json.loads((HERE/'preregistration-original.json').read_bytes())
    rt=contract['common_runtime']
    unit='p254-a46-partition-'+a.name.lower()+'-001.service'
    argv=rt['manager_env']+['systemd-run','--user','--wait','--pipe','--collect','--unit='+unit,
        '--property=RuntimeMaxSec=9000','--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity='+str(a.cpu),
        'env',rt['python_environment_assignment'],'PYTHONDONTWRITEBYTECODE=1']
    argv += [k+'='+v for k,v in rt['threads'].items()]
    argv += ['taskset','-c',str(a.cpu),rt['python'],'-u',str(HERE/'native_q_partition.py'),
        '--name',a.name,'--q-start',str(a.q_start),'--q-stop',str(a.q_stop),
        '--cpu',str(a.cpu),'--output',str(out)]
    (out/'actual-native-argv.json').write_text(json.dumps(dict(argv=argv,
        start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),unit=unit),indent=2))
    start=time.monotonic()
    proc=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,bufsize=1)
    def stream(pipe,path,target):
        with path.open('w',buffering=1) as f:
            for line in pipe:
                f.write(line); target.write(line); target.flush()
    threads=[threading.Thread(target=stream,args=(proc.stdout,out/'manager.rawstdout.txt',sys.stdout)),
        threading.Thread(target=stream,args=(proc.stderr,out/'manager.rawstderr.txt',sys.stderr))]
    for t in threads: t.start()
    code=proc.wait()
    for t in threads: t.join()
    record=dict(exit=code,seconds=time.monotonic()-start,argv=argv,
        end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (out/'manager-exit.json').write_text(json.dumps(record,indent=2))
    print('A46_PARTITION_'+a.name+'_NATIVE_TERMINAL',code,flush=True)
    return code


if __name__=='__main__': sys.exit(main())
