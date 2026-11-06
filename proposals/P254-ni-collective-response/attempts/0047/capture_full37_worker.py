"""Original all-exit capture for an explicitly owned CPU6/CPU7 A47 worker."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

HERE=Path(__file__).resolve().parent;A43=HERE.parent/'0043';ROOT=HERE.parents[3]
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'

def binding(path):
    path=Path(path);raw=path.read_bytes()
    return dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())

def publish(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,sort_keys=True);stream.write('\n');stream.flush();os.fsync(stream.fileno())

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu',type=int,choices=(6,7),required=True)
    parser.add_argument('--generation',type=int,default=1)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    name='full37-cpu%d-%03d'%(args.cpu,args.generation)
    unit='p254-a47-'+name;prefix=HERE/('actual-'+name)
    script=HERE/'complete_full37_sectors.py'
    output='/var/tmp/p254-a47-full37-cpu%d-001'%args.cpu
    numerical=[str(script),'--cpu',str(args.cpu),'--output',output]+(['--resume'] if args.resume else [])
    pythonpath=':'.join([str(ROOT/'src'),*[str(HERE.parent/name) for name in ('0039','0028','0027','0025','0024')],'/home/administrator/.cache/dft/p254-a28-runtime'])
    argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
        'systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec=2400','--property=MemoryMax=4294967296','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity='+str(args.cpu),'--property=WorkingDirectory='+str(ROOT),
        '/usr/bin/env','PYTHONPATH='+pythonpath,
        *[name+'=1' for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')],
        '/usr/bin/taskset','-c',str(args.cpu),PYTHON,'-u',str(HERE/'full37_kernel_entry.py'),
        str(prefix.with_suffix('.kernel-start.json')),unit,str(args.cpu),*numerical]
    started=time.monotonic()
    startup=dict(schema='p254-a47-actual-native-start-v1',unit=unit,argv=argv,numerical_argv=[PYTHON,'-u',*numerical],
        start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),runtime_max_seconds=2400,cpu=args.cpu,
        owned_sector_indices=list(range(20,74)) if args.cpu==6 else list(range(74,128)),
        amended_whole_native_budget_seconds=9000,original_native_seconds_already_spent=4123.160,
        remaining_shared_native_seconds_before_two_workers=4876.840,
        source_bindings=[binding(path) for path in (script,Path(__file__),HERE/'full37_kernel_entry.py',A43/'cubic_material_source.py',A43/'cubic_electronic_donor.py',A43/'cubic_position_moments.py')],
        original_A43_source=binding('/var/tmp/p254-a43-full-state-third-002/source.json'),output=output)
    publish(prefix.with_suffix('.startup.json'),startup)
    with prefix.with_suffix('.stdout.log').open('xb') as stdout,prefix.with_suffix('.stderr.log').open('xb') as stderr:
        process=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        def copy(source,destination,terminal):
            for line in iter(source.readline,b''):
                destination.write(line);destination.flush();terminal.buffer.write(line);terminal.flush()
            os.fsync(destination.fileno())
        threads=[threading.Thread(target=copy,args=(process.stdout,stdout,sys.stdout)),threading.Thread(target=copy,args=(process.stderr,stderr,sys.stderr))]
        for thread in threads:thread.start()
        code=process.wait()
        for thread in threads:thread.join()
    publish(prefix.with_suffix('.exit.json'),dict(schema='p254-a47-actual-native-exit-v1',unit=unit,
        manager_exit_code=code,actual_wall_seconds=time.monotonic()-started,
        end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),startup=binding(prefix.with_suffix('.startup.json')),
        stdout=binding(prefix.with_suffix('.stdout.log')),stderr=binding(prefix.with_suffix('.stderr.log'))))
    print('[p254-a47-'+name+'-native-exit] code='+str(code),flush=True)

if __name__=='__main__':
    main()
