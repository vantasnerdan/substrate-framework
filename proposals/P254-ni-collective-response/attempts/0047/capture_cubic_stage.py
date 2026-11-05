"""Capture source assembly or new-boundary consumer under remaining amended budget."""
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
    parser.add_argument('stage',choices=('assemble','feedback'))
    parser.add_argument('--charged-before',type=float,required=True)
    parser.add_argument('--runtime-max',type=int,required=True)
    parser.add_argument('--generation',type=int,default=1)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.charged_before<4123.160 or args.runtime_max<=0 or args.charged_before+args.runtime_max>9000:
        raise ValueError('actual cumulative Native budget amendment does not admit this cap')
    name='%s-%03d'%(args.stage,args.generation);unit='p254-a47-'+name;prefix=HERE/('actual-'+name)
    if args.stage=='assemble':
        script=HERE/'assemble_full37_source.py'
        arguments=['--output','/var/tmp/p254-a47-full37-source-001']+(['--resume'] if args.resume else [])
    else:
        if args.resume:raise ValueError('a new boundary consumer is an exclusive actual run')
        script=HERE/'consume_cubic_feedback.py'
        arguments=['--source','/var/tmp/p254-a47-full37-source-001/source.json','--output','/var/tmp/p254-a47-cubic-feedback-001']
    pythonpath=':'.join([str(ROOT/'src'),*[str(HERE.parent/name) for name in ('0039','0028','0027','0025','0024')],'/home/administrator/.cache/dft/p254-a28-runtime'])
    argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
        'systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec='+str(args.runtime_max),'--property=MemoryMax=4294967296','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=6','--property=WorkingDirectory='+str(ROOT),
        '/usr/bin/env','PYTHONPATH='+pythonpath,
        *[name+'=1' for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')],
        '/usr/bin/taskset','-c','6',PYTHON,'-u',str(HERE/'full37_kernel_entry.py'),
        str(prefix.with_suffix('.kernel-start.json')),unit,'6',str(script),*arguments]
    started=time.monotonic()
    startup=dict(schema='p254-a47-actual-native-start-v1',unit=unit,argv=argv,numerical_argv=[PYTHON,'-u',str(script),*arguments],
        start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),runtime_max_seconds=args.runtime_max,cpu=6,
        amended_whole_native_budget_seconds=9000,actual_cumulative_native_seconds_before_run=args.charged_before,
        source_bindings=[binding(path) for path in (script,Path(__file__),HERE/'full37_kernel_entry.py',HERE/'cubic_feedback.py',
            A43/'cubic_material_source.py',A43/'cubic_electronic_donor.py',A43/'cubic_position_moments.py',A43/'cubic_ionic_donor.py')])
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
