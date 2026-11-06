"""Capture each granted actual native invocation without closing its Herdr pane."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'
PYTHONPATH = ':'.join([str(ROOT/'src'), *[str(HERE.parent/x) for x in ('0039','0028','0027','0025','0024')], '/home/administrator/.cache/dft/p254-a28-runtime'])

def binding(path):
    raw = Path(path).read_bytes()
    return dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())

def publish(path, record):
    with path.open('w') as stream:
        json.dump(record,stream,indent=2,sort_keys=True); stream.write('\n')
        stream.flush(); os.fsync(stream.fileno())

def main():
    name, limit, command, *options = sys.argv[1:]
    unit = 'p254-a43-'+name
    prefix = HERE/('actual-'+name)
    numerical = [PYTHON,'-u',str(HERE/'run_material_cubic_kernel.py'),command,*options]
    if command == 'python-script':
        numerical = [PYTHON,'-u',*options]
    entry = [PYTHON,'-u',str(HERE/'actual_material_kernel_entry.py'),
             str(prefix.with_suffix('.kernel-start.json')),unit,*numerical[2:]]
    argv = ['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
        'systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec='+limit,'--property=MemoryMax=4294967296',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=6',
        '--property=WorkingDirectory='+str(ROOT),
        '/usr/bin/env','PYTHONPATH='+PYTHONPATH,
        *[x+'=1' for x in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')],
        '/usr/bin/taskset','-c','6',*entry]
    start = time.monotonic()
    startup = dict(schema='p254-a43-actual-native-start-v1',start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        argv=argv,numerical_argv=numerical,unit=unit,cwd=str(ROOT),cpu=6,memory_bytes=4294967296,
        swap_bytes=0,restart='no',runtime_max_seconds=int(limit),
        source_bindings=[binding(HERE/x) for x in ('cubic_position_moments.py','run_material_cubic_kernel.py','cubic_material_source.py','cubic_electronic_donor.py','cubic_ionic_donor.py','cubic_periodic_material.py','cubic_material_kernel.py','cubic_force_functional.py','cubic_cartesian_reference.py')])
    startup['kernel_entry'] = binding(HERE/'actual_material_kernel_entry.py')
    if command == 'python-script':
        startup['executed_script'] = binding(options[0])
    publish(prefix.with_suffix('.startup.json'),startup)
    with prefix.with_suffix('.stdout.log').open('xb') as out, prefix.with_suffix('.stderr.log').open('xb') as err:
        process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        def copy(source,destination,terminal):
            for line in iter(source.readline,b''):
                destination.write(line); destination.flush(); terminal.buffer.write(line); terminal.flush()
            os.fsync(destination.fileno())
        threads = [threading.Thread(target=copy,args=(process.stdout,out,sys.stdout)),threading.Thread(target=copy,args=(process.stderr,err,sys.stderr))]
        for thread in threads: thread.start()
        code = process.wait()
        for thread in threads: thread.join()
    receipt = dict(schema='p254-a43-actual-native-exit-v1',unit=unit,manager_exit_code=code,
        end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_seconds=time.monotonic()-start,
        startup=binding(prefix.with_suffix('.startup.json')),stdout=binding(prefix.with_suffix('.stdout.log')),stderr=binding(prefix.with_suffix('.stderr.log')))
    publish(prefix.with_suffix('.exit.json'),receipt)
    print('[p254-a43-'+name+'-native-exit] code='+str(code),flush=True)

if __name__ == '__main__':
    main()
