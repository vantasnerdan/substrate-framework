"""Owned CPU9 production transport build/exercise; every exit durable and visible."""
from pathlib import Path
import json
import os
import resource
import sys
import time
from production_transport import binding, run_recorded, atomic_json

HERE=Path(__file__).resolve().parent
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
PREFIX=Path('/home/administrator/.local/opt/nessi-6bb4c832')
ROOT=Path('/var/tmp/p254-a61-lean-001')


def native(root):
    started=time.monotonic()
    receipt=dict(affinity=sorted(os.sched_getaffinity(0)),herdr_env=os.environ.get('HERDR_ENV'),environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},sources=[binding(p) for p in sorted(HERE.iterdir()) if p.is_file()],upstream=[binding(p) for p in (HERE.parent/'0054'/'supplied_keldysh.hpp',HERE.parent/'0049'/'supplied_vie2.hpp',HERE.parent/'0049'/'component_io.hpp',HERE.parent/'0049'/'command_receipts.py',HERE.parent/'0040'/'nessi_gaussian_ivp.hpp',PREFIX/'lib'/'libcntr.so',PREFIX/'include'/'cntr'/'cntr_vie2_impl.hpp',PREFIX/'include'/'cntr'/'cntr_convolution_impl.hpp',PREFIX/'include'/'cntr'/'integration.hpp')],expected_material=[binding(Path('/var/tmp/p254-a58-common-feedback-step-002')/p) for p in ('result.json','actual-material-step.npz','fermions/slot-000/manifest.json','intrinsic/manifest.json')])
    atomic_json(root/'native-bindings.json',receipt)
    exitcode=1
    try:
        if receipt['affinity'] != [9] or receipt['herdr_env'] != '1' or any(v!='1' for v in receipt['environment'].values()):
            raise RuntimeError('actual grant/environment mismatch')
        run_recorded(['cmake','-S',str(HERE),'-B',str(root/'build'),'-DCMAKE_BUILD_TYPE=Release'],root,'configure',30)
        run_recorded(['/usr/bin/time','-v','-o',str(root/'build.time'),'cmake','--build',str(root/'build'),'--parallel','1'],root,'build',350)
        run_recorded([PYTHON,'-u',str(HERE/'exercise_transport.py'),str(root/'build'/'solve_result'),str(root)],root,'exercise',220)
        receipt['exit']=exitcode=0
    except BaseException as exc:
        receipt.update(exit=1,failure=repr(exc))
        raise
    finally:
        own=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(wall_seconds=time.monotonic()-started,cpu_seconds=own.ru_utime+own.ru_stime+children.ru_utime+children.ru_stime,owner_maxrss_kib=own.ru_maxrss,child_maxrss_kib=children.ru_maxrss,output_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        atomic_json(root/'native-outcome.json',receipt)
        print('[dft-done]' if exitcode==0 else '[dft-fail]', 'A61 native exit='+str(exitcode),flush=True)


def launch(root):
    if root != ROOT:
        raise ValueError('undeclared staging root')
    root.mkdir(exist_ok=False)
    argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus','systemd-run','--user','--wait','--pipe','--collect','--unit=p254-a61-lean-001.service','--property=RuntimeMaxSec=590','--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=StartLimitBurst=1','--property=CPUAffinity=9','env','HERDR_ENV=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','BLIS_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1','taskset','-c','9',PYTHON,'-u',str(Path(__file__).resolve()),'--native',str(root)]
    code=1
    try:
        run_recorded(argv,root,'systemd-native',595)
        code=0
    finally:
        atomic_json(root/'terminal.json',dict(exit=code,root=str(root),marker='A61_NATIVE_TERMINAL_001_EXIT_'+str(code)))
        print('A61_NATIVE_TERMINAL_001_EXIT_'+str(code),flush=True)


if __name__=='__main__':
    if sys.argv[1]=='--native':native(Path(sys.argv[2]))
    else:launch(Path(sys.argv[1]))
