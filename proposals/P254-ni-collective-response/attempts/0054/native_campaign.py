"""Owned CPU9 systemd-user runner with cumulative cost and exact raw receipts."""
from pathlib import Path
import json
import os
import resource
import sys
import time
from supplied_keldysh import A49
from command_receipts import binding, run_recorded

HERE = Path(__file__).resolve().parent
PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'
PREFIX = Path('/home/administrator/.local/opt/nessi-6bb4c832')


def native(root):
    start = time.monotonic()
    record = {'affinity':sorted(os.sched_getaffinity(0)),
              'environment':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
              'sources':[binding(p) for p in sorted(HERE.iterdir()) if p.is_file()],
              'upstream':[binding(p) for p in [A49/'supplied_vie2.hpp',A49/'component_io.hpp',
                  A49/'supplied_vie2.py',A49/'command_receipts.py',A49/'finite_gaussian_reference.py',
                  HERE.parent/'0040'/'nessi_gaussian_ivp.hpp',PREFIX/'lib'/'libcntr.so',
                  PREFIX/'include'/'cntr'/'cntr_vie2_impl.hpp',PREFIX/'include'/'cntr'/'cntr_convolution_impl.hpp',
                  PREFIX/'include'/'cntr'/'integration.hpp']]}
    (root/'native-bindings.json').write_text(json.dumps(record,indent=2)+'\n')
    try:
        if record['affinity'] != [9] or any(v!='1' for v in record['environment'].values()):
            raise RuntimeError('actual CPU/thread grant mismatch')
        prior = root.parent/'p254-a54-keldysh-001'
        if root!=prior and (prior/'symbolic-block.json').exists():
            previous = json.loads((prior/'native-bindings.json').read_text())
            old = next(p for p in previous['sources'] if p['path']==str(HERE/'symbolic_boundary.py'))
            if old!=binding(HERE/'symbolic_boundary.py'):
                raise RuntimeError('symbolic block changed; explicit new execution required')
            record['symbolic_carried'] = binding(prior/'symbolic-block.json')
        else:
            run_recorded([PYTHON,str(HERE/'symbolic_boundary.py'),str(root/'symbolic-block.json')],root,'symbolic',60)
        run_recorded(['cmake','-S',str(HERE),'-B',str(root/'build'),'-DCMAKE_BUILD_TYPE=Release'],root,'configure',60)
        run_recorded(['cmake','--build',str(root/'build'),'--parallel','1'],root,'build',400)
        binary = root/'build'/'solve_components'
        run_recorded([PYTHON,str(HERE/'finite_boundary_smoke.py'),str(binary),str(root)],root,'changed-oracle-and-actual',120)
        record['exit'] = 0
    except Exception as exc:
        record.update(exit=1,failure=repr(exc))
        raise
    finally:
        children,own = resource.getrusage(resource.RUSAGE_CHILDREN),resource.getrusage(resource.RUSAGE_SELF)
        record.update(wall_seconds=time.monotonic()-start,
                      cpu_seconds=children.ru_utime+children.ru_stime+own.ru_utime+own.ru_stime,
                      child_maxrss_kib=children.ru_maxrss,owner_maxrss_kib=own.ru_maxrss,
                      output_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        (root/'native-outcome.json').write_text(json.dumps(record,indent=2)+'\n')


def launch(root):
    root.mkdir(parents=True,exist_ok=False)
    prior = [json.loads(p.read_text()) for p in root.parent.glob('p254-a54-keldysh-*/systemd-native.json') if p.parent!=root]
    charged = sum(p['wall_seconds'] for p in prior)
    if charged+600>1800:
        raise RuntimeError('remaining cumulative Native allocation below attempt ceiling')
    (root/'cumulative-grant.json').write_text(json.dumps({'charged_wall_before':charged,'cumulative_limit':1800,'attempt_limit':600},indent=2)+'\n')
    argv = ['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
            'systemd-run','--user','--wait','--pipe','--collect',
            '--unit=p254-a54-keldysh-'+root.name.rsplit('-',1)[-1]+'.service',
            '--property=RuntimeMaxSec=600','--property=MemoryMax=3221225472',
            '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=9','env',
            'OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1',
            'BLIS_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1','taskset','-c','9',
            PYTHON,'-u',str(Path(__file__).resolve()),'--native',str(root)]
    try:
        run_recorded(argv,root,'systemd-native',610)
        print('A54_NATIVE_TERMINAL_'+root.name.rsplit('-',1)[-1]+'_EXIT_0',flush=True)
    except Exception:
        print('A54_NATIVE_TERMINAL_'+root.name.rsplit('-',1)[-1]+'_EXIT_NONZERO',flush=True)
        raise

if __name__ == '__main__':
    if sys.argv[1]=='--native':
        native(Path(sys.argv[2]))
    else:
        launch(Path(sys.argv[1]))
