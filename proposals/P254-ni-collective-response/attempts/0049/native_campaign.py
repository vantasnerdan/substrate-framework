"""Owned systemd-user Native campaign, durable receipts and terminal completion."""
from pathlib import Path
import json
import os
import resource
import sys
import time
from command_receipts import binding, run_recorded

HERE = Path(__file__).resolve().parent
PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'
PREFIX = Path('/home/administrator/.local/opt/nessi-6bb4c832')


def native(root):
    start = time.monotonic()
    record = {'affinity':sorted(os.sched_getaffinity(0)),
              'environment':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
              'scope':'new A49 symbolic/build/affected-boundary finite oracle only; original A40 evidence carried',
              'sources':[binding(p) for p in sorted(HERE.iterdir()) if p.is_file()],
              'upstream':[binding(HERE.parent/'0040'/'nessi_gaussian_ivp.hpp'),
                          binding(PREFIX/'lib'/'libcntr.so'),
                          binding(PREFIX/'include'/'cntr'/'cntr_vie2_impl.hpp'),
                          binding(PREFIX/'include'/'cntr'/'cntr_convolution_impl.hpp')]}
    (root/'native-bindings.json').write_text(json.dumps(record,indent=2)+'\n')
    try:
        if record['affinity'] != [8] or any(v != '1' for v in record['environment'].values()):
            raise RuntimeError('actual CPU/thread grant mismatch')
        prior = root.parent/'p254-a49-supplied-vie2-001'
        if root != prior and (prior/'symbolic-blocks.json').exists():
            previous = json.loads((prior/'native-bindings.json').read_text())
            old = next(p for p in previous['sources'] if p['path'] == str(HERE/'symbolic_boundary.py'))
            if old['sha256'] != binding(HERE/'symbolic_boundary.py')['sha256']:
                raise RuntimeError('symbolic source changed; explicit new block required')
            record['symbolic_carried'] = binding(prior/'symbolic-blocks.json')
        else:
            run_recorded([PYTHON,str(HERE/'symbolic_boundary.py'),str(root/'symbolic-blocks.json')],root,'symbolic',60)
        built = root.parent/'p254-a49-supplied-vie2-002'
        binary = built/'build'/'solve_components'
        if root != built and binary.exists():
            previous = json.loads((built/'native-bindings.json').read_text())
            for name in ('CMakeLists.txt','supplied_vie2.hpp','component_io.hpp','solve_components.cpp'):
                old = next(p for p in previous['sources'] if p['path'] == str(HERE/name))
                if old['sha256'] != binding(HERE/name)['sha256']:
                    raise RuntimeError('native source changed; explicit rebuild required')
            if previous['upstream'] != record['upstream']:
                raise RuntimeError('upstream bindings changed; explicit rebuild required')
            record['build_carried'] = binding(binary)
        else:
            run_recorded(['cmake','-S',str(HERE),'-B',str(root/'build'),'-DCMAKE_BUILD_TYPE=Release'],root,'configure',60)
            run_recorded(['cmake','--build',str(root/'build'),'--parallel','1'],root,'build',400)
            binary = root/'build'/'solve_components'
        run_recorded([PYTHON,str(HERE/'finite_boundary_smoke.py'),str(binary),str(root)],root,'smoke',120)
        record['exit'] = 0
    except Exception as exc:
        record['exit'] = 1
        record['failure'] = repr(exc)
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        own = resource.getrusage(resource.RUSAGE_SELF)
        record.update(wall_seconds=time.monotonic()-start,
                      cpu_seconds=usage.ru_utime+usage.ru_stime+own.ru_utime+own.ru_stime,
                      child_maxrss_kib=usage.ru_maxrss,owner_maxrss_kib=own.ru_maxrss)
        (root/'native-outcome.json').write_text(json.dumps(record,indent=2)+'\n')


def launch(root):
    root.mkdir(parents=True,exist_ok=False)
    prior = [json.loads(p.read_text()) for p in root.parent.glob('p254-a49-supplied-vie2-*/systemd-native.json') if p.parent != root]
    charged = sum(p['wall_seconds'] for p in prior)
    if charged+600 > 1800:
        raise RuntimeError('remaining cumulative Native grant below per-attempt ceiling')
    (root/'cumulative-grant.json').write_text(json.dumps({'charged_wall_before':charged,'cumulative_limit':1800,'attempt_limit':600},indent=2)+'\n')
    argv = ['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
            'systemd-run','--user','--wait','--pipe','--collect',
            '--unit=p254-a49-supplied-vie2-'+root.name.rsplit('-',1)[-1]+'.service',
            '--property=RuntimeMaxSec=600','--property=MemoryMax=3221225472',
            '--property=MemorySwapMax=0','--property=Restart=no','env',
            'OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',
            'NUMEXPR_NUM_THREADS=1','BLIS_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1',
            'taskset','-c','8',PYTHON,'-u',str(Path(__file__).resolve()),'--native',str(root)]
    try:
        run_recorded(argv,root,'systemd-native',610)
        print('A49_'+'NATIVE_TERMINAL_'+root.name.rsplit('-',1)[-1]+'_EXIT_0',flush=True)
    except Exception:
        print('A49_'+'NATIVE_TERMINAL_'+root.name.rsplit('-',1)[-1]+'_EXIT_NONZERO',flush=True)
        raise

if __name__ == '__main__':
    if sys.argv[1] == '--native':
        native(Path(sys.argv[2]))
    else:
        launch(Path(sys.argv[1]))
