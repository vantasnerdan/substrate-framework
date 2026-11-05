"""Owned A75 CPU10 native production and actual saved-path consumer, with all-exit costs."""
from pathlib import Path
import json
import os
import resource
import subprocess
import sys
import time
import traceback
HERE = Path(__file__).resolve().parent


def dump(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def native():
    started = time.monotonic()
    result = dict(schema='p254-a75-native-findings-v1',status='running')
    try:
        view = json.loads((HERE/'owned-view.json').read_bytes())
        context = {key:value for key,value in os.environ.items() if key.startswith('HERDR_') or key in
            ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')}
        cg = Path('/sys/fs/cgroup')/Path('/proc/self/cgroup').read_text().strip().split('::')[-1].lstrip('/')
        runtime = dict(affinity=sorted(os.sched_getaffinity(0)),environment=context,
            cgroup={key:(cg/key).read_text().strip() for key in ('memory.max','memory.swap.max')},
            python=sys.executable,pythonpath=os.environ.get('PYTHONPATH'))
        result['runtime'] = runtime
        if runtime['affinity'] != [10] or context.get('HERDR_PANE_ID') != view['pane_id'] or context.get('HERDR_ENV') != '1':
            raise RuntimeError('owned HERDR/CPU10 mismatch')
        if runtime['cgroup'] != {'memory.max':'3221225472','memory.swap.max':'0'}:
            raise RuntimeError('resource constraints differ')
        if any(context.get(k) != '1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')):
            raise RuntimeError('one-thread contract mismatch')
        dump(HERE/'native-startup.json',runtime)
        from produce_joint import produce
        from component_reader import consume
        catalog = produce(HERE/'component-001',started+1750)
        record = json.loads(catalog.read_bytes())
        exercise = consume(catalog,require_complete=record['status']=='complete-global-joint-component')
        dump(HERE/'saved-path-exercise.json',exercise)
        from history_banks import binding
        result.update(status=record['status'],catalog=binding(catalog),saved_path_exercise=binding(HERE/'saved-path-exercise.json'),
            completed_pairs=exercise['completed_pairs'],failure=record['failure'],scope=record['scope'])
        if record['status'] != 'complete-global-joint-component':
            raise RuntimeError('global joint component incomplete; saved actual partial boundary')
    except BaseException as exc:
        result.update(status='failed',failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        raise
    finally:
        result.update(native_elapsed_seconds=time.monotonic()-started,native_cpu_seconds=time.process_time(),
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
        dump(HERE/'native-findings.json',result)
        print('A75_NATIVE_TERMINAL '+json.dumps(result),flush=True)


def launch():
    receipts = HERE/'receipts-001'
    receipts.mkdir(exist_ok=False)
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    context = {key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    view = json.loads((HERE/'owned-view.json').read_bytes())
    if context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != view['pane_id']:
        raise RuntimeError('launch must originate in actual owned pane')
    unit = 'p254-a75-global-joint-full31-001.service'
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec=1795','--property=TimeoutStopSec=5',
        '--property=MemoryMax=3221225472','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=10','env',contract['python_environment_assignment'],
        *[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in context.items()],
        'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),
        'taskset','-c','10',contract['python'],'-u',str(Path(__file__).resolve()),'--native']
    started, wall = time.monotonic(), time.time()
    dump(receipts/'launch.json',dict(argv=argv,start_unix=wall,unit=unit,context=context,owned_view=view['pane_id'],
        grant_seconds=1800,systemd_cap_seconds=1795,production_boundary_seconds=1750,
        reserve='45s finalization/source hashing and actual saved-path consumer',
        checkpoint_granularity='one complete original128-sector/114-coordinate pair; all six numeric NPY slices fsync before coverage',
        expected_useful_output='961 actual global jointGamma correction Sigma and full jointGamma Pi pairs once'))
    rc = None
    try:
        with (receipts/'stdout.txt').open('xb') as stdout,(receipts/'stderr.txt').open('xb') as stderr:
            proc = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in proc.stdout:
                stdout.write(line)
                stdout.flush()
                sys.stdout.buffer.write(line)
                sys.stdout.buffer.flush()
            rc = proc.wait()
        print((receipts/'stderr.txt').read_text(),flush=True)
    finally:
        query = contract['manager_env']+['systemctl','--user','show',unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,CPUAffinity,ActiveState,SubState']
        manager = subprocess.run(query,capture_output=True,check=False)
        dump(receipts/'manager-terminal.json',dict(argv=query,exitcode=manager.returncode,stdout=manager.stdout.decode(),stderr=manager.stderr.decode()))
        dump(receipts/'exit.json',dict(argv=argv,exitcode=rc,start_unix=wall,end_unix=time.time(),elapsed_seconds=time.monotonic()-started))
        print('A75_OWNED_LAUNCH_TERMINAL rc='+str(rc),flush=True)
    return rc


if __name__ == '__main__':
    if '--native' in sys.argv:
        native()
    else:
        sys.exit(launch())
