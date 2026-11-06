"""One-shot CPU8/3GiB native slot with all-exit transport receipts."""
from pathlib import Path
import json
import os
import resource
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'
ROOT = Path('/var/tmp/p254-a65-time-factor-001')
UNIT = 'p254-a65-time-factor-001.service'


def save(path, record):
    with Path(path).open('w') as stream:
        json.dump(record,stream,indent=2,sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def native():
    resource.setrlimit(resource.RLIMIT_CPU,(232,233))
    started = time.monotonic()
    code = 1
    receipt = dict(start_wall_time=time.time(),symbolic_preceding_wall_seconds=1.97)
    try:
        import material_pilot
        frozen = json.loads((HERE/'preregistration.json').read_bytes())
        for entry in frozen['source_inventory']:
            if material_pilot.binding(entry['path']) != entry:
                raise RuntimeError('frozen upstream source changed: '+entry['path'])
        cgroup_text = Path('/proc/self/cgroup').read_text()
        cg = Path('/sys/fs/cgroup') / cgroup_text.strip().split('::',1)[1].lstrip('/')
        environment = {k:v for k,v in os.environ.items() if k.startswith('HERDR_') or k in (
            'PYTHONPATH','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')}
        constraints = {name:(cg/name).read_text().strip() for name in ('memory.max','memory.swap.max')}
        cpuset = cg/'cpuset.cpus.effective'
        receipt.update(affinity=sorted(os.sched_getaffinity(0)),environment=environment,cgroup=cgroup_text,
            constraints=constraints,cpuset_effective=cpuset.read_text() if cpuset.exists() else None,
            native_sources=[material_pilot.binding(p) for p in HERE.glob('*.py')],
            interpreter=material_pilot.binding(PYTHON),view=material_pilot.binding(HERE/'native-view.json'),
            symbolic=material_pilot.binding(HERE/'symbolic-receipt.json'),
            preregistration=material_pilot.binding(HERE/'preregistration.json'))
        save(ROOT/'native-startup.json',receipt)
        if receipt['affinity'] != [8] or constraints != {'memory.max':'3221225472','memory.swap.max':'0'}:
            raise RuntimeError('actual bounded resource grant mismatch')
        for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
            if environment[key] != '1':
                raise RuntimeError('one-thread environment mismatch')
        if environment.get('HERDR_ENV') != '1' or not environment.get('HERDR_PANE_ID'):
            raise RuntimeError('actual Herdr pane context absent')
        # Reuse the immediately executed exact check; its transcribed receipt is bound.
        material_pilot.run(ROOT)
        receipt['result'] = material_pilot.binding(ROOT/'result.json')
        code = 0
    except BaseException as error:
        receipt['failure'] = repr(error)
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(exit=code,wall_seconds=time.monotonic()-started,
            cpu_seconds=usage.ru_utime+usage.ru_stime,maxrss_kib=usage.ru_maxrss)
        save(ROOT/'native-outcome.json',receipt)
        print('A65_NATIVE_DONE' if code == 0 else 'A65_NATIVE_FAIL',flush=True)


def launch():
    ROOT.mkdir(exist_ok=False)
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())
    runtime = contract['common_runtime']
    forwarded = {key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    if forwarded.get('HERDR_ENV') != '1' or not forwarded.get('HERDR_PANE_ID'):
        raise RuntimeError('launcher is not in the owned Herdr view')
    argv = runtime['manager_env']+['systemd-run','--user','--wait','--pipe','--collect',
        '--unit='+UNIT,'--property=RuntimeMaxSec=235','--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=StartLimitBurst=1',
        '--property=CPUAffinity=8','env',runtime['python_environment_assignment'],
        'PYTHONDONTWRITEBYTECODE=1']+[key+'='+value for key,value in runtime['threads'].items()]+[
        key+'='+value for key,value in forwarded.items()]+['taskset','-c','8',PYTHON,'-u',str(Path(__file__).resolve()),'--native']
    record = dict(argv=argv,unit=UNIT,herdr=forwarded,started=time.time(),automatic_retries=0,
                  cap_semantics='235s wall / 232s CPU hard process slot plus 1.97s preceding symbolic run; prior global costs retained')
    save(ROOT/'launch.json',record)
    code = 1
    try:
        with (ROOT/'stdout').open('wb') as stdout, (ROOT/'stderr').open('wb') as stderr:
            completed = subprocess.run(argv,stdout=stdout,stderr=stderr,timeout=245,check=False)
            code = completed.returncode
    except BaseException as error:
        record['transport_failure'] = repr(error)
        raise
    finally:
        record.update(exit=code,finished=time.time())
        save(ROOT/'terminal.json',record)
        print('A65_TERMINAL_EXIT_'+str(code),flush=True)
    if code:
        raise SystemExit(code)


if __name__ == '__main__':
    native() if sys.argv[1:] == ['--native'] else launch()
