"""One authorized q0 nine-cell production/consumer, never a mesh launcher."""
from pathlib import Path
import json
import os
import resource
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent


def save(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def native():
    from component_emitter import produce_component, binding
    from component_consumer import consume_component
    began = time.monotonic()
    group = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
        limits={key:(cg/key).read_text().strip() for key in ('memory.max','memory.swap.max')},
        environment={key:value for key,value in os.environ.items() if key.startswith('HERDR_') or key in (
            'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
        executable=sys.executable,pythonpath=os.environ.get('PYTHONPATH'))
    save(HERE/'native-startup.json',runtime)
    view = json.loads((HERE/'owned-view.json').read_bytes())
    if (runtime['affinity'] != [8] or runtime['limits'] != {'memory.max':'2147483648','memory.swap.max':'0'}
            or runtime['environment'].get('HERDR_PANE_ID') != view['pane_id']
            or runtime['environment'].get('HERDR_ENV') != '1'
            or any(runtime['environment'].get(k) != '1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'))):
        raise RuntimeError('actual assigned Native runtime differs')
    result = dict(schema='p254-a76-native-findings-v1',status='running',runtime=runtime,components=[])
    try:
        for component in ('ee','mixed_ei','cubic'):
            tick = time.monotonic()
            catalog = produce_component(component,0,HERE/('q0-'+component+'-001'),times=[0.,15.,30.])
            consumer = consume_component(catalog)
            consumer_path = catalog.parent/'consumer.json'
            save(consumer_path,consumer)
            result['components'].append(dict(component=component,catalog=binding(catalog),consumer=binding(consumer_path),
                completed_pairs=consumer['completed_pairs'],scalar_elements_per_order=consumer['scalar_elements_consumed_per_order'],
                elapsed_seconds=time.monotonic()-tick))
            save(HERE/'native-findings.json',result)
            print('A76_COMPONENT_COMPLETE '+json.dumps(result['components'][-1]),flush=True)
        result['status'] = 'complete-owned-q0-nine-pairs-all-three-components'
    except BaseException as exc:
        result['status'] = 'failed'
        result['failure'] = dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
        raise
    finally:
        result.update(elapsed_seconds=time.monotonic()-began,native_cpu_seconds=time.process_time(),
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            full_mesh=False,full31_claim=False)
        save(HERE/'native-findings.json',result)
        print('A76_NATIVE_TERMINAL '+json.dumps(result),flush=True)


def launch():
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    context = {key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    view = json.loads((HERE/'owned-view.json').read_bytes())
    if context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != view['pane_id']:
        raise RuntimeError('launch must originate in actual owned pane')
    receipts = HERE/'receipts-001'
    receipts.mkdir(exist_ok=False)
    unit = 'p254-a76-q0-three-components-001.service'
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec=590','--property=MemoryMax=2147483648','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=8','env',contract['python_environment_assignment'],
        *[k+'='+v for k,v in contract['threads'].items()],*[k+'='+v for k,v in context.items()],
        'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(receipts/'time.txt'),
        contract['python'],'-u',str(HERE/'run_native.py'),'--native']
    began = time.monotonic()
    save(receipts/'launch.json',dict(argv=argv,unit=unit,actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        context=context,native_cap_seconds=600,scope='q0 all three, genuine 0/15/30 all nine pairs, no retries'))
    rc = None
    try:
        with (receipts/'stdout.txt').open('wb') as stdout, (receipts/'stderr.txt').open('wb') as stderr:
            proc = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in iter(proc.stdout.readline,b''):
                stdout.write(line)
                stdout.flush()
                sys.stdout.buffer.write(line)
                sys.stdout.buffer.flush()
            rc = proc.wait()
    finally:
        save(receipts/'exit.json',dict(argv=argv,returncode=rc,elapsed_seconds=time.monotonic()-began))
        shown = subprocess.run(contract['manager_env']+['systemctl','--user','show',unit,
            '--property=Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,ActiveState,SubState'],capture_output=True,text=True)
        save(receipts/'manager-terminal.json',dict(argv=shown.args,returncode=shown.returncode,stdout=shown.stdout,stderr=shown.stderr))
        print('A76_OWNED_LAUNCH_TERMINAL rc='+str(rc),flush=True)
    if rc:
        raise SystemExit(rc)


if __name__ == '__main__':
    native() if sys.argv[1:] == ['--native'] else launch()
