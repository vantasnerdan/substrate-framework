"""A69 bounded CPU9 single-thread original projected control and saved-path exercise."""
from pathlib import Path
import json
import os
import resource
import subprocess
import sys
import time
import traceback
HERE=Path(__file__).resolve().parent


def dump(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def native():
    started=time.monotonic()
    result=dict(schema='p254-a69-native-findings-v1',status='running',findings=[],blockers=[])
    try:
        from projected_gamma_control import produce,consume,jsonable,save_json
        result['runtime']=dict(affinity=sorted(os.sched_getaffinity(0)),
            environment={key:value for key,value in os.environ.items() if key.startswith('HERDR_') or key in (
                'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')})
        cg=Path('/sys/fs/cgroup')/Path('/proc/self/cgroup').read_text().strip().split('::')[-1].lstrip('/')
        result['runtime']['cgroup']={key:(cg/key).read_text().strip() for key in ('memory.max','memory.swap.max')}
        if result['runtime']['affinity'] != [9] or os.environ.get('HERDR_PANE_ID') != 'w1:p78' or os.environ.get('HERDR_ENV') != '1':
            raise RuntimeError('verified original owned view/CPU9 required')
        if result['runtime']['cgroup'] != {'memory.max':'3221225472','memory.swap.max':'0'}:
            raise RuntimeError('assigned cgroup resource constraints differ')
        if any(os.environ.get(k) != '1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')):
            raise RuntimeError('single-thread constraint differs')
        dump(HERE/'native-startup.json',result['runtime'])
        path=produce(HERE/'control-001')
        exercise=consume(path)
        save_json(HERE/'saved-path-exercise.json',exercise)
        from material_electron_ion_inputs import file_binding
        result['control']=file_binding(path)
        result['saved_path_exercise']=file_binding(HERE/'saved-path-exercise.json')
        result['findings']=[dict(status='complete',scope='original represented full coupled Gaussian one-probe covariance at genuine0,30'),
            dict(status='exercised',scope='actual saved NPZ read through canonical normalized_ordered, both raw cross orders and variances retained')]
        result['status']='complete'
    except BaseException as exc:
        result['status']='failed'
        result['blockers'].append(dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        raise
    finally:
        result['native_elapsed_seconds']=time.monotonic()-started
        result['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        dump(HERE/'native-findings.json',result)
        print('A69_NATIVE_TERMINAL '+result['status'],flush=True)


def launch():
    receipts=HERE/'receipts-001';receipts.mkdir(exist_ok=False)
    contract=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    context={key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    if os.environ.get('HERDR_ENV') != '1' or os.environ.get('HERDR_PANE_ID') != 'w1:p78' or not all(context.get(key) for key in ('HERDR_TAB_ID','HERDR_WORKSPACE_ID')):
        raise RuntimeError('actual owned view context missing; no synthetic Herdr environment')
    unit='p254-a69-projected-gamma-control-001.service'
    argv=contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec=295','--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=9',
        'env',contract['python_environment_assignment'],*[k+'='+v for k,v in contract['threads'].items()],
        *[k+'='+v for k,v in context.items()], 'PYTHONDONTWRITEBYTECODE=1',
        '/usr/bin/time','-v','-o',str(receipts/'time.txt'),'taskset','-c','9',
        contract['python'],'-u',str(Path(__file__).resolve()),'--native']
    started=time.monotonic();wall=time.time()
    dump(receipts/'launch.json',dict(argv=argv,start_unix=wall,context=context,unit=unit,
        grant_seconds=300,systemd_cap_seconds=295,owned_view='w1:p78',label='A69 original projected Gamma coupled control',
        expected_useful_output='one full coupled one-probe adjoint0<-30 and exact2x2 raw covariance',
        prior_costs='carried unchanged from original preregistered campaign; no replay'))
    rc=None
    try:
        with (receipts/'stdout.txt').open('xb') as stdout,(receipts/'stderr.txt').open('xb') as stderr:
            rc=subprocess.run(argv,stdout=stdout,stderr=stderr,check=False).returncode
    finally:
        query=contract['manager_env']+['systemctl','--user','show',unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        manager=subprocess.run(query,capture_output=True,check=False)
        dump(receipts/'manager-terminal.json',dict(argv=query,exitcode=manager.returncode,stdout=manager.stdout.decode(),stderr=manager.stderr.decode()))
        dump(receipts/'exit.json',dict(argv=argv,exitcode=rc,start_unix=wall,end_unix=time.time(),elapsed_seconds=time.monotonic()-started))
        print('A69_OWNED_LAUNCH_TERMINAL rc='+str(rc),flush=True)
    return rc


if __name__=='__main__':
    if '--native' in sys.argv:native()
    else:sys.exit(launch())
