"""Launch durable owner, inspect the actual waiting worker, release and observe."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
import traceback
from receipts import save,observe
from source_graph import HERE,ROOT,load,digest,own_sources

PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
ENV={**os.environ,'XDG_RUNTIME_DIR':'/run/user/1000','DBUS_SESSION_BUS_ADDRESS':'unix:path=/run/user/1000/bus'}


def command(argv):
    result=subprocess.run(argv,env=ENV,text=True,capture_output=True)
    return {'argv':argv,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr}


def launch(generation):
    directory=HERE/('generation-'+generation);manifest=load(directory/'source-manifest.json');unit=manifest['unit']
    if own_sources()!=manifest['oracle_sources'] or (directory/'owner-launch.json').exists():
        raise ValueError('Changed source or existing launch; create a distinct generation')
    args=['systemd-run','--user','--unit='+unit,'--property=Type=exec','--property=CPUAffinity=5',
          '--property=MemoryMax=2147483648','--property=MemorySwapMax=0','--property=Restart=no',
          '--property=RemainAfterExit=yes','--property=KillMode=control-group','--working-directory='+str(ROOT)]
    for key,value in {'PYTHONPATH':str(ROOT/'src'),'PYTHONDONTWRITEBYTECODE':'1','PYTHONUNBUFFERED':'1',
                      'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'}.items():
        args+=['--setenv='+key+'='+value]
    args += [PYTHON,'-B',str(HERE/'run_oracle.py'),'--generation',str(directory)]
    result=command(args)
    save(directory/'owner-launch.json',{'launch':result,'unit':unit,'manifest_sha256':digest(directory/'source-manifest.json'),
         'pane_id':os.environ.get('HERDR_PANE_ID'),'tab_id':os.environ.get('HERDR_TAB_ID')})
    print('SENSITIVITY_OWNER_LAUNCH',unit,json.dumps(result),flush=True)
    try:
        if result['returncode']:
            raise RuntimeError('systemd launch failed')
        deadline=time.monotonic()+60
        while not (directory/'worker-start.json').exists():
            if time.monotonic()>deadline:
                raise TimeoutError('Worker never provided startup handshake')
            time.sleep(.1)
        worker=load(directory/'worker-start.json')['runtime']
        observed=observe(worker['pid'])
        mainpid=command(['systemctl','--user','show',unit,'--property=MainPID','--value'])
        save(directory/'worker-postexec-observed.json',{'worker':observed,'systemd_mainpid':mainpid,
             'worker_start_receipt_sha256':digest(directory/'worker-start.json')})
        keys=('pid','start_ticks','executable','executable_sha256','affinity','cgroup','cmdline')
        if any(observed[k]!=worker[k] for k in keys) or mainpid['stdout'].strip()!=str(worker['pid']):
            raise ValueError('External postexec identity mismatch')
        if observed['affinity']!=[5] or observed['live_kernel']['memory.max']!='2147483648' or observed['live_kernel']['memory.swap.max']!='0':
            raise ValueError('External live kernel resource mismatch')
        save(directory/'worker-release.json',{'pid':worker['pid'],'start_ticks':worker['start_ticks'],
             'postexec_observed_sha256':digest(directory/'worker-postexec-observed.json')})
        print('ACTUAL_WORKER_OBSERVED',worker['pid'],worker['start_ticks'],observed['executable'],
              observed['affinity'],observed['live_kernel']['memory.max'],observed['live_kernel']['memory.swap.max'],flush=True)
        while True:
            state=command(['systemctl','--user','show',unit,'--property=SubState','--value'])
            if state['stdout'].strip() in ('exited','failed','dead') or state['returncode']:
                break
            time.sleep(.25)
    except BaseException:
        save(directory/'owner-failure.json',{'traceback':traceback.format_exc(),'unit':unit})
        print(traceback.format_exc(),flush=True)
    properties=command(['systemctl','--user','show',unit,
        '--property=ActiveState,SubState,Result,ExecMainCode,ExecMainStatus,ExecMainStartTimestamp,ExecMainExitTimestamp,CPUUsageNSec,MemoryCurrent,MemoryPeak,MemoryMax,MemorySwapMax,CPUAffinity,ControlGroup,MainPID,NRestarts'])
    journal=command(['journalctl','--user','-u',unit,'--no-pager','-o','short-iso-precise'])
    with (directory/'owner-journal.txt').open('x') as stream:
        stream.write(journal['stdout']);stream.write(journal['stderr']);stream.flush();os.fsync(stream.fileno())
    save(directory/'owner-outcome.json',{'unit':unit,'properties':properties,'journal_argv':journal['argv'],
         'journal_sha256':digest(directory/'owner-journal.txt'),'journal_returncode':journal['returncode']})
    print(journal['stdout'],flush=True)
    print('SENSITIVITY_OWNER_TERMINAL',unit,properties['stdout'],flush=True)
    return 0 if (directory/'receipt.json').exists() else 2


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--generation',required=True)
    args=parser.parse_args()
    raise SystemExit(launch(args.generation))
