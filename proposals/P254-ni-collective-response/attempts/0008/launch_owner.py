#!/usr/bin/env python3
"""Launch one durable systemd owner and preserve actual terminal/resource evidence."""
import json
import os
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
UNIT='p254-w1-primary-boundary-0008-v1.service'
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
ENV={**os.environ,'XDG_RUNTIME_DIR':'/run/user/1000','DBUS_SESSION_BUS_ADDRESS':'unix:path=/run/user/1000/bus'}


def command(argv):
    completed=subprocess.run(argv,env=ENV,text=True,capture_output=True)
    return {'argv':argv,'returncode':completed.returncode,'stdout':completed.stdout,'stderr':completed.stderr}


if __name__=='__main__':
    if (HERE/'owner-launch.json').exists():
        raise SystemExit('Distinct owner generation required; never blindly relaunch')
    prepare=command([PYTHON,str(HERE/'run_oracle.py'),'--prepare'])
    print(prepare['stdout'],flush=True)
    args=['systemd-run','--user','--unit='+UNIT,'--property=Type=exec','--property=CPUAffinity=5',
          '--property=MemoryMax=2G','--property=Restart=no','--property=RemainAfterExit=yes',
          '--property=KillMode=control-group','--working-directory='+str(ROOT)]
    for key,value in {'PYTHONPATH':str(ROOT/'src'),'PYTHONUNBUFFERED':'1','OMP_NUM_THREADS':'1',
                      'MKL_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'}.items():
        args+=['--setenv='+key+'='+value]
    args += [PYTHON,str(HERE/'run_oracle.py')]
    launch=command(args) if prepare['returncode']==0 else {'returncode':2,'stderr':'Source registration failed'}
    with (HERE/'owner-launch.json').open('x') as stream:
        json.dump({'prepare':prepare,'launch':launch,'unit':UNIT},stream,indent=2)
    print('PRIMARY_OWNER_LAUNCH',UNIT,launch,flush=True)
    if launch['returncode']:
        raise SystemExit(2)
    while True:
        observed=command(['systemctl','--user','show',UNIT,'--property=SubState','--value'])
        if observed['stdout'].strip() in ('exited','failed','dead') or observed['returncode']:
            break
        time.sleep(0.2)
    properties=command(['systemctl','--user','show',UNIT,'--property=ActiveState,SubState,Result,ExecMainCode,ExecMainStatus,ExecMainStartTimestamp,ExecMainExitTimestamp,CPUUsageNSec,MemoryCurrent,MemoryPeak,MemoryMax,CPUAffinity,ControlGroup,MainPID,NRestarts'])
    journal=command(['journalctl','--user','-u',UNIT,'--no-pager','-o','short-iso-precise'])
    with (HERE/'owner-outcome.json').open('x') as stream:
        json.dump({'unit':UNIT,'properties':properties,'journal_command':journal['argv'],
                   'journal_returncode':journal['returncode'],'post_gc_memory_policy':'Only observed fields retained; no reconstructed post-GC values.'},stream,indent=2)
    with (HERE/'owner-journal.txt').open('x') as stream:
        stream.write(journal['stdout']);stream.write(journal['stderr'])
    print(journal['stdout'],flush=True)
    print('PRIMARY_OWNER_BOUNDARY_SNAPSHOT',UNIT,properties['stdout'],flush=True)
