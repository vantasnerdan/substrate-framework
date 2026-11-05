"""Owned native startup, raw streams and terminal receipt for the granted oracle."""
import contextlib
import datetime
import json
import os
from pathlib import Path
import resource
import runpy
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
UNIT = 'p254-a46-sensitive-all-slot-oracle-001.service'

def save(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2,sort_keys=True))

def main():
    started = time.monotonic()
    manager_argv = ['systemctl','--user','show',UNIT,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState,InvocationID']
    manager = subprocess.run(manager_argv,capture_output=True,text=True)
    properties = dict(line.split('=',1) for line in manager.stdout.splitlines() if '=' in line)
    kernel = {}
    cg = Path('/sys/fs/cgroup')/properties.get('ControlGroup','').lstrip('/')
    for name in ('memory.max','memory.swap.max','cpuset.cpus.effective'):
        try:
            kernel[name] = dict(availability='available',value=(cg/name).read_text().strip())
        except OSError as e:
            kernel[name] = dict(availability='unavailable',reason=str(e))
    save('native-startup.json',dict(unit=UNIT,observed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        actual_affinity=sorted(os.sched_getaffinity(0)),manager=properties,
        manager_argv=manager_argv,manager_exit=manager.returncode,
        manager_rawstdout=manager.stdout,manager_rawstderr=manager.stderr,
        own_kernel_interfaces=kernel,pid=os.getpid()))
    code = 0
    with (HERE/'oracle.rawstdout.txt').open('w',buffering=1) as out, (HERE/'oracle.rawstderr.txt').open('w',buffering=1) as err:
        with contextlib.redirect_stdout(out),contextlib.redirect_stderr(err):
            try:
                runpy.run_path(str(HERE/'oracle_all_slot_coulomb.py'),run_name='__main__')
            except BaseException:
                code = 1
                traceback.print_exc()
    receipt = dict(unit=UNIT,exit=code,seconds=time.monotonic()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        ended_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        stdout=str(HERE/'oracle.rawstdout.txt'),stderr=str(HERE/'oracle.rawstderr.txt'))
    save('native-exit.json',receipt)
    print(json.dumps(receipt),flush=True)
    return code

if __name__ == '__main__':
    sys.exit(main())
