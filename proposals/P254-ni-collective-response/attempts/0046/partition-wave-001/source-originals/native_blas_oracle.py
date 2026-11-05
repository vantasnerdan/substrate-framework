"""Durable capture for the one newly changed BLAS source oracle."""
import contextlib
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
MODE = os.environ.get('A46_OPTIMIZATION','blas')
OUT = HERE/({'blas':'optimization-native-001','shared':'shared-native-001','partition':'partition-native-001'}[MODE])
UNIT = 'p254-a46-'+MODE+'-oracle-001.service'
SCRIPT = {'blas':'oracle_all_slot_blas.py','shared':'oracle_all_slot_shared.py','partition':'oracle_q_partition.py'}[MODE]

def main():
    OUT.mkdir(exist_ok=True)
    start = time.monotonic()
    m = subprocess.run(['systemctl','--user','show',UNIT,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,CPUAffinity,InvocationID'],capture_output=True,text=True)
    props = dict(x.split('=',1) for x in m.stdout.splitlines() if '=' in x)
    cg = Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    (OUT/'startup.json').write_text(json.dumps(dict(manager=props,stdout=m.stdout,stderr=m.stderr,
        affinity=sorted(os.sched_getaffinity(0)),kernel_memory_max=(cg/'memory.max').read_text(),
        kernel_swap_max=(cg/'memory.swap.max').read_text()),indent=2))
    code = 0
    with (OUT/'science.rawstdout.txt').open('w') as out,(OUT/'science.rawstderr.txt').open('w') as err:
        with contextlib.redirect_stdout(out),contextlib.redirect_stderr(err):
            try:
                runpy.run_path(str(HERE/SCRIPT),run_name='__main__')
            except BaseException:
                code = 1
                traceback.print_exc()
    receipt = dict(exit=code,seconds=time.monotonic()-start,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (OUT/'native-exit.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt),flush=True)
    return code
if __name__ == '__main__':
    sys.exit(main())
