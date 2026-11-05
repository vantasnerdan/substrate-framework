"""Native final merger of actual disjoint complete saved outer-q contributions."""
import argparse
import contextlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import traceback


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--partition',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    m=subprocess.run(['systemctl','--user','show','p254-a46-full-q-merge-001.service',
        '--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,CPUAffinity,InvocationID'],capture_output=True,text=True)
    props=dict(x.split('=',1) for x in m.stdout.splitlines() if '=' in x)
    cg=Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    (out/'actual-startup.json').write_text(json.dumps(dict(manager=props,manager_stdout=m.stdout,
        manager_stderr=m.stderr,affinity=sorted(os.sched_getaffinity(0)),
        kernel_memory_max=(cg/'memory.max').read_text(),kernel_swap_max=(cg/'memory.swap.max').read_text()),indent=2))
    code=0
    with (out/'science.rawstdout.txt').open('w') as stdout,(out/'science.rawstderr.txt').open('w') as stderr:
        with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):
            try:
                from all_slot_coulomb_partition import merge_partitions
                result=merge_partitions(a.partition,out/'result')
                receipt=dict(status='FULL_ALL_Q_COMPLETE',archive=result['record']['archive'],
                    counts=result['record']['counts'],decoded_norms=result['decoded_norms'])
                (out/'complete-saved-consumption.json').write_text(json.dumps(receipt,indent=2))
                print(json.dumps(receipt),flush=True)
            except BaseException:
                code=1; traceback.print_exc()
    terminal=dict(exit=code,seconds=time.monotonic()-start,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'actual-native-exit.json').write_text(json.dumps(terminal,indent=2))
    print(json.dumps(terminal),flush=True)
    return code
if __name__=='__main__': sys.exit(main())
