"""One actual CPU7 all31 bare force production and saved-path consumer."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import traceback
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in
    ('0063','0062','0060','0057','0055','0052','0051','0050','0047','0046','0044','0045','0041','0030','0029')]
from produce_force_cache import produce
from history_storage import atomic_json, binding
OUT = Path('/var/tmp/p254-a74-full31-bare-force-001')


def main():
    group = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
        threads={key:os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS',
            'MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        limits={key:(cg/key).read_text().strip() for key in ('memory.max','memory.swap.max')},
        HERDR={key:value for key,value in os.environ.items() if key.startswith('HERDR_')},
        argv=sys.argv,pid=os.getpid(),cgroup=group)
    atomic_json(HERE/'receipts-001/runtime.json',runtime)
    if (runtime['affinity'] != [7] or set(runtime['threads'].values()) != {'1'} or
            runtime['limits'] != {'memory.max':'3221225472','memory.swap.max':'0'} or
            runtime['HERDR'].get('HERDR_PANE_ID') != 'w1:p7F' or os.environ.get('HERDR_ENV') != '1'):
        raise ValueError('actual owned one-thread CPU7/3GiB/noSwap context differs')
    result = produce(OUT)
    atomic_json(HERE/'production-result.json',result)
    print('A74_FORCE_NATIVE_COMPLETE '+json.dumps(result),flush=True)


if __name__ == '__main__':
    code = 0
    started = time.monotonic()
    try:
        main()
    except BaseException:
        code = 1
        traceback.print_exc()
    finally:
        atomic_json(HERE/'receipts-001/native-exit.json',dict(returncode=code,
            elapsed_seconds=time.monotonic()-started,
            maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        print('A74_FORCE_NATIVE_TERMINAL_'+str(code),flush=True)
    raise SystemExit(code)
