"""One bounded paid-byte source algebra Native, with costs never reset."""
from pathlib import Path
import json
import sys
from supplied_keldysh import A49
from command_receipts import binding, run_recorded
from native_campaign import HERE, PYTHON

root = Path(sys.argv[1])
root.mkdir(parents=True,exist_ok=False)
prior = [json.loads(p.read_text()) for p in root.parent.glob('p254-a54-keldysh-*/systemd-native.json')]
charged = sum(p['wall_seconds'] for p in prior)
if charged+60>1800:
    raise RuntimeError('cumulative Native grant exhausted')
record = {'charged_before':charged,'attempt_limit':60,'cumulative_limit':1800,
          'source':binding(HERE/'coarse_forensics.py'),'owned_pane':'w1:p6S',
          'prior_observed_marker':'A54_NATIVE_TERMINAL_002_EXIT_0'}
(root/'grant.json').write_text(json.dumps(record,indent=2)+'\n')
argv = ['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
        'systemd-run','--user','--wait','--pipe','--collect','--unit=p254-a54-keldysh-003.service',
        '--property=RuntimeMaxSec=60','--property=MemoryMax=3221225472','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=9','env','OPENBLAS_NUM_THREADS=1',
        'OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','BLIS_NUM_THREADS=1',
        'PYTHONDONTWRITEBYTECODE=1','taskset','-c','9',PYTHON,'-u',str(HERE/'coarse_forensics.py'),str(root)]
try:
    run_recorded(argv,root,'systemd-native',70)
    print('A54_NATIVE_TERMINAL_003_EXIT_0',flush=True)
except Exception:
    print('A54_NATIVE_TERMINAL_003_EXIT_NONZERO',flush=True)
    raise
