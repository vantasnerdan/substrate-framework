"""One capped CPU6 Native attempt; continue the original A43+A47 ledger."""
import json
import os
from pathlib import Path
import sys
from capture_command import capture,binding,save
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'

def main():
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('actual inherited Herdr context required')
    sources=[binding(p) for p in sorted(HERE.glob('*.py'))]
    prereg=json.loads((HERE/'computational-preregistration-002.json').read_bytes())
    if prereg['source_bindings']!=sources:raise RuntimeError('implementation changed after preregistration')
    unit='p254-a51-gamma-contour-002';prefix=HERE/'actual-native-002'
    pythonpath=':'.join([str(ROOT/'src'),*[str(HERE.parent/n) for n in ('0041','0043','0047','0039','0028','0027','0025','0024')],'/home/administrator/.cache/dft/p254-a28-runtime'])
    argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
        'systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec=800','--property=MemoryMax=3221225472','--property=MemorySwapMax=0',
        '--property=Restart=no','--property=CPUAffinity=6','--property=WorkingDirectory='+str(ROOT),
        '/usr/bin/env','PYTHONPATH='+pythonpath,
        *[k+'=1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')],
        '/usr/bin/taskset','-c','6',PYTHON,'-u',str(HERE.parent/'0047/full37_kernel_entry.py'),
        str(prefix)+'.kernel-start.json',unit,'6',str(HERE/'run_gamma_contour.py'),'/var/tmp/p254-a51-gamma-contour-002','780']
    save(HERE/'actual-native-source-input-bindings-002.json',dict(sources=sources,preregistration=binding(HERE/'computational-preregistration-002.json'),input_inventory=binding(HERE/'input-inventory-expanded-001.json'),
        inherited_herdr_env=os.environ['HERDR_ENV'],original_cumulative_seconds=7973.897332867935,remaining_before_seconds=1026.1026671320647))
    code=capture('actual-native-002',argv)
    print('A51_NATIVE_TERMINAL_002_EXIT_'+str(code),flush=True)
    return code

if __name__=='__main__':sys.exit(main())
