"""One-shot A62 native supervisor; exact argv and all exit streams retained."""
from pathlib import Path
import os
import sys
import json
import time
import subprocess
import hashlib
HERE=Path(__file__).resolve().parent
if os.environ.get('HERDR_ENV') != '1':
    raise RuntimeError('actual inherited Herdr context required')
herdr = [k+'='+v for k,v in os.environ.items() if k.startswith('HERDR_')]
suffix=sys.argv[1] if len(sys.argv)>1 else '001'
if suffix not in ('001','002','003'):
    raise ValueError('only the original invocation and explicit runtime/status amendments exist')
RECEIPTS=HERE/('receipts-'+suffix)
RECEIPTS.mkdir(exist_ok=False)
UNIT='p254-a62-history-smoke-'+suffix+'.service'
runtime=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']['python_environment_assignment']
argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
      'systemd-run','--user','--wait','--pipe','--collect','--unit='+UNIT,
      '--property=RuntimeMaxSec=570','--property=MemoryMax=3221225472','--property=MemorySwapMax=0',
      '--property=Restart=no','--property=CPUAffinity=9','env',runtime,
      'A62_NATIVE_UNIT='+UNIT,'A62_RECEIPTS_PATH='+str(RECEIPTS)]+herdr+[
      'OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','BLIS_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1',
      'PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(RECEIPTS/'native-time.txt'),
      '/home/administrator/projects/substrate-framework/.venv/bin/python','-u',str(HERE/'material_smoke.py')]
started=time.monotonic()
(RECEIPTS/'launch.json').write_text(json.dumps(dict(argv=argv,owner='A62 StreamedHistoryFunctional',
    unit=UNIT,cpu=9,memory_bytes=3221225472,swap_bytes=0,threads=1,runtime_cap_seconds=570,
    wall_start=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),stdout='stdout.txt',stderr='stderr.txt'),indent=2)+'\n')
print('A62_OWNED_BEGIN '+UNIT+' CPU9 BLAS1 3GiB noSwap',flush=True)
with (RECEIPTS/'stdout.txt').open('wb') as out,(RECEIPTS/'stderr.txt').open('wb') as err:
    result=subprocess.run(argv,cwd=HERE,stdout=out,stderr=err)
record=dict(argv=argv,exitcode=result.returncode,elapsed_seconds=time.monotonic()-started,
    wall_end=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
for name in ('stdout.txt','stderr.txt','native-time.txt'):
    p=RECEIPTS/name
    if p.exists():
        record[name]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
(RECEIPTS/'exit.json').write_text(json.dumps(record,indent=2)+'\n')
print('A62_OWNED_EXIT_'+suffix+'_'+str(result.returncode)+' seconds='+str(record['elapsed_seconds']),flush=True)
sys.exit(result.returncode)
