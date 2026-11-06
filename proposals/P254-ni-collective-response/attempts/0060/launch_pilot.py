"""Durable native argv/stdout/stderr/exit supervisor for the owned A60 unit."""
from pathlib import Path
import os,sys,json,time,subprocess,hashlib
HERE=Path(__file__).resolve().parent
cache = len(sys.argv)>1 and sys.argv[1]=='cache64'
suffix = '002' if cache else '001'
script = 'run_cache64.py' if cache else 'run_pilot.py'
unit = 'p254-a60-bounded-ee-'+suffix+'.service'
RECEIPTS=HERE/('receipts-'+suffix)
RECEIPTS.mkdir(exist_ok=False)
argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus','systemd-run','--user','--wait','--pipe','--collect','--unit='+unit,'--property=RuntimeMaxSec='+('180' if cache else '570'),'--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=8','env','HERDR_ENV=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1','/usr/bin/time','-v','-o',str(RECEIPTS/'native-time.txt'),'/home/administrator/projects/substrate-framework/.venv/bin/python','-u',str(HERE/script)]
start=time.monotonic()
(RECEIPTS/'launch.json').write_text(json.dumps(dict(argv=argv,wall_start=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),owner='A60 BoundedTimeEE',cpu=8,memory_bytes=3221225472,swap_bytes=0,prior_ee_seconds=15722.9054025717266,native_cumulative_ceiling_seconds=600,stdout='stdout.txt',stderr='stderr.txt'),indent=2)+'\n')
print('A60_OWNED_NATIVE_BEGIN unit='+unit+' CPU8 BLAS1 3GiB noSwap',flush=True)
with (RECEIPTS/'stdout.txt').open('wb') as out,(RECEIPTS/'stderr.txt').open('wb') as err:
    proc=subprocess.run(argv,cwd=HERE,stdout=out,stderr=err)
record=dict(argv=argv,exitcode=proc.returncode,elapsed_seconds=time.monotonic()-start,wall_end=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
for name in ('stdout.txt','stderr.txt','native-time.txt'):
    p=RECEIPTS/name
    if p.exists():record[name]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
(RECEIPTS/'exit.json').write_text(json.dumps(record,indent=2)+'\n')
print('A60_OWNED_NATIVE_EXIT_'+suffix+' code='+str(proc.returncode)+' elapsed='+str(record['elapsed_seconds']),flush=True)
raise SystemExit(proc.returncode)
