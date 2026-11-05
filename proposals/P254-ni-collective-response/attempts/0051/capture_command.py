"""Durable command receipts; require the actual inherited Herdr context."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent

def binding(path):
    path=Path(path);h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest())

def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())

def capture(label,argv):
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('actual HERDR_ENV=1 required; no synthetic context')
    prefix=HERE/label;start=time.monotonic()
    save(str(prefix)+'.start.json',dict(argv=argv,source=binding(__file__),herdr_env=os.environ['HERDR_ENV'],
        context={k:v for k,v in os.environ.items() if k.startswith(('HERDR_','MIDTERM_'))},start_unix=time.time()))
    with open(str(prefix)+'.stdout','xb') as out,open(str(prefix)+'.stderr','xb') as err:
        result=subprocess.run(argv,stdout=out,stderr=err,check=False)
    receipt=dict(argv=argv,exit_code=result.returncode,elapsed_seconds=time.monotonic()-start,
        stdout=binding(str(prefix)+'.stdout'),stderr=binding(str(prefix)+'.stderr'))
    save(str(prefix)+'.exit.json',receipt)
    print(json.dumps(receipt),flush=True)
    return result.returncode

if __name__=='__main__':sys.exit(capture(sys.argv[1],sys.argv[2:]))
