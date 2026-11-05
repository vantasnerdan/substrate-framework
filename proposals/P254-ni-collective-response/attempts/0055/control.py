"""Capture actual Herdr control calls without synthetic terminal context."""
import os,sys,subprocess,time,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main(name,*argv):
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('actual Herdr required')
    output=HERE/'control';output.mkdir(exist_ok=True)
    source=dict(argv=list(argv),start_unix=time.time(),HERDR={k:v for k,v in os.environ.items() if k.startswith('HERDR_')})
    (output/(name+'.start.json')).write_text(json.dumps(source,indent=2)+'\n')
    start=time.monotonic()
    with (output/(name+'.stdout')).open('xb') as out,(output/(name+'.stderr')).open('xb') as err:
        p=subprocess.run(argv,stdout=out,stderr=err)
    (output/(name+'.exit.json')).write_text(json.dumps(dict(exitcode=p.returncode,elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
    print((output/(name+'.stdout')).read_text(),end='');print((output/(name+'.stderr')).read_text(),end='',file=sys.stderr)
    return p.returncode
if __name__=='__main__':sys.exit(main(*sys.argv[1:]))
