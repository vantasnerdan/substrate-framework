"""A52 durable control-command evidence; never synthesize Herdr scope."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
HERE=Path(__file__).resolve().parent
OUT=Path('/var/tmp/p254-a52-control-001')
def binding(path):
    p=Path(path);h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def save(path,obj):
    with Path(path).open('x') as f:
        json.dump(obj,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def capture(label,argv):
    OUT.mkdir(exist_ok=True)
    save(OUT/(label+'.argv.json'),dict(argv=argv,start_unix=time.time(),HERDR_ENV=os.environ.get('HERDR_ENV'),source=binding(__file__)))
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('actual HERDR_ENV=1 required')
    began=time.monotonic()
    with (OUT/(label+'.stdout')).open('xb') as out,(OUT/(label+'.stderr')).open('xb') as err:
        p=subprocess.run(argv,stdout=out,stderr=err,check=False)
    save(OUT/(label+'.exit.json'),dict(exitcode=p.returncode,elapsed_seconds=time.monotonic()-began,stdout=binding(OUT/(label+'.stdout')),stderr=binding(OUT/(label+'.stderr'))))
    print((OUT/(label+'.stdout')).read_text(),end='');print((OUT/(label+'.stderr')).read_text(),end='',file=sys.stderr)
    return p.returncode
if __name__=='__main__':sys.exit(capture(sys.argv[1],sys.argv[2:]))
