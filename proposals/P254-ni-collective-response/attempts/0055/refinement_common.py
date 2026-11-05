"""Three-node source-bound midpoint reference and literal old-pair mapping."""
from pathlib import Path
import hashlib,json,os,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'src')]+[str(HERE.parent/n) for n in ('0052','0051','0050','0048','0046','0045','0044','0043','0042','0041','0030','0029','0028','0027','0025','0024')]
import numpy as np
TIMES=np.array([0.,.5,1.])
NEW_PAIRS=[(0,1),(1,0),(1,1),(1,2),(2,1)]
OLD_MAP={0:0,2:1}
def binding(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def save(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_suffix(p.suffix+'.tmp')
    with temp.open('w') as f:
        json.dump(value,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    temp.replace(p)
def catalog():
    return json.loads((HERE.parent/'0053/actual-input-catalog-001.json').read_bytes())
def authenticate(entry):
    if binding(entry['path'])!=entry:raise ValueError('original byte binding changed: '+entry['path'])
    return Path(entry['path'])
def reference(c,coefficient_binding):
    from general_moving_reference import SavedMovingReference
    return SavedMovingReference(c,HERE/'reference/bundle.json',coefficient_binding)
def old_merge(old,new,slot_first=False):
    """Merge five genuine entries with four literal final paid entries."""
    out=np.empty((old.shape[0],3,3,*old.shape[3:]) if slot_first else (3,3,*old.shape[2:]),dtype=old.dtype)
    for a,oa in OLD_MAP.items():
        for b,ob in OLD_MAP.items():
            if slot_first:out[:,a,b]=old[:,oa,ob]
            else:out[a,b]=old[oa,ob]
    for k,(a,b) in enumerate(NEW_PAIRS):
        if slot_first:out[:,a,b]=new[:,k]
        else:out[a,b]=new[k]
    return out
