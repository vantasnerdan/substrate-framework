"""Full saved two-time dressed legs; original raw non-Gamma Gaussian bath."""
from pathlib import Path
import sys, json, hashlib, os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'src')]+[str(HERE.parent/n) for n in ('0055','0054','0053','0052','0051','0050','0048','0046','0045','0044','0043','0042','0041','0030','0029','0028','0027','0025','0024')]
import numpy as np
from general_moving_reference import SavedMovingReference
ORDERS=('greater','lesser','retarded')
PAIRS=tuple((a,b) for a in range(3) for b in range(3))
DEFAULT_DATASET=Path('/var/tmp/p254-a53-fixed-kernel-material-003/actual-material-step.npz')
DEFAULT_CATALOG=HERE.parent/'0053/actual-input-catalog-003.json'
def binding(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def save(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    with tmp.open('w') as f:json.dump(value,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(p)
def authenticate(entry):
    if binding(entry['path'])!=entry:raise ValueError('binding changed: '+entry['path'])
    return Path(entry['path'])
def spectral(ret,a,b):
    # theta(0)=1: the diagonal is the raw retarded jump, NOT R-R^dagger.
    return ret[a,b] if a>=b else -ret[b,a].swapaxes(-1,-2).conj()
def greater(ret,less,a,b):return less[a,b]+spectral(ret,a,b)
class DressedReference:
    def __init__(self,c,dataset,baseline_catalog,coefficient_binding):
        self.c=c;self.dataset=binding(dataset);self.catalog_binding=binding(baseline_catalog)
        cat=json.loads(Path(baseline_catalog).read_bytes());self.catalog=cat
        self.bare=SavedMovingReference(c,authenticate(cat['reference_bundle']),coefficient_binding)
        with np.load(dataset,allow_pickle=False) as z:
            self.times=z['times'];self.er=z['electron_ret'];self.el=z['electron_less'];self.dr=z['gamma_ret'];self.dl=z['gamma_less']
        nt=len(self.times)
        if not np.array_equal(self.times,self.bare.times):raise ValueError('exact source times required')
        if self.er.shape!=(nt,nt,128,12,12) or self.el.shape!=self.er.shape or self.dr.shape!=(nt,nt,114,114) or self.dl.shape!=self.dr.shape:raise ValueError('full saved two-time axes required')
        if not all(np.isfinite(x).all() for x in (self.er,self.el,self.dr,self.dl)):raise ValueError('nonfinite raw dressed state')
        self.bindings=[self.dataset,self.catalog_binding,*self.bare.bindings]
    def electron(self,a,b):return greater(self.er,self.el,a,b),self.el[a,b]
    def ordered(self,q,a,b,selection='positions'):
        if int(q)!=0:return self.bare.ordered(q,a,b,selection)
        ix=self.bare.position_indices(q,selection)
        if selection=='phase':raise ValueError('saved dressed dataset contains positions, not phase reconstruction')
        ix=np.arange(111) if selection=='ionic' else np.arange(111,114) if selection=='field' else np.arange(114)
        return greater(self.dr,self.dl,a,b)[np.ix_(ix,ix)],self.dl[a,b][np.ix_(ix,ix)]
    def boson(self,q,a,b):return self.ordered(q,a,b,'ionic')
    def close(self):self.bare.close()
def inputs(dataset,baseline_catalog):
    from joint_material_inputs import load_inputs
    c,baseline,old,bindings=load_inputs();old.close()
    return c,DressedReference(c,dataset,baseline_catalog,bindings['manifest']),bindings
