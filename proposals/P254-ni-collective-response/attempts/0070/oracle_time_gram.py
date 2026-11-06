"""Immediate exact independent h/mode-sum oracle and exposing mutations."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0064'))
from time_gram import time_gram


def check():
    nt,m,nh,modes=2,2,2,3
    left=np.empty((nt,m,nh,modes),object)
    right=np.empty((nh,modes,nt,m),object)
    for index in np.ndindex(left.shape):left[index]=sp.Symbol('L_'+'_'.join(map(str,index)))
    for index in np.ndindex(right.shape):right[index]=sp.Symbol('R_'+'_'.join(map(str,index)))
    expected=np.empty((nt,nt,m,m),object)
    for t,s,i,j in np.ndindex(expected.shape):
        expected[t,s,i,j]=sum(left[t,i,h,a]*right[h,a,s,j]
            for h in range(nh) for a in range(modes))
    actual=time_gram(left,right)
    differences=[sp.expand(a-b) for a,b in zip(actual.flat,expected.flat)]
    if any(value!=0 for value in differences):raise AssertionError('ordinary bilinear time identity failed')
    absent=right.copy();absent[-1,-1]=0
    mutations={
        'manufactured_adjoint':time_gram(left,np.vectorize(sp.conjugate,otypes=[object])(right)),
        'right_time_output_swapped':time_gram(left,right.swapaxes(2,3)),
        'missing_inner_mode':time_gram(left,absent),
        'extra_half_normalization':actual/2,
    }
    sensitive={name:any(sp.expand(a-b)!=0 for a,b in zip(value.flat,expected.flat))
               for name,value in mutations.items()}
    if not all(sensitive.values()):raise AssertionError('non-exposing time-Gram mutation')
    records=[]
    for path in (Path(sys.modules['time_gram'].__file__),Path(__file__)):
        payload=path.read_bytes();records.append(dict(path=str(path),bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest()))
    return dict(schema='p254-time-gram-exact-block-oracle-v1',status='passed',
        identity='Independent explicit sum over every h/mode equals actual ordinary reshape/GEMM output at every t,s,i,j',
        dimensions=dict(times=nt,outputs=m,inner_h=nh,modes=modes),
        exact_entrywise_differences_zero=True,sensitivity=sensitive,source_bindings=records,
        limit='Finite exact layout/bilinear identity, not physical calibration or floating point enclosure.')


if __name__=='__main__':
    result=check();(HERE/'symbolic-findings.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
