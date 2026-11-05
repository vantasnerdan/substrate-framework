"""Literal two-node order1 source algebra on paid bytes; no solver replay."""
from pathlib import Path
import json
import sys
import os
import resource
import numpy as np
from supplied_keldysh import A49
from command_receipts import binding


def maximum(a):
    return float(np.max(np.abs(a)))


def main(root):
    paid = Path('/var/tmp/p254-a53-fixed-kernel-material-001/fermions/slot-000')
    new = Path('/var/tmp/p254-a54-keldysh-002/actual-slot000/solution.npz')
    with np.load(paid/'input.npz',allow_pickle=False) as data:
        q = data['source_less']
        h = float(data['times'][1])
        assert int(data['order'])==1 and np.array_equal(data['times'],[0.,1.])
    with np.load(paid/'solution.npz',allow_pickle=False) as data:
        r,l,f,fc,fl,fcl = [data[name] for name in ('ret','less','f_ret','fcc_ret','f_less','fcc_less')]
    fl = fl.copy()
    fl[1,0] = -fcl[0,1].conj().T  # actual general companion, not old export's self continuation
    # cntr_vie2_impl.hpp:847-948 and convolution_impl's order1 Gregory weights.
    # At two nodes w(0,j)=0, w(1,0)=w(1,1)=1/2; F^R diagonal is zero.
    predicted01 = q[0,1]-.5*h*(fl[0,0]@r[1,0].conj().T+fl[0,1]@r[1,1].conj().T)
    predicted11 = q[1,1]-.5*h*(f[1,0]@l[0,1]+fl[1,0]@r[1,0].conj().T+fl[1,1]@r[1,1].conj().T)
    right11 = l[1,1]+.5*h*(r[1,0]@fcl[0,1]+r[1,1]@fcl[1,1]+l[1,0]@fc[1,0].conj().T)-q[1,1]
    with np.load(new,allow_pickle=False) as data:
        new_l,new_r = data['less'],data['ret']
    report = {'original_input':binding(paid/'input.npz'),'paid_failure':binding(paid/'solution.npz'),
              'new_solution':binding(new),
              'saved_F_R_diagonal_max':max(maximum(f[n,n]) for n in (0,1)),
              'actual_affinity':sorted(os.sched_getaffinity(0)),
              'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'thread_environment':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
              'old_source_algebra_G01_reconstruction_error':maximum(predicted01-l[0,1]),
              'old_source_algebra_G11_reconstruction_error':maximum(predicted11-l[1,1]),
              'predicted_one_sided_diagonal_exchange':maximum(predicted11+predicted11.conj().T),
              'paid_one_sided_diagonal_exchange':maximum(l[1,1]+l[1,1].conj().T),
              'old_right_Dyson_diagonal_residual':maximum(right11),
              'new_vs_old_retarded_difference':maximum(new_r-r),
              'new_vs_old_lesser_difference':maximum(new_l-l),
              'mechanism':'Literal upstream one-sided trapezoid equation reconstructs the paid output, including its non-anti-Hermitian diagonal. Its right companion equation is incompatible at the measured finite-step defect. No invalid-input hypothesis or reader tolerance is needed. Two-sided Langreth discretization preserves the invariant but is not the same finite-step equation; its exported residual remains accuracy debt.',
              'scope':'Exact source-backed two-node matrix arithmetic on paid bytes, double precision; not a rerun of the failed solver or an independent material exact evolution.'}
    Path(root,'coarse-source-forensics.json').write_text(json.dumps(report,indent=2)+'\n')
    assert report['old_source_algebra_G01_reconstruction_error']<1e-12
    assert report['old_source_algebra_G11_reconstruction_error']<1e-12
    print('A54_PAID_SOURCE_ALGEBRA_COMPLETE',flush=True)

if __name__ == '__main__':
    main(Path(sys.argv[1]))
