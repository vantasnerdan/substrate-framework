"""Read actual saved numeric boundary and retain CAR/positivity, never clip."""
from pathlib import Path
import json
import sys
import numpy as np
from supplied_keldysh import load_solution
from command_receipts import binding


def consume(path):
    result = load_solution(path)
    ret, less = result['ret'], result['less']
    stat = int(result['statistics'])
    diagonal = np.array([less[n,n] for n in range(len(less))])
    state = (-1j if stat == -1 else 1j)*diagonal
    eigenvalues = np.linalg.eigvalsh(state)
    report = {'solution':binding(path), 'times':result['times'].tolist(), 'statistics':stat,
              'all_finite':bool(np.isfinite(ret).all() and np.isfinite(less).all()),
              'raw_lesser_exchange_max':float(np.max(np.abs(less+less.swapaxes(0,1).swapaxes(2,3).conj()))),
              'causal_upper_retarded_max':float(np.max(np.abs(ret[np.triu_indices(len(ret),1)]))),
              'state_eigenvalue_min':float(eigenvalues.min()), 'state_eigenvalue_max':float(eigenvalues.max()),
              'initial_state_real':state[0].real.tolist(), 'initial_state_imag':state[0].imag.tolist(),
              'final_state_real':state[-1].real.tolist(), 'final_state_imag':state[-1].imag.tolist(),
              'initial_trace_real':float(np.trace(state[0]).real), 'final_trace_real':float(np.trace(state[-1]).real),
              'observable':'density -iG<' if stat == -1 else 'physical Cartesian position covariance +iD<',
              'accuracy_scope':'Invariant diagnostics only; material finite-step accuracy requires independent refinement.'}
    if stat == -1:
        spectral = np.array([1j*ret[n,n] for n in range(len(ret))])
        holes = spectral-state
        report.update(raw_CAR_defect_max=float(np.max(np.abs(spectral-np.eye(ret.shape[-1])))),
                      hole_eigenvalue_min=float(np.linalg.eigvalsh(holes).min()),
                      occupation_above_one_max=float(max(0.,eigenvalues.max()-1.)))
    else:
        report['equal_time_position_retarded_max'] = float(max(np.max(np.abs(ret[n,n])) for n in range(len(ret))))
        report['uncertainty_scope'] = 'Position block alone cannot establish full x/p Heisenberg uncertainty.'
    return report

if __name__ == '__main__':
    Path(sys.argv[2]).write_text(json.dumps(consume(sys.argv[1]),indent=2)+'\n')
