"""Consume the actual saved q1 cubic component, never a six-array bank."""
from pathlib import Path
import json
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0066'))
from column_common import CHANNELS, authenticate, file_binding, atomic_json


def consume(catalog):
    spec = json.loads(Path(catalog).read_bytes())
    if (spec['schema'] != 'p254-additive-history-component-v1' or
            spec['status'] != 'complete-owned-component-q1' or
            spec['component'] != 'cubic' or spec['completed_q'] != [1] or
            spec['full_mesh'] or spec['times_inverse_ry'] != list(range(31))):
        raise ValueError('complete genuine31 q1 cubic component required')
    values = {}
    for name in CHANNELS:
        entry = spec['arrays'][name]
        values[name] = np.load(authenticate({k:entry[k] for k in ('path','bytes','sha256')}),
                               mmap_mode='r', allow_pickle=False)
        if values[name].shape != (31,31,111,111) or values[name].dtype != np.dtype('complex128'):
            raise ValueError('canonical full111 complex128 history required')
    coverage = np.load(authenticate(spec['coverage']['array']), allow_pickle=False)
    if coverage.shape != (31,31) or not coverage.all():
        raise ValueError('actual persisted961 coverage required')
    probe = np.arange(1,112)/111 + 1j*np.arange(111)/113
    norms = {name:0. for name in CHANNELS}
    action = {name:0. for name in CHANNELS}
    causal_max = 0.
    for a in range(31):
        for b in range(31):
            for name in CHANNELS:
                cell = values[name][a,b]
                if not np.isfinite(cell).all():
                    raise ValueError('nonfinite saved numerical cell')
                norms[name] += float(np.vdot(cell,cell).real)
                action[name] += float(np.linalg.norm(cell@probe)**2)
            expected = values['greater'][a,b]-values['lesser'][a,b] if a>=b else 0.
            causal_max = max(causal_max,float(np.max(np.abs(values['retarded'][a,b]-expected))))
    reused_equal = True
    with np.load(authenticate(spec['paid_reuse']['column']), allow_pickle=False) as column:
        for name in CHANNELS:
            reused_equal = reused_equal and np.array_equal(values[name][:,0],column[name])
    with np.load(authenticate(spec['paid_reuse']['endpoints']), allow_pickle=False) as endpoint:
        for name in CHANNELS:
            reused_equal = reused_equal and np.array_equal(values[name][0,30],endpoint['old_'+name][0,1])
            reused_equal = reused_equal and np.array_equal(values[name][30,30],endpoint['old_'+name][1,1])
    if not reused_equal or causal_max != 0.:
        raise ValueError('literal paid reuse or theta(t>=s) violated')
    final_action = values['retarded'][30,29]@probe
    result = dict(schema='p254-a72-actual-saved-component-consumer-v1',
        catalog=file_binding(catalog),pairs_consumed=int(coverage.sum()),
        full_history_frobenius={k:float(np.sqrt(v)) for k,v in norms.items()},
        all_pair_probe_action_l2={k:float(np.sqrt(v)) for k,v in action.items()},
        final_new_pair=[30,29],final_new_retarded_action_norm=float(np.linalg.norm(final_action)),
        final_new_retarded_trace=[float(np.trace(values['retarded'][30,29]).real),float(np.trace(values['retarded'][30,29]).imag)],
        literal_paid33_cells_equal=reused_equal,causal_difference_max=causal_max,
        eligible_scope='externalGamma/full111/internalq1-partnerq3-only; not full64 or six-array functional')
    atomic_json(Path(catalog).parent/'consumer.json', result)
    return result


if __name__ == '__main__':
    print(json.dumps(consume(Path(sys.argv[1]))),flush=True)
