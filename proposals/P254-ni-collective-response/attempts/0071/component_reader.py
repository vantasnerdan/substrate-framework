"""Read only the q1 raw EI Sigma component; never admit it as a full functional."""
from pathlib import Path
import json
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0062'))
from history_banks import ORDERS, authenticate, binding

SCHEMA = 'p254-additive-history-component-v1'
SHAPE = (31, 31, 128, 12, 12)


class EIComponent:
    def __init__(self, catalog, *, require_complete=True):
        self.record = json.loads(Path(catalog).read_bytes())
        r = self.record
        if (r['schema'] != SCHEMA or r['component'] != 'mixed_ei' or
                r['owned_q'] != [1] or r['full_mesh'] is not False or
                r['times_inverse_ry'] != list(range(31)) or
                r['axes'] != ['time','time','slot','band','band']):
            raise ValueError('not the original q1 raw EI component')
        pairs = list(map(tuple, r['coverage']['completed_pairs']))
        if len(pairs) != len(set(pairs)) or any(not (0 <= a < 31 and 0 <= b < 31) for a,b in pairs):
            raise ValueError('duplicate or invalid pair coverage')
        complete = set(pairs) == {(a,b) for a in range(31) for b in range(31)}
        if (r['status'] != ('complete-owned-component-q1' if complete else 'partial-owned-component-q1') or
                r['completed_q'] != ([1] if complete else [])):
            raise ValueError('coverage/status mismatch')
        if require_complete and not complete:
            raise ValueError('incomplete owned component')
        if set(r['sigma']) != set(ORDERS):
            raise ValueError('three canonical contour orders required')
        self.pairs = pairs
        self.arrays = {}
        try:
            for order in ORDERS:
                item = r['sigma'][order]
                arr = np.load(authenticate(item), mmap_mode='r', allow_pickle=False)
                self.arrays[order] = arr
                if arr.shape != SHAPE or arr.dtype != np.complex128 or item['shape'] != list(SHAPE) or item['dtype'] != 'complex128':
                    raise ValueError('canonical shape/dtype mismatch')
        except BaseException:
            self.close()
            raise

    def close(self):
        for arr in self.arrays.values():
            arr._mmap.close()
        self.arrays.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def consume(catalog, *, require_complete=True):
    with EIComponent(catalog, require_complete=require_complete) as saved:
        totals = {order: 0.0 for order in ORDERS}
        maxima = {order: 0.0 for order in ORDERS}
        causal_defect = 0.0
        for a,b in saved.pairs:
            values = {order: arr[a,b] for order,arr in saved.arrays.items()}
            for order, value in values.items():
                if not np.isfinite(value).all():
                    raise ValueError('nonfinite persisted component pair')
                totals[order] += float(np.vdot(value,value).real)
                maxima[order] = max(maxima[order],float(np.max(np.abs(value))))
            expected = values['greater']-values['lesser'] if a >= b else np.zeros_like(values['retarded'])
            causal_defect = max(causal_defect,float(np.max(np.abs(values['retarded']-expected))))
        sample = {}
        if (30,0) in saved.pairs:
            for order,arr in saved.arrays.items():
                z = complex(arr[30,0,65,0,0])
                sample[order] = [z.real,z.imag]
        return dict(schema='p254-a71-saved-component-consumer-v1', catalog=binding(catalog),
            status=saved.record['status'], completed_pairs=len(saved.pairs),
            actual_scalar_elements_per_order=len(saved.pairs)*128*12*12,
            all128_slots_consumed=True, both_spins=True, all12_bands=True,
            times_inverse_ry=saved.record['times_inverse_ry'], owned_q=[1], full_mesh=False,
            frobenius_norm={key:float(np.sqrt(value)) for key,value in totals.items()},
            max_abs=maxima, causal_max_abs_defect=causal_defect,
            sample_30_0_slot65_band0=sample,
            scope='Actual authenticated three-array q1 raw EI Sigma; not EE, cubic, joint Gamma, all64q or a full solve')
