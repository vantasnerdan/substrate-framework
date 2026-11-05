"""Authenticated global joint component reader; never admits a full six-array functional."""
from pathlib import Path
import json
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0062'))
from history_banks import ORDERS, authenticate, binding


class JointComponent:
    def __init__(self,catalog,*,require_complete=True):
        self.record = json.loads(Path(catalog).read_bytes())
        r = self.record
        if (r['schema'] != 'p254-additive-history-component-v1' or r['component'] != 'joint' or
            r['q_additive'] is not False or r['full_solve_eligible'] is not False or
            r['full_mesh'] is not False or r['times_inverse_ry'] != list(range(31)) or
            r['sigma_convention'] != 'joint_minus_ionic_Gamma' or r['pi_convention'] != 'full_joint_Gamma_once' or
            (r['slots'],r['bands'],r['spins'],r['gamma_coordinates']) != (128,12,2,114)):
            raise ValueError('not the original global jointGamma component')
        self.pairs = list(map(tuple,r['coverage']['completed_pairs']))
        all_pairs = {(a,b) for a in range(31) for b in range(31)}
        if len(self.pairs) != len(set(self.pairs)) or not set(self.pairs) <= all_pairs:
            raise ValueError('duplicate or invalid pair coverage')
        complete = set(self.pairs) == all_pairs
        if r['status'] != ('complete-global-joint-component' if complete else 'partial-global-joint-component') or r['full_time_square'] != complete:
            raise ValueError('coverage/status mismatch')
        if require_complete and not complete:
            raise ValueError('incomplete global component')
        self.banks = {}
        try:
            for kind,tail,axes in (('sigma',(128,12,12),['time','time','slot','band','band']),
                                   ('pi',(114,114),['time','time','coordinate','coordinate'])):
                if set(r[kind]) != set(ORDERS) or r['axes'][kind] != axes:
                    raise ValueError('canonical contour orders/axes required')
                self.banks[kind] = {}
                for order in ORDERS:
                    item = r[kind][order]
                    arr = np.load(authenticate(item),mmap_mode='r',allow_pickle=False)
                    self.banks[kind][order] = arr
                    if arr.shape != (31,31,*tail) or arr.dtype != np.complex128 or item['shape'] != list(arr.shape) or item['dtype'] != 'complex128':
                        raise ValueError('canonical shape/dtype mismatch')
        except BaseException:
            self.close()
            raise

    def close(self):
        for bank in self.banks.values():
            for arr in bank.values():
                arr._mmap.close()
        self.banks.clear()

    def __enter__(self):
        return self

    def __exit__(self,*exc):
        self.close()


def consume(catalog,*,require_complete=True):
    with JointComponent(catalog,require_complete=require_complete) as saved:
        norms, maxima, defects, samples = {}, {}, {}, {}
        for kind,bank in saved.banks.items():
            totals = {order:0.0 for order in ORDERS}
            maxima[kind] = {order:0.0 for order in ORDERS}
            defects[kind] = 0.0
            for a,b in saved.pairs:
                values = {order:arr[a,b] for order,arr in bank.items()}
                for order,value in values.items():
                    if not np.isfinite(value).all():
                        raise ValueError('nonfinite persisted pair')
                    totals[order] += float(np.vdot(value,value).real)
                    maxima[kind][order] = max(maxima[kind][order],float(np.max(np.abs(value))))
                residual = values['retarded']-(values['greater']-values['lesser']) if a>=b else values['retarded']
                defects[kind] = max(defects[kind],float(np.max(np.abs(residual))))
            norms[kind] = {order:float(np.sqrt(total)) for order,total in totals.items()}
            if (30,0) in saved.pairs:
                indices = (30,0,65,0,0) if kind=='sigma' else (30,0,111,113)
                samples[kind] = {order:[float(arr[indices].real),float(arr[indices].imag)] for order,arr in bank.items()}
        return dict(schema='p254-a75-saved-global-joint-consumer-v1',catalog=binding(catalog),
            status=saved.record['status'],component='joint',q_additive=False,completed_pairs=len(saved.pairs),
            actual_scalar_elements_per_order=dict(sigma=len(saved.pairs)*128*12*12,pi=len(saved.pairs)*114*114),
            all128_slots_consumed=True,both_spins=True,all12_bands=True,all111_plus3_coordinates_consumed=True,
            times_inverse_ry=saved.record['times_inverse_ry'],full_solve_eligible=False,
            frobenius_norm=norms,max_abs=maxima,causal_max_abs_defect=defects,sample_30_0=samples,
            scope=saved.record['scope'])
