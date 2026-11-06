"""Actual literal q1-full plus original q0-nine summed-bank consumer."""
from pathlib import Path
import json
import os
import sys
import time
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0077'))
from component_accumulator import reopen
from history_banks import binding, authenticate, save, ORDERS


def main(case_path):
    case = json.loads(Path(case_path).read_bytes())
    if case['cpu'] not in range(6,11) or sorted(os.sched_getaffinity(0)) != [case['cpu']] or os.environ.get('HERDR_PANE_ID') != case['pane_id']:
        raise ValueError('owned authorized actual runtime required')
    started = time.monotonic()
    opened = reopen(case['bank'])
    state = opened['state']
    # This specific check is before any additional q admission; reject accidental
    # scope drift rather than silently comparing only a subset of a newer sum.
    for component in ('ee','mixed_ei','cubic'):
        if any(item['component']==component and item['q'] not in (0,1) for item in state['imports']):
            raise ValueError('new q already admitted; literal paid-nine case must enumerate it')
    prereg = json.loads((HERE.parent/'0077/preregistration.json').read_bytes())
    q1 = {'ee':'/var/tmp/p254-a73-ee-q1-canonical-001/catalog.json','mixed_ei':str(authenticate(prereg['EI_catalog'])),'cubic':str(authenticate(prereg['cubic_catalog']))}
    findings=[]
    for component in ('ee','mixed_ei','cubic'):
        zero_path=HERE.parent/('0076/q0-'+component+'-001/catalog.json')
        zero=json.loads(zero_path.read_bytes())
        one=json.loads(Path(q1[component]).read_bytes())
        key,kind=('arrays','cubic') if component=='cubic' else ('sigma','sigma')
        per_order={}
        grid=state['times_inverse_ry']
        for order in ORDERS:
            source0=np.load(authenticate(zero[key][order]),mmap_mode='r',allow_pickle=False)
            source1=np.load(authenticate(one[key][order]),mmap_mode='r',allow_pickle=False)
            bank=np.load(authenticate(state['families'][component]['arrays'][kind][order]),mmap_mode='r',allow_pickle=False)
            try:
                maximum=0.0
                count=nonzero=0
                for t in (0,15,30):
                    for u in (0,15,30):
                        actual=bank[grid.index(t),grid.index(u)]
                        literal=source1[one['times_inverse_ry'].index(t),one['times_inverse_ry'].index(u)]+source0[zero['times_inverse_ry'].index(t),zero['times_inverse_ry'].index(u)]
                        difference=float(np.max(np.abs(actual-literal)))
                        maximum=max(maximum,difference)
                        if not np.array_equal(actual,literal):
                            raise ValueError('literal q1+q0 paid sum differs: '+str((component,order,t,u,difference)))
                        count+=actual.size
                        nonzero+=int(np.count_nonzero(actual))
                per_order[order]=dict(compared_entries=count,nonzero_entries=nonzero,max_abs_difference=maximum,bitwise_equal=True)
            finally:
                source0._mmap.close()
                source1._mmap.close()
                bank._mmap.close()
        findings.append(dict(component=component,q0_catalog=binding(zero_path),q1_catalog=binding(q1[component]),orders=per_order))
    final=reopen(case['bank'])
    if final['pointer']!=opened['pointer']:
        raise ValueError('sole-writer pointer changed during saved consumer')
    save(case['result'],dict(schema='p254-a82-literal-paid-nine-consumer-v1',status='actual-literal-sum-pass',pointer=opened['pointer'],scope=dict(times=[0,15,30],ordered_pairs=9,components=['ee','mixed_ei','cubic'],orders=list(ORDERS),CAR_slots=128,bands=12,cubic_full_coordinates=111,physical_q=[0,1],contractions_recomputed=False),findings=findings,elapsed_seconds=time.monotonic()-started))
    print('A82_LITERAL_PAID_NINE_PASS',flush=True)


if __name__=='__main__':
    main(sys.argv[1])
