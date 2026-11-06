"""Sole-writer actual settled-component admission; no physics evaluation."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0077'))
from component_accumulator import accumulate, reopen, publish, _normalize, _reference, _digest, _summary
from history_banks import binding, authenticate, save, ORDERS


def consume(root, catalog_path):
    catalog = json.loads(Path(catalog_path).read_bytes())
    component, q, grid, pairs, selectors, kinds = _normalize(catalog, _reference()[0])
    opened = reopen(root)
    family = opened['state']['families'][component]
    coverage = np.load(authenticate(family['coverage']), mmap_mode='r', allow_pickle=False)
    digests = np.load(authenticate(family['digests']), mmap_mode='r', allow_pickle=False)
    sources = {kind:{order:np.load(authenticate(desc[order]), mmap_mode='r', allow_pickle=False) for order in ORDERS} for kind,desc in kinds.items()}
    checked = 0
    try:
        for a,b in pairs:
            for kind,slot,external in selectors:
                key = q,grid[a],grid[b],external
                expected = _digest([sources[kind][order][a,b,slot] for order in ORDERS])
                if not coverage[key] or digests[key].tobytes() != expected:
                    raise ValueError('saved contribution digest/coverage mismatch: '+str(key))
                checked += 1
        # Authenticate and actually open/read each committed summed contour array.
        arrays = []
        for kind,orders in family['arrays'].items():
            for order,descriptor in orders.items():
                array = np.load(authenticate(descriptor), mmap_mode='r', allow_pickle=False)
                try:
                    finite = all(bool(np.isfinite(array[i]).all()) for i in range(array.shape[0]))
                    if not finite:
                        raise ValueError('nonfinite committed summed array')
                    arrays.append(dict(kind=kind,order=order,binding=descriptor,shape=list(array.shape),finite=finite))
                finally:
                    array._mmap.close()
        return dict(status='actual-saved-bank-consumer-pass',checked_cells=checked,catalog=binding(catalog_path),pointer=opened['pointer'],arrays=arrays,coverage=_summary(opened['state']))
    finally:
        coverage._mmap.close()
        digests._mmap.close()
        for orders in sources.values():
            for array in orders.values():
                array._mmap.close()


def main(case_path):
    case = json.loads(Path(case_path).read_bytes())
    out = Path(case['receipts'])
    started, cpu = time.monotonic(), time.process_time()
    context = {k:v for k,v in os.environ.items() if k.startswith('HERDR_')}
    cpu_id = case.get('cpu',10)
    if cpu_id not in range(6,11) or context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != case['pane_id'] or sorted(os.sched_getaffinity(0)) != [cpu_id]:
        raise ValueError('owned actual authorized-CPU/Herdr context required')
    save(out/'runtime.json',dict(context=context,affinity=sorted(os.sched_getaffinity(0)),executable=sys.executable,pythonpath=os.environ.get('PYTHONPATH'),threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},cgroup=Path('/proc/self/cgroup').read_text()))
    results = []
    for index,path in enumerate(case['catalogs']):
        receipt = accumulate(case['bank'],path)
        save(out/(str(index)+'-admission.json'),receipt)
        consumer = consume(case['bank'],path)
        save(out/(str(index)+'-consumer.json'),consumer)
        results.append(dict(admission=binding(out/(str(index)+'-admission.json')),consumer=binding(out/(str(index)+'-consumer.json'))))
        print(json.dumps(dict(event='admitted-and-consumed',catalog=path,component=receipt['component'],added_cells=receipt['added_cells'],coverage=consumer['coverage'])),flush=True)
    final = reopen(case['bank'])
    coverage = _summary(final['state'])
    complete = all(item['complete'] for item in coverage.values())
    publication_path = case.get('canonical_output',case['publication']) if complete else case['publication']
    publication = publish(case['bank'],publication_path)
    final_consumer = None
    if complete:
        from history_banks import OwnedKernels
        reopened = reopen(case['bank'],require_complete=True)
        arrays = []
        with OwnedKernels(publication_path) as kernels:
            for kind,orders in kernels.banks.items():
                for order,array in orders.items():
                    finite = all(bool(np.isfinite(array[i]).all()) for i in range(array.shape[0]))
                    if not finite:
                        raise ValueError('nonfinite complete A62 saved array')
                    arrays.append(dict(kind=kind,order=order,shape=list(array.shape),dtype=str(array.dtype),finite=finite,binding=publication[kind][order]))
        final_consumer = out/'complete-A62-consumer.json'
        save(final_consumer,dict(status='actual-complete-A62-consumer-pass',catalog=binding(publication_path),pointer=reopened['pointer'],coverage=coverage,arrays=arrays,nonclaims=['trajectory','fixed point','physical closure']))
    save(out/'result.json',dict(schema='p254-a82-actual-admission-result-v1',status='complete-admission',case=binding(case_path),results=results,pointer=final['pointer'],coverage=coverage,publication=binding(publication_path),complete_A62_consumer=binding(final_consumer) if final_consumer else None,elapsed_seconds=time.monotonic()-started,native_cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,nonclaims=['trajectory','physical closure']))
    print('A82_ADMISSION_COMPLETE',flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
