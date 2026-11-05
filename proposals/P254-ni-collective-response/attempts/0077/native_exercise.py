"""One real A71/A72 admission and saved-bank consumer, not a physics replay."""
from pathlib import Path
import json
import os
import resource
import time
import numpy as np
from component_accumulator import accumulate, publish, reopen, HERE
from history_banks import binding, authenticate, save, ORDERS, OwnedKernels


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=('001', '002'), default='001')
    args = parser.parse_args()
    started = time.monotonic()
    cpu_start = time.process_time()
    output = HERE / ('exercise-' + args.case)
    output.mkdir(exist_ok=False)
    root = output / 'bank'
    context = {k: v for k, v in os.environ.items() if k.startswith('HERDR_')}
    if context.get('HERDR_ENV') != '1' or not all(context.get(k) for k in ('HERDR_WORKSPACE_ID', 'HERDR_TAB_ID', 'HERDR_PANE_ID')):
        raise ValueError('actual inherited owned Herdr context required')
    save(output / 'startup.json', dict(pid=os.getpid(), context=context, affinity=sorted(os.sched_getaffinity(0)),
        argv=__import__('sys').argv, python=__import__('sys').executable,
        thread_environment={k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
        cgroup=Path('/proc/self/cgroup').read_text()))
    prereg = json.loads((HERE / 'preregistration.json').read_bytes())
    receipts = []
    for name in ('EI_catalog', 'cubic_catalog'):
        descriptor = prereg[name]
        authenticate(descriptor)
        receipt = accumulate(root, descriptor['path'])
        save(output / (name + '-admission.json'), receipt)
        reopened = reopen(root)
        if reopened['pointer'] != receipt['pointer']:
            raise ValueError('reopened pointer differs from committed admission')
        receipts.append(receipt)
        print(json.dumps(dict(event='admitted', component=receipt['component'], cells=receipt['added_cells'])), flush=True)
    state = reopen(root)['state']
    consumer = []
    # Every saved element of every contour order is consumed and compared, one
    # source/bank mapping at a time; neither an all-component RAM bank nor mocks.
    for name, component, kind, catalog_key in (('EI_catalog','mixed_ei','sigma','sigma'),
                                               ('cubic_catalog','cubic','cubic','arrays')):
        catalog = json.loads(Path(prereg[name]['path']).read_bytes())
        family = state['families'][component]
        per_order = {}
        for order in ORDERS:
            source = np.load(authenticate(catalog[catalog_key][order]), mmap_mode='r', allow_pickle=False)
            bank = np.load(authenticate(family['arrays'][kind][order]), mmap_mode='r', allow_pickle=False)
            try:
                difference = 0.0
                nonzero = 0
                for a in range(31):
                    for b in range(31):
                        if not np.array_equal(bank[a,b], source[a,b]):
                            raise ValueError('first real q1 accumulated numeric values differ')
                        nonzero += int(np.count_nonzero(bank[a,b]))
                        difference = max(difference, float(np.max(np.abs(bank[a,b]-source[a,b]))))
                per_order[order] = dict(max_abs_difference=difference, nonzero_entries=nonzero,
                                       compared_entries=int(bank.size))
            finally:
                source._mmap.close()
                bank._mmap.close()
        consumer.append(dict(component=component, orders=per_order))
    before = binding(root / 'current.json')
    duplicates = []
    for name in ('EI_catalog', 'cubic_catalog'):
        receipt = accumulate(root, prereg[name]['path'])
        if receipt['status'] != 'idempotent' or receipt['added_cells'] != 0:
            raise ValueError('identical paid import was not idempotent')
        if binding(root / 'current.json') != before:
            raise ValueError('duplicate import changed atomic pointer bytes')
        duplicates.append(receipt)
    final = reopen(root)
    previous = json.loads(authenticate(final['pointer']['previous']).read_bytes())
    if previous['generation'] != 1 or previous['families']['mixed_ei'] != final['state']['families']['mixed_ei']:
        raise ValueError('previous state or unchanged-family reference was not retained')
    published_path = output / 'partial-publication.json'
    publication = publish(root, published_path)
    if publication['full_solve_eligible'] or publication['status'] != 'partial' or 'sigma' in publication or 'pi' in publication:
        raise ValueError('partial storage was upgraded to final six-array eligibility')
    rejected = {}
    try:
        reopen(root, require_complete=True)
    except ValueError as error:
        rejected['accumulator'] = str(error)
    else:
        raise ValueError('complete solve admitted missing EE/joint/restq')
    try:
        with OwnedKernels(published_path):
            pass
    except ValueError as error:
        rejected['A62_OwnedKernels'] = str(error)
    else:
        raise ValueError('partial receipt admitted by A62 complete solver')
    findings = dict(schema='p254-a77-real-saved-path-exercise-v1', status='complete-honestly-partial-bank',
        input_catalogs={name: prereg[name] for name in ('EI_catalog','cubic_catalog')},
        input_dataset=state['input_dataset'], original_B_identity=state['original_B_identity'],
        source_code=[binding(HERE / name) for name in ('component_accumulator.py','native_exercise.py','launch_native.py')],
        admissions=receipts, saved_path_consumer=consumer, duplicates=duplicates,
        pointer_reopen=final['pointer'], previous_generation=previous['generation'],
        unchanged_EI_family_shared=True, retained_state_count=len(list((root/'states').iterdir())),
        retained_family_count=len(list((root/'objects').iterdir())),
        partial_publication=binding(published_path), coverage=publication['coverage'], complete_solve_rejected=rejected,
        elapsed_seconds=time.monotonic()-started, native_cpu_seconds=time.process_time()-cpu_start,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        nonclaims=['EE admission','joint admission','remaining63 EI/cubic momenta','all64/all961 full functional','physical solve'])
    save(output / 'findings.json', findings)
    print(json.dumps(dict(event='exercise-terminal',status=findings['status'],elapsed_seconds=findings['elapsed_seconds'],
        native_cpu_seconds=findings['native_cpu_seconds'],findings=binding(output/'findings.json'))),flush=True)


if __name__ == '__main__':
    main()
