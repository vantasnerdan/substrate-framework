"""Execute one explicitly allocated packing or real component-production case."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/'0063'), str(HERE.parent/'0076')]
from history_storage import authenticate, binding
from production_transport import atomic_json


def run(case_path):
    started = time.monotonic()
    case_path = Path(case_path).resolve()
    case = json.loads(case_path.read_bytes())
    receipt = Path(case['receipts'])
    group = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
        memory_max=(cg/'memory.max').read_text().strip(),
        swap_max=(cg/'memory.swap.max').read_text().strip(),
        environment={key:value for key,value in os.environ.items() if key.startswith('HERDR_') or key in (
            'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
        executable=sys.executable,pythonpath=os.environ.get('PYTHONPATH'))
    atomic_json(receipt/'runtime.json',runtime)
    if (runtime['affinity'] != [case['cpu']] or runtime['memory_max'] != str(case['memory_bytes']) or
            runtime['swap_max'] != '0' or runtime['environment'].get('HERDR_ENV') != '1' or
            runtime['environment'].get('HERDR_PANE_ID') != case['pane_id'] or
            any(runtime['environment'].get(key) != '1' for key in (
                'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'))):
        raise RuntimeError('allocated actual Native runtime differs')
    for source in case['source_bindings']:
        authenticate(source)
    result = dict(schema='p254-a80-real-production-result-v1',case=binding(case_path),
                  operation=case['operation'],status='running',runtime=runtime)
    atomic_json(receipt/'result.json',result)
    try:
        if case['operation'] == 'pack-q1':
            from pack_paid_ee import pack
            catalog = pack(case['output'])
        elif case['operation'] == 'component':
            from component_emitter import produce_component
            catalog = produce_component(case['component'],case['q'],case['output'],
                times=case.get('times'),owned_pairs=case.get('owned_pairs'))
            from component_consumer import consume_component
            consumer = catalog.parent/'consumer.json'
            atomic_json(consumer,consume_component(catalog))
            result['consumer'] = binding(consumer)
        else:
            raise ValueError('unknown explicitly allocated operation')
        result.update(status='complete-owned-production',catalog=binding(catalog))
    except BaseException as error:
        result.update(status='failed',failure=dict(type=type(error).__name__,message=str(error),
                                                  traceback=traceback.format_exc()))
        raise
    finally:
        result.update(elapsed_seconds=time.monotonic()-started,
            native_cpu_seconds=time.process_time(),
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            full_mesh_claim=False)
        atomic_json(receipt/'result.json',result)
        print('A80_NATIVE_RESULT '+json.dumps(result),flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: run_production.py ACTUAL_BOUND_CASE.json')
    run(sys.argv[1])
