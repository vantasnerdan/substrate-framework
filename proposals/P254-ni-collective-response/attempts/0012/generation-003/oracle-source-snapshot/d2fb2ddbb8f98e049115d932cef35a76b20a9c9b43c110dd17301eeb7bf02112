"""Freeze every source before durable scientific execution; no oracle calls."""
import argparse
import os
from pathlib import Path
import platform
import sys
import time
import traceback
from receipts import save
from source_graph import HERE,closed_graph,own_sources,digest


def prepare(generation):
    directory=HERE/('generation-'+generation)
    directory.mkdir(exist_ok=False)
    try:
        graph=closed_graph()
        sources=own_sources()
        snapshot=directory/'oracle-source-snapshot';snapshot.mkdir()
        for path,sha in sources.items():
            target=snapshot/sha
            if not target.exists():
                with target.open('xb') as stream:
                    stream.write(Path(path).read_bytes());stream.flush();os.fsync(stream.fileno())
        import numpy,scipy,sympy
        dependencies={module.__name__:{'version':module.__version__,'path':module.__file__,
                      'sha256':digest(module.__file__)} for module in (numpy,scipy,sympy)}
        unit='p254-w1-force-sensitivity-0012-g'+generation+'-'+str(time.time_ns())+'.service'
        manifest={'generation':str(directory),'unit':unit,'oracle_sources':sources,
                  'upstream_verified':graph['verified'],'prepared_at_ns':str(time.time_ns()),
                  'scientific_runtime':{'python':sys.executable,'version':sys.version,'platform':platform.platform(),
                    'executable':str(Path('/proc/self/exe').resolve()),'executable_sha256':digest('/proc/self/exe'),
                    'dependencies':dependencies},
                  'registered_cpu':5,'registered_memory_bytes':2147483648,'registered_swap_bytes':0,
                  'registered_blas_threads':1,'preregistration':'banked 464504ff before sensitivity implementation/execution; native outputs predate analysis'}
        save(directory/'source-manifest.json',manifest)
        save(directory/'prepare-receipt.json',{'status':'ready','manifest_sha256':digest(directory/'source-manifest.json'),
             'verified_graph_records':len(graph['verified']),'oracle_sources':sources,'unit':unit})
        print('SENSITIVITY_SOURCES_READY',directory,len(graph['verified']),unit,flush=True)
        return 0
    except BaseException:
        save(directory/'prepare-failure.json',{'status':'failed','traceback':traceback.format_exc(),
             'oracle_sources':own_sources(),'time_ns':str(time.time_ns())})
        print(traceback.format_exc(),flush=True)
        return 2


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--generation',required=True)
    args=parser.parse_args()
    if not args.generation.isdigit():
        raise SystemExit('Numeric generation required')
    raise SystemExit(prepare(args.generation))
