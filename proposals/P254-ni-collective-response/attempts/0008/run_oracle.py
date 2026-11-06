#!/usr/bin/env python3
"""Durable generation-owned primary audit, source registration and checkpoints."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
GENERATION='primary-boundary-0008-v1'


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def manifest():
    paths=list(HERE.glob('*.py'))+[HERE/'registration.json',HERE.parents[1]/'W1-loaded-primary-audit-plan.json',ROOT/'src/substrate_framework/material_primary_boundary.py',ROOT/'src/substrate_framework/verification.py']
    return {str(p):digest(p) for p in sorted(paths)}


def runtime():
    cg=Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    root=Path('/sys/fs/cgroup')/cg.lstrip('/')
    observed={name:(root/name).read_text().strip() for name in ('memory.current','memory.peak','memory.max','memory.events','cpu.stat') if (root/name).exists()}
    return {'pid':os.getpid(),'affinity':sorted(os.sched_getaffinity(0)),
            'environment':{k:os.getenv(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','PYTHONPATH')},
            'cgroup':cg,'cgroup_observed':observed,'proc_status':Path('/proc/self/status').read_text(),
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def context():
    plan=json.loads((HERE.parents[1]/'W1-loaded-primary-audit-plan.json').read_text())
    supplement=json.loads((HERE/'registration.json').read_text())
    published_path=next(p for p in plan['input_sources_sha256'] if p.endswith('loaded-ni31h6/result.json'))
    published=json.loads(Path(published_path).read_text())
    read_json=lambda p:json.loads(Path(p).read_text())
    spec=read_json(published['data']['artifacts']['force_spec']['path'])
    scf=read_json(published['data']['artifacts']['force_result']['path'])
    attempt=Path(scf['data']['attempt_directory'])
    return {'plan':plan,'supplement':supplement,'published':published,'scf':scf,'spec':spec,
            'relax':read_json(Path(published_path).parent.parent/'loaded-ni31h6.json'),
            'optimizer':read_json(Path(published_path).parent/'optimizer.json'),
            'settings':read_json(next(p for p in plan['input_sources_sha256'] if p.endswith('w1-settings.json'))),
            'attempt_spec':read_json(attempt/'spec.json'),'command':read_json(attempt/'command.json')}


def execute():
    start=time.monotonic(); run_id=f'{GENERATION}-{time.time_ns()}'
    report={'schema_version':1,'generation':GENERATION,'run_id':run_id,'status':'running','blocks':{}}
    try:
        registered=json.loads((HERE/'source-manifest.json').read_text())
        if registered['files']!=manifest():
            raise ValueError('Source generation changed; register a distinct generation')
        report['source_manifest_sha256']=digest(HERE/'source-manifest.json')
        report['runtime_start']=runtime(); save(HERE/f'{run_id}-start.json',report['runtime_start'])
        print('PRIMARY_BOUNDARY_STARTED',run_id,flush=True)
        proof=report['runtime_start']
        if proof['affinity']!=[5] or proof['cgroup_observed']['memory.max']!='2147483648' or any(proof['environment'][k]!='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')):
            raise ValueError('Actual CPU5/2GiB/BLAS1 owner contract violated')
        # Only ASE cross-checks resolve here; SciPy/SymPy use the canonical oracle environment.
        import numpy, scipy, sympy
        sys.path.append(str(Path(json.loads((HERE/'registration.json').read_text())['ase_root']).parent))
        from substrate_framework.verification import CheckLedger
        import check_sources,check_units,check_geometry,check_forces,check_observables
        ledger=CheckLedger('P254-loaded-primary-boundary')
        c=context()
        report['dependencies']={'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,'sympy':sympy.__version__}
        for name,module in [('sources',check_sources),('units',check_units),('geometry',check_geometry),('forces',check_forces),('observables',check_observables)]:
            offset=len(ledger.passed)
            result=module.run(ledger,c)
            checkpoint={'generation':GENERATION,'source_manifest_sha256':report['source_manifest_sha256'],'result':result,'checks':ledger.passed[offset:]}
            path=HERE/f'{run_id}-{name}-checkpoint.json'; save(path,checkpoint)
            report['blocks'][name]={'checkpoint':str(path),'sha256':digest(path),'result':result,'checks':checkpoint['checks']}
            print('PRIMARY_BOUNDARY_CHECKPOINT',run_id,name,digest(path),flush=True)
        ledger.finish()
        if registered['files']!=manifest():
            raise ValueError('Source generation changed during execution')
        report.update(status='complete',route_verdict='established_as_stated_primary_resolution_publication_and_raw_stationarity',checks=ledger.passed,check_count=len(ledger.passed),limitations=c['plan']['limits'])
        report['runtime_end']=runtime(); report['elapsed_seconds']=time.monotonic()-start
        save(HERE/'receipt.json',report)
        print('PRIMARY_BOUNDARY_COMPLETE',run_id,digest(HERE/'receipt.json'),flush=True)
        return 0
    except BaseException:
        report.update(status='failed',traceback=traceback.format_exc(),runtime_end=runtime(),elapsed_seconds=time.monotonic()-start)
        save(HERE/f'{run_id}-failure.json',report)
        traceback.print_exc()
        print('PRIMARY_BOUNDARY_TERMINAL_FAILURE',run_id,flush=True)
        return 2


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--prepare',action='store_true'); args=parser.parse_args()
    if args.prepare:
        save(HERE/'source-manifest.json',{'schema_version':1,'generation':GENERATION,'registered_ns':time.time_ns(),'files':manifest()})
        print('PRIMARY_BOUNDARY_REGISTERED',GENERATION,digest(HERE/'source-manifest.json'),flush=True)
    else:
        raise SystemExit(execute())
