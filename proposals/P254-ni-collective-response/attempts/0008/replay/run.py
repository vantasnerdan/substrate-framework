#!/usr/bin/env python3
"""Replay the frozen complete audit into a fresh output directory; originals stay intact."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
AUDIT=HERE.parent
sys.path.insert(0,str(AUDIT))
from run_oracle import context,digest,runtime,save


def run(output):
    started=time.monotonic()
    output.mkdir(parents=True,exist_ok=False)
    receipt={'status':'running','blocks':{},'runtime_start':runtime()}
    try:
        frozen=json.loads((AUDIT/'source-manifest.json').read_text())
        if any(digest(p)!=h for p,h in frozen['files'].items()):
            raise ValueError('Original source manifest mismatch')
        reg=json.loads((HERE/'registration.json').read_text())
        if digest(__file__)!=reg['replay_source_sha256']:
            raise ValueError('Replay source registration mismatch')
        actual=receipt['runtime_start']
        if actual['affinity']!=[5] or actual['cgroup_observed']['memory.max']!='2147483648' or any(actual['environment'][k]!='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')):
            raise ValueError('CPU5/2GiB/BLAS1 ownership violated')
        import numpy,scipy,sympy
        c=context()
        sys.path.append(str(Path(c['supplement']['ase_root']).parent))
        from substrate_framework.verification import CheckLedger
        import check_sources,check_units,check_geometry,check_forces,check_observables
        ledger=CheckLedger('P254-primary-frozen-replay')
        print('PRIMARY_REPLAY_STARTED',output,flush=True)
        for name,module in [('sources',check_sources),('units',check_units),('geometry',check_geometry),('forces',check_forces),('observables',check_observables)]:
            offset=len(ledger.passed); result=module.run(ledger,c)
            receipt['blocks'][name]=result
            save(output/f'{name}.json',{'source_manifest_sha256':digest(AUDIT/'source-manifest.json'),'checks':ledger.passed[offset:],'result':result})
        original=json.loads((AUDIT/'receipt.json').read_text())
        ledger.check('all independently replayed numerical blocks exactly reproduce first receipt',all(value==original['blocks'][name]['result'] for name,value in receipt['blocks'].items()))
        # Exercise omitted-claim failure channels against full physical-identity predicate.
        final=c['final']; identity=lambda g:g==final
        for key,replacement in [('cell_angstrom',(numpy.eye(3)*7.1).tolist()),('pbc',[False]*3),('symbols',final['symbols'][:-1]),('constraint_count',1)]:
            ledger.mutation_sensitive('physical identity mutation '+key,identity,final,[{**final,key:replacement}])
        # Alter the raw native text, reparse it, and prove published force comparison rejects it.
        from native_parse import parse_stdout
        text=Path(c['scf']['data']['artifacts']['pw.stdout']['path']).read_text()
        mutated=text.replace('force =     0.00000000    0.00012328    0.00012328','force =     0.00000000    0.00022328    0.00012328',1)
        raw=parse_stdout(mutated)['forces']*c['units']['force_factor']
        ledger.check('native stdout force mutation is actually changed and publication rejects it',mutated!=text and numpy.max(abs(raw-numpy.array(c['published']['data']['forces_ev_per_angstrom'])))>1e-3)
        ledger.finish()
        receipt.update(status='complete',checks=ledger.passed,check_count=len(ledger.passed),source_manifest_sha256=digest(AUDIT/'source-manifest.json'),registration_sha256=digest(HERE/'registration.json'),runtime_end=runtime(),elapsed_seconds=time.monotonic()-started)
        save(output/'receipt.json',receipt)
        print('PRIMARY_REPLAY_COMPLETE',digest(output/'receipt.json'),flush=True)
        return 0
    except BaseException:
        receipt.update(status='failed',traceback=traceback.format_exc(),runtime_end=runtime())
        save(output/'failure.json',receipt);traceback.print_exc()
        print('PRIMARY_REPLAY_TERMINAL_FAILURE',flush=True)
        return 2


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',required=True,type=Path)
    raise SystemExit(run(parser.parse_args().output_dir))
