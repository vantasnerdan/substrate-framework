"""One durable worker: handshake first, then small verified scientific blocks."""
import argparse
import os
from pathlib import Path
import sys
import time
import traceback
from receipts import save,observe,checkpoint
from source_graph import HERE,load,own_sources,digest


def execute(directory):
    start=time.monotonic();blocks={};ledger=None
    manifest=load(directory/'source-manifest.json')
    try:
        if own_sources()!=manifest['oracle_sources']:
            raise ValueError('Worker code changed since source preparation')
        worker=observe(os.getpid())
        save(directory/'worker-start.json',{'runtime':worker,'manifest_sha256':digest(directory/'source-manifest.json'),
             'oracle_sources':own_sources(),'blas_environment':{k:os.environ.get(k) for k in
             ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')}})
        print('WORKER_POSTEXEC_READY',os.getpid(),worker['start_ticks'],flush=True)
        deadline=time.monotonic()+90
        while not (directory/'worker-release.json').exists():
            if time.monotonic()>deadline:
                raise TimeoutError('No external postexec worker proof/release')
            time.sleep(.1)
        release=load(directory/'worker-release.json')
        if release['pid']!=os.getpid() or release['start_ticks']!=worker['start_ticks']:
            raise ValueError('Release not bound to this worker identity')
        if worker['affinity']!=[5] or worker['live_kernel']['memory.max']!='2147483648' or worker['live_kernel']['memory.swap.max']!='0':
            raise ValueError('Actual worker violates registered caps')
        from substrate_framework.verification import CheckLedger
        ledger=CheckLedger('P254-W1-force-sensitivity-0012')
        import source_graph,check_units,check_native,check_intervals,check_differences
        context={'manifest':manifest}
        sequence=[('sources',source_graph.run),('units',check_units.run),
                  ('native',check_native.run),('intervals',check_intervals.run),
                  ('refinements',check_differences.run_refinements),
                  ('high-start',check_differences.run_high_start)]
        for name,check in sequence:
            print('BLOCK_BEGIN',name,flush=True)
            before=len(ledger.passed); result=check(ledger,context)
            blocks[name]=result
            checkpoint(directory,name,{'result':result,'checks':ledger.passed[before:],
                       'manifest_sha256':digest(directory/'source-manifest.json'),'worker_pid':os.getpid()})
            print('BLOCK_COMPLETE',name,len(ledger.passed)-before,flush=True)
        tally=ledger.finish()
        report={'status':'complete','route_verdict':'established as stated: registered retrospective numerical diagnostics',
                'parent_admission':False,'blocks':blocks,'checks_passed':tally.passed_count,
                'printing_uncertainty':'Exact closed token-rounding boxes only, separate from all unbounded native numerical/physical errors.',
                'limits':context['plan']['limits'],
                'small_ratio_budget':{'source_rounding':'exact rational boxes','conversion':'exact source decimal CODATA2006',
                  'SCF_force_forward_error':'not bounded by printed estimated energy accuracy',
                  'cutoff_k_MV':'registered one-step diagnostics; no extrapolated asymptotic error bound',
                  'branch_identity':'cold-start diagnostic, not proof of same state or global ground state',
                  'zero_mode_spectrum_conditioning_jitter':'not evaluated: no Hessian/eigenvalue or projected force claim',
                  'remaining':'all-four-primary Gamma barrier; mode/q/size/strain/density/isotope/disorder/matched controls; G1-G5 and authorized final merger'}}
        save(directory/'report.json',report)
        save(directory/'receipt.json',{'status':'complete','exit_code':0,'elapsed_seconds':time.monotonic()-start,
             'worker':worker,'final_worker':observe(os.getpid()),'source_manifest_sha256':digest(directory/'source-manifest.json'),
             'report_sha256':digest(directory/'report.json'),'checks':ledger.passed,
             'checkpoint_sha256':digest(directory/'checkpoint.json')})
        print('SENSITIVITY_ORACLE_COMPLETE',directory/'report.json',digest(directory/'report.json'),flush=True)
        return 0
    except BaseException as error:
        failure={'status':'failed','error':repr(error),'traceback':traceback.format_exc(),
                 'elapsed_seconds':time.monotonic()-start,'completed_blocks':list(blocks),
                 'checks':ledger.passed if ledger else [],'source_manifest_sha256':digest(directory/'source-manifest.json')}
        save(directory/'failure.json',failure)
        print(failure['traceback'],flush=True)
        print('SENSITIVITY_ORACLE_FAILED',directory,flush=True)
        return 2


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--generation',required=True,type=Path)
    arguments=parser.parse_args()
    if arguments.generation.parent!=HERE or not arguments.generation.name.startswith('generation-'):
        raise SystemExit('Generation must be owned immediate child')
    raise SystemExit(execute(arguments.generation))
