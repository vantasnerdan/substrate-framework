"""Produce real all-sector Gamma coefficients, then the nonstationary reduced mean.

Source chunks are distinct work, not retries. Completed immutable cache items
survive interruption; the original production deadline never resets on resume.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import perf_counter
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0028')]
from reduced_projection_cache import atomic_json
from matched_reference_io import runtime
from substrate_framework.streamed_material_state import file_binding


def command(directory, argv, timeout):
    directory.mkdir(parents=True, exist_ok=False)
    intent = dict(argv=argv, cwd=str(ROOT), actual_started_utc=datetime.now(timezone.utc).isoformat(),
                  timeout_seconds=timeout)
    atomic_json(directory/'intent.json', intent)
    started = perf_counter()
    code = None
    failure = None
    process = None
    with (directory/'stdout').open('wb') as out, (directory/'stderr').open('wb') as err:
        try:
            process = subprocess.Popen(argv, cwd=ROOT, stdout=out, stderr=err)
            code = process.wait(timeout=timeout)
        except BaseException as exc:
            failure = dict(kind=type(exc).__name__, message=str(exc))
            if process is not None:
                process.kill()
                code = process.wait()
    outcome = dict(returncode=code, failure=failure, elapsed_seconds=perf_counter()-started,
                   stdout=file_binding(directory/'stdout'), stderr=file_binding(directory/'stderr'))
    atomic_json(directory/'outcome.json', outcome)
    if failure is not None or code != 0:
        raise RuntimeError('scientific command failed; actual outcome '+str(directory/'outcome.json'))
    return outcome


def run(args):
    runtime(args.output)
    now = datetime.now(timezone.utc)
    config = dict(case=args.case, boundary=file_binding(args.boundary), cache=str(args.cache.resolve()),
                  geometry=file_binding(args.geometry_archive), seconds=args.seconds,
                  bands=args.bands, features=args.features, times=args.times, chunk=args.chunk)
    contract_path = args.output/'production-contract.json'
    if contract_path.exists():
        if not args.resume:
            raise ValueError('existing production requires explicit --resume')
        contract = json.loads(contract_path.read_text())
        if contract['configuration'] != config:
            raise ValueError('original production identity/deadline differs')
    else:
        contract = dict(configuration=config, original_started_utc=now.isoformat(),
                        deadline_unix=now.timestamp()+args.seconds,
                        no_automatic_retry=True,
                        source_files=[file_binding(HERE/name) for name in (
                            'run_reduced_gamma_campaign.py', 'reduced_material_projection.py',
                            'run_reduced_material_dynamics.py')])
        atomic_json(contract_path, contract)
    run_dir = args.output/'attempts'/uuid4().hex
    run_dir.mkdir(parents=True)
    result = dict(schema='P254-reduced-Gamma-production-v1', status='running',
                  actual_invocation_id=os.environ.get('INVOCATION_ID'),
                  original_contract=file_binding(contract_path), commands=[],
                  scientific_output='all128 explicit-sector Gamma inputs and actual reduced mean at requested times',
                  scope='No full64-q spatial response, finite-dose result, material accuracy or campaign completion claim')
    result_path = run_dir/'result.json'
    atomic_json(result_path, result)
    print(json.dumps(dict(event='reduced-Gamma-start', case=args.case, run=str(result_path))), flush=True)
    remaining = lambda: contract['deadline_unix']-datetime.now(timezone.utc).timestamp()
    try:
        for index, first in enumerate(range(0, 128, args.chunk)):
            # Final export additionally assembles both all-sector passive splits.
            # It shares the original case deadline, without a renewed budget.
            cap = remaining()-5. if first+args.chunk >= 128 else 600.
            budget = min(cap, remaining()-5.)
            if budget <= 0:
                raise TimeoutError('original production deadline exhausted; completed scientific cache preserved')
            last = min(127, first+args.chunk-1)
            argv = [sys.executable, str(HERE/'reduced_material_projection.py'), 'export',
                    '--case', args.case, '--boundary', str(args.boundary), '--cache', str(args.cache),
                    '--slots', f'{first}-{last}', '--q', '0', '--bands', str(args.bands),
                    '--features', str(args.features), '--seconds', str(budget), '--maxrss-mib', '3072']
            if index == 0:
                argv.append('--kernels')
            outcome = command(run_dir/('source-chunk-%02d' % index), argv,
                              min(budget+20., remaining()))
            result['commands'].append(outcome)
            manifest = json.loads((args.cache/'manifest.json').read_text())
            sources = sum(('edge/00/%03d' % s) in manifest['items'] for s in range(128))
            print(json.dumps(dict(event='source-chunk-complete', chunk=index, case=args.case,
                                  completed_Gamma_sources=sources, manifest=str(args.cache/'manifest.json'),
                                  actual_chunk_seconds=outcome['elapsed_seconds'])), flush=True)
            atomic_json(result_path, result)
        manifest = json.loads((args.cache/'manifest.json').read_text())
        if manifest['status'] not in ('gamma-complete', 'complete'):
            raise RuntimeError('all requested commands finished without a complete Gamma input bank')
        print(json.dumps(dict(event='Gamma-inputs-complete', case=args.case,
                              manifest=str(args.cache/'manifest.json'))), flush=True)
        if remaining() <= 0:
            raise TimeoutError('original deadline reached before reduced mean; inputs retained')
        baseline_output = args.output/'baseline'
        argv = [sys.executable, str(HERE/'run_reduced_material_dynamics.py'),
                '--manifest', str(args.cache/'manifest.json'), '--geometry-archive', str(args.geometry_archive),
                '--output', str(baseline_output), '--bands', '8', '--features', '7',
                '--times', args.times, '--baseline-only']
        if (baseline_output/'result.json').exists():
            argv.append('--resume')
        outcome = command(run_dir/'actual-reduced-mean', argv, remaining())
        result['commands'].append(outcome)
        result['coefficient_manifest'] = file_binding(args.cache/'manifest.json')
        result['actual_reduced_mean'] = file_binding(baseline_output/'result.json')
        result['status'] = 'complete'
    except BaseException as exc:
        result['status'] = 'failed'
        result['failure'] = dict(kind=type(exc).__name__, message=str(exc))
        raise
    finally:
        result['actual_finished_utc'] = datetime.now(timezone.utc).isoformat()
        atomic_json(result_path, result)
        print(('[dft-done]' if result['status'] == 'complete' else '[dft-fail]')+
              ' reduced-Gamma '+args.case+' '+str(result_path), flush=True)


def interrupted(signum, frame):
    del frame
    raise InterruptedError('production interrupted by signal '+str(signum))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--boundary', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--geometry-archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=7200.)
    parser.add_argument('--bands', type=int, choices=(8, 12), default=12)
    parser.add_argument('--features', type=int, choices=(7, 27), default=27)
    parser.add_argument('--times', default='0,1,10,30')
    parser.add_argument('--chunk', type=int, default=8)
    parser.add_argument('--resume', action='store_true')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    args = parser.parse_args()
    if not 1 <= args.chunk <= 128 or args.seconds <= 0:
        raise ValueError('positive original budget and source-chunk size required')
    run(args)


if __name__ == '__main__':
    main()
