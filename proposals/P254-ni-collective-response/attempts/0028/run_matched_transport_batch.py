#!/usr/bin/env python3
"""Execute existing initial-coefficient transport for disjoint source orbits.

No contraction mathematics is implemented here. Completed direct representative
receipts are the only inputs; each child publishes a complete source orbit.
Completed child results survive a batch timeout. An interrupted source is not
resumable by this driver and requires a distinct output directory.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

from matched_reference_io import binding, load_case, registration, runtime, save_json
from matched_transport_io import transport_registration
from run_matched_representative_batch import PILOTS, execute

HERE = Path(__file__).resolve().parent


def jobs_for(args, rows):
    jobs = []
    for spin, ik in sorted({(s.ispin, s.source_ik) for _, s, _ in rows}):
        pilot = PILOTS.get((args.case, spin, ik))
        source = HERE/pilot if pilot else args.representative_root/args.case/f'source-{spin}-{ik}'
        output = args.output/args.case/f'source-{spin}-{ik}'
        argv = [sys.executable, str(HERE/'run_matched_transport_material.py'),
                '--case', args.case, '--spin', str(spin), '--source-ik', str(ik),
                '--target', 'all', '--input-dir', str(source), '--oracle', str(args.oracle),
                '--accept-declared-geometry', '--array-root', str(args.array_root),
                '--output', str(output)]
        jobs.append(dict(spin=spin, source_ik=ik, argv=argv,
                         direct_input=str(source), receipt=str(output/'result.json')))
    return jobs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', required=True)
    p.add_argument('--worker-index', type=int, default=0)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--wall-seconds', type=int, required=True)
    p.add_argument('--source-seconds', type=int, required=True)
    p.add_argument('--representative-root', type=Path, required=True)
    p.add_argument('--oracle', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--array-root', type=Path, required=True)
    args = p.parse_args()
    if not 0 <= args.worker_index < args.workers or not 0 < args.source_seconds <= args.wall_seconds:
        p.error('invalid disjoint worker or time allocation')
    runtime(args.array_root)
    started = time.monotonic()
    direct_binding, direct = registration()
    transport_binding, _ = transport_registration()
    _, _, _, rows, _ = load_case(args.case, direct, background=False)
    selected = jobs_for(args, rows)[args.worker_index::args.workers]
    worker = args.output/args.case/f'worker-{args.worker_index}'
    worker.mkdir(parents=True, exist_ok=False)
    save_json(worker/'plan.json', dict(source=binding(__file__),
        executor=binding(HERE/'run_matched_representative_batch.py'),
        direct_registration=direct_binding, transport_registration=transport_binding,
        oracle=binding(args.oracle), jobs=selected, cpu_affinity=sorted(os.sched_getaffinity(0)),
        wall_seconds=args.wall_seconds, source_seconds=args.source_seconds,
        complete_case=False, positive_time_symmetry=False))
    print(json.dumps(dict(event='prepared', case=args.case, worker=args.worker_index,
                          sources=len(selected))), flush=True)
    completed, status = [], 'complete-source-orbits'
    try:
        for index, job in enumerate(selected):
            if args.wall_seconds-(time.monotonic()-started) < args.source_seconds:
                status = 'bounded-stop'
                break
            outcome = execute(job, worker/f'command-{index:03d}', args.source_seconds)
            report = json.loads(Path(job['receipt']).read_bytes())
            if (report['status'] != 'complete-selected-source-orbit'
                    or not report['complete_source_orbit'] or report['source_coefficients_recomputed']):
                raise ValueError('child did not publish the complete original source orbit')
            completed.append(outcome)
            print(json.dumps(dict(event='source-complete', spin=job['spin'],
                source_ik=job['source_ik'], receipt=outcome['receipt'])), flush=True)
    except Exception:
        status = 'failed'
        raise
    finally:
        save_json(worker/'result.json', dict(status=status, completed=completed,
            planned_sources=len(selected), complete_case=False, response_comparison=False,
            elapsed_seconds=time.monotonic()-started))
        print(json.dumps(dict(event='terminal', status=status, completed=len(completed))), flush=True)


if __name__ == '__main__':
    main()
