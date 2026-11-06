#!/usr/bin/env python3
"""Run unchanged direct coefficient workers for original identity representatives.

No tensor transport or response is calculated here. Each source has one owner,
its own immutable receipt and resumable source state. Existing pilot directories
are consumed in place. A timeout preserves completed sources, not an unfinished
contraction or an unreceipted final archive.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from matched_reference_io import binding, load_case, registration, runtime, save_json

HERE = Path(__file__).resolve().parent
PILOTS = {
    ('loaded-ni31h6', 1, 1): 'matched-loaded-pilot-001',
    ('pristine-ni32', 1, 1): 'matched-pristine-pilot-001',
    ('empty-vacancy-ni31', 1, 1): 'matched-empty-pilot-001',
    ('dispersed-ni32h6', 1, 11): 'matched-dispersed-pilot-001',
}


def representative_jobs(case, rows, output, array_root):
    import numpy as np
    groups = {}
    for _, slot, _ in rows:
        groups.setdefault((slot.ispin, slot.source_ik), []).append(slot)
    jobs = []
    for (spin, ik), slots in sorted(groups.items()):
        originals = [s for s in slots if not s.orbital_conjugation
                     and np.array_equal(s.operation.rotation, np.eye(3))
                     and not np.any(s.operation.fractional_translation)
                     and not np.any(s.reciprocal_shift)]
        if len(originals) != 1:
            raise ValueError(f'{case} {spin}:{ik}: ambiguous original representative')
        slot = originals[0]
        pilot = PILOTS.get((case, spin, ik))
        directory = HERE/pilot if pilot else output/case/f'source-{spin}-{ik}'
        arrays = Path('/var/tmp')/('p254-a28-'+pilot) if pilot else array_root
        argv = [sys.executable, str(HERE/'run_matched_reference_source.py'),
                '--case', case, '--spin', str(spin), '--source-ik', str(ik),
                '--target', ','.join(map(str, slot.mesh_index)),
                '--array-root', str(arrays), '--output', str(directory)]
        receipt = directory/('slot-'+'-'.join(map(str, slot.mesh_index))+'.json')
        jobs.append(dict(spin=spin, source_ik=ik, mesh_index=list(slot.mesh_index),
                         argv=argv, receipt=str(receipt)))
    return jobs


def execute(job, directory, timeout):
    directory.mkdir(parents=True, exist_ok=False)
    save_json(directory/'intent.json', dict(argv=job['argv'], timeout_seconds=timeout))
    started = time.monotonic()
    with (directory/'stdout').open('wb') as out, (directory/'stderr').open('wb') as err:
        try:
            result = subprocess.run(job['argv'], stdout=out, stderr=err, timeout=timeout)
            code, expired = result.returncode, False
        except subprocess.TimeoutExpired:
            code, expired = 124, True
    record = dict(returncode=code, timed_out=expired,
                  elapsed_seconds=time.monotonic()-started,
                  stdout=binding(directory/'stdout'), stderr=binding(directory/'stderr'))
    save_json(directory/'outcome.json', record)
    if code:
        raise RuntimeError(f'direct source worker failed: {directory} ({code})')
    record['receipt'] = binding(job['receipt'])
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--worker-index', type=int, default=0)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--wall-seconds', type=int, required=True)
    parser.add_argument('--source-seconds', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--array-root', type=Path, required=True)
    args = parser.parse_args()
    runtime(args.array_root)
    if not (0 <= args.worker_index < args.workers and
            0 < args.source_seconds <= args.wall_seconds):
        parser.error('invalid worker or bounded time allocation')
    started = time.monotonic()
    reg_binding, reg = registration()
    _, _, _, rows, _ = load_case(args.case, reg, background=False)
    jobs = representative_jobs(args.case, rows, args.output, args.array_root)
    selected = jobs[args.worker_index::args.workers]
    worker = args.output/args.case/f'worker-{args.worker_index}'
    worker.mkdir(parents=True, exist_ok=False)
    save_json(worker/'plan.json', dict(source=binding(__file__), registration=reg_binding,
              jobs=selected, complete_case=False, cpu_affinity=sorted(os.sched_getaffinity(0)),
              wall_seconds=args.wall_seconds, source_seconds=args.source_seconds))
    print(json.dumps(dict(event='prepared', case=args.case, worker=args.worker_index,
                          sources=len(selected))), flush=True)
    completed = []
    status = 'complete-representatives'
    try:
        for index, job in enumerate(selected):
            if args.wall_seconds-(time.monotonic()-started) < args.source_seconds:
                status = 'bounded-stop'
                break
            outcome = execute(job, worker/f'command-{index:03d}', args.source_seconds)
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
