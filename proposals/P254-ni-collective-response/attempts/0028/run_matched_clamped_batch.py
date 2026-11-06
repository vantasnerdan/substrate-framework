#!/usr/bin/env python3
"""Run the fixed clamped-control plan through the existing command executor.

Per-case stdout, stderr and exit receipts remain durable. Independent cases
continue after an observed failure; a comparison starts only after both inputs
complete. No retries, new physics, resource discovery or result substitution.
"""
import argparse
import json
import time
from pathlib import Path

from material_run_contract import binding, bound_json, read_bound_bytes
from matched_reference_io import save_json
from run_matched_representative_batch import execute


def load_result(record):
    """Authenticate scientific output, not a source document (16MiB limit)."""
    raw, _ = read_bound_bytes(record, max_bytes=256 << 20)
    return json.loads(raw)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--wall-seconds', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan_binding = binding(args.plan)
    plan = bound_json(plan_binding)
    for record in plan['implementation_sources']:
        if binding(record['path']) != record:
            raise ValueError('batch implementation changed: '+record['path'])
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    completed, failed, held = {}, [], []
    print(json.dumps(dict(event='batch-start', plan=plan_binding)), flush=True)
    try:
        for index, job in enumerate(plan['jobs']):
            if any(key not in completed for key in job.get('requires', [])):
                held.append(dict(key=job['key'], reason='incomplete dependency'))
                continue
            timeout = job['timeout_seconds']
            if time.monotonic()-started+timeout > args.wall_seconds:
                held.append(dict(key=job['key'], reason='remaining cumulative wall budget'))
                continue
            print(json.dumps(dict(event='case-start', key=job['key'])), flush=True)
            try:
                outcome = execute(job, args.output/f'command-{index:03d}', timeout)
                result = load_result(outcome['receipt'])
                if result['status'] != 'complete':
                    raise ValueError('producer did not complete: '+result['status'])
                completed[job['key']] = outcome
                print(json.dumps(dict(event='case-complete', key=job['key'],
                                      receipt=outcome['receipt'])), flush=True)
            except Exception as error:
                failure = dict(key=job['key'], error=repr(error))
                failed.append(failure)
                save_json(args.output/f'failure-{index:03d}.json', failure)
                print(json.dumps(dict(event='case-failed', **failure)), flush=True)
    finally:
        complete = len(completed) == len(plan['jobs']) and not failed and not held
        save_json(args.output/'result.json', dict(
            status='complete' if complete else 'incomplete', plan=plan_binding,
            completed=completed, failed=failed, held=held,
            elapsed_seconds=time.monotonic()-started,
            scope='clamped controls and mass-only comparisons; coupled E/D/C remain open'))
        print('[dft-done] A28 matched clamped batch' if complete else
              '[dft-fail] A28 matched clamped batch', flush=True)
    return 0 if complete else 1


if __name__ == '__main__':
    raise SystemExit(main())
