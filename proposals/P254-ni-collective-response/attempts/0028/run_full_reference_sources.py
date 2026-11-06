#!/usr/bin/env python3
"""Run the fixed independent source list through the existing material units.

No retries, resource discovery or new scientific admission protocol. Each
command gets its exact durable argv/stdout/stderr/outcome before the next one.
An unobserved launcher outcome stops new launches; existing unit caps remain.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import threading
import time


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def publish(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def authenticated(record):
    raw = Path(record['path']).read_bytes()
    if len(raw) != int(record['bytes']) or hashlib.sha256(raw).hexdigest() != record['sha256']:
        raise ValueError('batch input changed: '+record['path'])
    return raw


def worker(cases, began, cap, stopped):
    outcomes = []
    for case in cases:
        if stopped.is_set() or time.monotonic()-began+case['timeout_seconds'] > cap:
            outcomes.append(dict(case=case['key'], status='not-launched'))
            continue
        authenticated(case['grant'])
        prefix = Path(case['command_receipt_prefix'])
        launch = dict(case=case['key'], argv=case['argv'], started_at=stamp())
        publish(str(prefix)+'.start.json', launch)
        print(json.dumps(dict(event='source-start', **launch)), flush=True)
        start = time.monotonic()
        try:
            result = subprocess.run(case['argv'], capture_output=True,
                                    timeout=case['timeout_seconds'], check=False)
            stdout, stderr, code = result.stdout, result.stderr, result.returncode
            state = 'launcher-complete'
        except subprocess.TimeoutExpired as error:
            stdout, stderr, code = error.stdout or b'', error.stderr or b'', None
            state = 'launcher-outcome-unobserved'
        Path(str(prefix)+'.stdout').write_bytes(stdout)
        Path(str(prefix)+'.stderr').write_bytes(stderr)
        outcome = dict(launch, status=state, exit_code=code, finished_at=stamp(),
                       elapsed_seconds=time.monotonic()-start,
                       stdout_path=str(prefix)+'.stdout', stderr_path=str(prefix)+'.stderr')
        publish(str(prefix)+'.json', outcome)
        # The command outcome is durable before reading its scientific output.
        if code == 0:
            summary = json.loads(Path(case['summary_path']).read_bytes())
            outcome['science_status'] = summary['status']
            if case.get('summary_kind') in ('field-moments', 'mixed-field'):
                schema = {'field-moments': 'P254-full-field-moments-material-v1',
                          'mixed-field': 'P254-joint-field-material-v1'}[case['summary_kind']]
                if (summary['schema'] != schema or summary['mode'] != 'source'
                        or [summary['spin'], summary['source_ik']] != case['source_key']):
                    raise ValueError('field checkpoint differs from allocated original source')
                outcome['completed_keys'] = [case['key']]
                outcome['first_useful_block_seconds'] = summary['elapsed_seconds']
            else:
                outcome['completed_keys'] = [item['key'] for item in summary['outputs']]
                outcome['first_useful_block_seconds'] = summary['first_useful_block_seconds']
        if code != 0 or outcome.get('science_status') != 'complete':
            stopped.set()
        outcomes.append(outcome)
        print(json.dumps(dict(event='source-outcome', **outcome)), flush=True)
    return outcomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    args = parser.parse_args()
    batch = json.loads(args.batch.read_bytes())
    for record in batch['implementation_sources']:
        authenticated(record)
    began, stopped = time.monotonic(), threading.Event()
    results = []
    print(json.dumps(dict(event='batch-start', batch=str(args.batch), at=stamp())), flush=True)
    with ThreadPoolExecutor(max_workers=len(batch['workers'])) as pool:
        futures = [pool.submit(worker, cases, began, batch['wall_cap_seconds'], stopped)
                   for cases in batch['workers']]
        for future in as_completed(futures):
            try:
                results.extend(future.result())
            except BaseException as error:
                stopped.set()
                results.append(dict(status='worker-failed', error=repr(error)))
    complete = all(item.get('science_status') == 'complete' for item in results)
    report = dict(status='complete' if complete else 'incomplete', outcomes=results,
                  elapsed_seconds=time.monotonic()-began, finished_at=stamp(),
                  scope='original-source checkpoints, not the full weighted material aggregate')
    publish(batch['result_path'], report)
    print('[dft-done] A28 full-source batch' if complete else '[dft-fail] A28 full-source batch', flush=True)
    return 0 if complete else 1


if __name__ == '__main__':
    raise SystemExit(main())
