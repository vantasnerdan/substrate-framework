#!/usr/bin/env python3
"""Owned-output and native fork shutdown regressions; no material physics claim."""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
from mixed_material_driver_support import ContinuationBudget, owned_tree
from run_mixed_q_material import watchdog
from saved_case_mean_owner import (
    TERMINAL_BYTES, original_owner, continuation_budget, finish_receipt,
    terminal_signals, reap_owned_workers)
from substrate_framework.streamed_material_state import file_binding


def native_block(ready):
    os.write(ready, b'1')
    os.close(ready)
    # CPU-bound C defers Python handlers, unlike an interruptible sleep.
    hashlib.pbkdf2_hmac('sha256', b'owned-worker', b'P254', 100_000_000)


def ready_worker():
    read, write = os.pipe()
    worker = mp.get_context('fork').Process(target=native_block, args=(write,))
    worker.start()
    os.close(write)
    try:
        if os.read(read, 1) != b'1':
            raise AssertionError('owned worker did not enter its operation')
        deadline = time.monotonic()+3
        while int(Path(f'/proc/{worker.pid}/stat').read_text().split()[13]) == 0:
            if time.monotonic() > deadline:
                raise AssertionError('owned native worker never consumed one CPU tick')
            time.sleep(.005)
    except BaseException:
        if worker.is_alive():
            worker.kill()
        worker.join()
        raise
    finally:
        os.close(read)
    return worker


def check_native_shutdown(directory):
    with terminal_signals():
        worker = ready_worker()
        try:
            budget = ContinuationBudget(directory, started_at=time.time(), max_actions=1,
                max_seconds=.05, max_bytes=10_000_000, max_memory_bytes=1_000_000_000,
                worker_memory_bytes=500_000_000, existing_bytes=0)
            try:
                with watchdog(budget):
                    try:
                        while worker.is_alive():
                            worker.join(timeout=.01)
                    finally:
                        if worker.is_alive():
                            worker.terminate()
                        worker.join()
            except RuntimeError:
                pass
            if worker.is_alive() or worker.exitcode != -signal.SIGTERM:
                raise AssertionError('deadline cleanup returned with a live or non-default-terminated native worker')
        finally:
            reap_owned_workers()
        stranded = ready_worker()
        cleanup = reap_owned_workers()
        if stranded.is_alive() or stranded.exitcode != -signal.SIGKILL:
            raise AssertionError('outer ownership cleanup did not reap its remaining native worker')
    return dict(deadline_worker_exitcode=worker.exitcode,
                outer_cleanup_worker_exitcode=stranded.exitcode, cleanup=cleanup)


def check_output_cap(directory):
    caps = dict(max_actions=1, max_seconds=30., max_output_bytes=1,
                max_memory_bytes=1_000_000_000, worker_memory_bytes=500_000_000)
    tiny = directory/'tiny'
    tiny.mkdir()
    try:
        original_owner(tiny, dict(caps=caps), resume=False, started_unix=time.time())
    except ValueError:
        pass
    else:
        raise AssertionError('one-byte cap admitted owner metadata')
    if sum(owned_tree(tiny).values()) > caps['max_output_bytes'] or (tiny/'owner.json').exists():
        raise AssertionError('refused owner consumed its output cap')
    full = directory/'terminal'
    full.mkdir()
    caps['max_output_bytes'] = 5_000_000
    owner, binding = original_owner(full, dict(caps=caps), resume=False, started_unix=time.time())
    budget = continuation_budget(full, owner)
    finish_receipt(full, 'first', dict(status='partial', scope='finite ownership regression'), budget)
    used = sum(owned_tree(full).values())
    (full/'owned-payload').write_bytes(b'x'*(caps['max_output_bytes']-used-TERMINAL_BYTES+1))
    before = owned_tree(full)
    try:
        continuation_budget(full, owner)
    except RuntimeError:
        pass
    else:
        raise AssertionError('resume admitted more terminal bytes than remained')
    try:
        finish_receipt(full, 'unadmitted', dict(status='partial'), None)
    except RuntimeError:
        pass
    else:
        raise AssertionError('unadmitted terminal was published')
    if owned_tree(full) != before or file_binding(full/'owner.json') != binding:
        raise AssertionError('failed terminal reservation changed immutable owner or output')
    return dict(one_byte_refusal=True, exhausted_resume_stderr_only=True,
                original_owner_unchanged=True, actual_owned_bytes=str(sum(before.values())),
                declared_cap_bytes=str(caps['max_output_bytes']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='p254-mean-lifecycle-') as temporary:
        directory = Path(temporary)
        native = check_native_shutdown(directory)
        output = check_output_cap(directory)
    result = dict(schema='P254-saved-mean-lifecycle-oracle-v1', status='passed',
        scope='Actual owned files, cumulative budget and C-blocked forked processes; no new material mean trajectory',
        native_shutdown=native, output_cap=output,
        sources=[file_binding(path) for path in (
            Path(__file__), HERE/'saved_case_mean_owner.py', HERE/'run_saved_case_mean.py',
            HERE/'run_mixed_q_material.py', HERE/'mixed_material_driver_support.py',
            ROOT/'src/substrate_framework/partitioned_joint_slots.py')])
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
