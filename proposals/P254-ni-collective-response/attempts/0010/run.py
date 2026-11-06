#!/usr/bin/env python3
"""Run only the preregistered explicit energy-gauge oracle; never DFT."""
import argparse
import json
from pathlib import Path

import check_algebra
import check_consumer
import check_projection
from oracle_evidence import atomic_json, digest, runtime, source_guard

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest_path = HERE / 'source-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest_sha = digest(manifest_path)
    source_guard(manifest['files'])
    args.output.mkdir(parents=True, exist_ok=False)
    observed_start = runtime()
    if observed_start['affinity'] != [5] or observed_start['cgroup_observed']['memory.max'] != '2147483648':
        raise ValueError('Actual CPU/memory reservation differs from preregistration')
    if observed_start['cgroup_observed']['memory.swap.max'] != '0':
        raise ValueError('Actual swap cap differs from preregistration')
    if any(observed_start['environment'][name] != '1' for name in
           ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')):
        raise ValueError('Actual BLAS thread reservation differs from preregistration')
    atomic_json(args.output / 'runtime-start.json', observed_start)
    ledger, blocks = [], {}
    tasks = [('operator_algebra', lambda: check_algebra.run(ledger)),
             ('complex_projection', lambda: check_projection.run(ledger)),
             ('production_consumer', lambda: check_consumer.run(ledger, HERE / 'inputs/production-gauge-consumer.py'))]
    for name, work in tasks:
        blocks[name] = work()
        atomic_json(args.output / (name + '.json'), {'source_manifest_sha256': manifest_sha, 'result': blocks[name]})
        print('FERMI_GAUGE_BLOCK_COMPLETE', name, flush=True)
    source_guard(manifest['files'])
    if digest(manifest_path) != manifest_sha:
        raise ValueError('Source manifest changed during execution')
    observed_end = runtime()
    events = dict(line.split() for line in observed_end['cgroup_observed']['memory.events'].splitlines())
    if any(int(events.get(name, 0)) for name in ('max', 'oom', 'oom_kill', 'oom_group_kill')):
        raise ValueError('Kernel observed memory cap/OOM event during oracle')
    receipt = {'schema_version': 1, 'status': 'complete', 'source_manifest_sha256': manifest_sha,
               'checks': ledger, 'check_count': len(ledger), 'blocks': blocks,
               'consumer': manifest['production_consumer'], 'runtime_start': observed_start, 'runtime_end': observed_end,
               'scope': 'Exact gauge algebra and actual frozen production consumer only; failedG8 stays failed. No corrected native, physical convergence, continuous resource monitoring or tested restoration claim.'}
    atomic_json(args.output / 'receipt.json', receipt)
    print('FERMI_GAUGE_ORACLE_COMPLETE', len(ledger), digest(args.output / 'receipt.json'), flush=True)


if __name__ == '__main__':
    main()
