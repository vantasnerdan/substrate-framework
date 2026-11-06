#!/usr/bin/env python3
"""Build current case metadata, reusing only byte-equal original frame/D arrays.

Different basis labels never license reuse: every candidate is compared with
freshly reconstructed current-source values, in current full row/column order.
Different arrays are saved, not relabeled. Original states remain unchanged.
"""
import argparse
import json
import os
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import matched_reference_io as io
from matched_joint_case import load_case, model_identity
from matched_joint_state import original_sectors
from substrate_framework.streamed_joint_state import JointBudget
from substrate_framework.streamed_material_state import DiskMean, authenticate, file_binding


def make(args, budget):
    binding = file_binding(args.original)
    if binding['sha256'] != args.original_sha256:
        raise ValueError('original reuse source changed')
    context = load_case(args.case)
    old = DiskMean(args.original, verify=False)
    by_id = {slot.identity: slot for slot in old.slots}
    if set(by_id) != set(context.expected_slots):
        raise ValueError('complete identical slot names required for reuse candidates')
    directory = args.output/'original'
    directory.mkdir()
    rows, comparisons = [], []
    reused_bytes = 0
    for index, (metadata, psi, d) in enumerate(original_sectors(context, budget)):
        slot = by_id[metadata['identity']]
        old_values = None
        if slot.shape == psi.shape:
            authenticate(slot.frame)
            authenticate(slot.coefficient)
            old_values = old.arrays(slot)
        record = dict(metadata, shape=list(psi.shape))
        reused = {}
        for name, value, candidate_index, candidate in (
                ('frame', psi, 0, slot.frame), ('coefficient', d, 1, slot.coefficient)):
            same = old_values is not None and np.array_equal(value, old_values[candidate_index])
            if same:
                record[name] = candidate
                reused_bytes += int(candidate['bytes'])
            else:
                record[name] = budget.array(directory/('%s-%04d.npy' % (name, index)), value)
            reused[name] = bool(same)
        rows.append(record)
        comparisons.append(dict(slot=metadata['identity'], reused=reused,
            old_metadata=slot.metadata(), current_metadata=dict(metadata, shape=list(psi.shape))))
    if tuple(row['identity'] for row in rows) != context.expected_slots:
        raise ValueError('current complete original slot coverage differs')
    evidence = dict(case=args.case, context_identity=context.identity, source=file_binding(__file__),
        candidate_boundary=binding, source_exact_array_comparisons=comparisons,
        state_transition='Current original source metadata and values; prior labels not adopted',
        no_gram_inverse=True, no_refilling=True, full_raw_columns=True)
    budget.reserve(len(json.dumps(evidence).encode())+128*4096)
    zero = np.zeros_like(context.reference.masses)
    identity = model_identity(dict(kind='source-exact-rebound-original', context=context.identity,
                                  candidate=binding, source=file_binding(__file__)))
    DiskMean._publish(directory, rows, identity, 0., zero, zero, evidence)
    new = DiskMean(directory/'boundary.json', expected_identity=identity)
    return dict(status='complete', case=args.case, context_identity=context.identity,
        original_boundary=file_binding(new.path), candidate_boundary=binding,
        source=file_binding(__file__), comparisons=comparisons,
        reused_logical_array_bytes=str(reused_bytes), full_slots=len(rows),
        scope='Exact current-source raw electronic state; no coefficient, H0 or trajectory computation')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', required=True)
    p.add_argument('--original', type=Path, required=True)
    p.add_argument('--original-sha256', required=True)
    p.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    io.runtime()
    args.output.mkdir(parents=True, exist_ok=False)
    budget = JointBudget(max_actions=1, max_seconds=600, max_bytes=20*1024**3)
    budget.reserve(1024**2)
    result = dict(status='partial', source=file_binding(__file__), case=args.case)
    try:
        result = make(args, budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    (args.output/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
