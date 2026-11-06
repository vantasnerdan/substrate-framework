#!/usr/bin/env python3
"""Locate a reported original-state reuse mismatch without changing either state."""
import argparse
import json
import os
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import matched_reference_io as io
from matched_joint_case import load_case
from matched_joint_state import original_sectors
from substrate_framework.streamed_joint_state import JointBudget
from substrate_framework.streamed_material_state import DiskMean, authenticate, file_binding


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
    budget = JointBudget(max_actions=1, max_seconds=180, max_bytes=1024**2)
    original_binding = file_binding(args.original)
    if original_binding['sha256'] != args.original_sha256:
        raise ValueError('original boundary changed')
    context = load_case(args.case)
    original = DiskMean(args.original, verify=False)
    rows = []
    for (metadata, psi, d), slot in zip(original_sectors(context, budget), original.slots):
        authenticate(slot.frame)
        authenticate(slot.coefficient)
        old_psi, old_d = original.arrays(slot)
        expected = dict(metadata, shape=psi.shape)
        changed = {key: dict(saved=slot.metadata()[key], expected=value)
                   for key, value in expected.items() if slot.metadata()[key] != value}
        row = dict(slot=slot.identity, metadata_changes=changed,
                   frame_equal=np.array_equal(psi, old_psi), D_equal=np.array_equal(d, old_d),
                   frame_max_difference=float(np.max(np.abs(psi-old_psi))),
                   D_max_difference=float(np.max(np.abs(d-old_d))))
        rows.append(row)
        if changed or not row['frame_equal'] or not row['D_equal']:
            break
    args.output.mkdir(parents=True, exist_ok=False)
    budget.json(args.output/'result.json', dict(case=args.case, original=original_binding,
        source=file_binding(__file__), context_identity=context.identity, first_difference=rows[-1],
        compared_slots=len(rows), budget=budget.snapshot(), scope='diagnostic only; neither state changed'))


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
