#!/usr/bin/env python3
"""Exact full-shape zero values and unthresholded real complex128 payloads."""
import argparse
import json
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
from mixed_payload_output import write_payload
from mixed_response_storage import COMPONENTS, MappedMixedState, publish_checkpoint
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.streamed_material_state import authenticate, file_binding


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    work = args.output.with_name(args.output.stem+'-work')
    work.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1, max_seconds=30., max_bytes=64<<20,
                     max_memory_bytes=512<<20, worker_memory_bytes=256<<20)
    start = perf_counter()
    checks, outputs = [], []

    def require(name, condition, **details):
        checks.append(dict(name=name, passed=bool(condition), **details))
        if not condition:
            raise AssertionError(name)

    shape = (768, 1024)
    zero = np.full(shape, complex(-0.,0.), dtype=np.complex128)
    binding = write_payload(work/'computed-zero.npy', zero, budget)
    outputs.append(binding)
    small = {c: {name: np.zeros(3 if name in ('u', 'p') else 1, complex)
                 for name in ('u', 'p', 'coordinate', 'momentum')} for c in COMPONENTS}
    checkpoint = publish_checkpoint(work/'zero-consumer.json',
        legs={c: [binding] for c in COMPONENTS}, small=small, q=(0., 0., 0.),
        target_plus=[0], target_minus=[0], time=0., identity='computed-zero-consumer',
        kind='computed-zero', provenance=dict(operation='payload-output-oracle'), budget=budget)
    outputs.append(checkpoint)
    with MappedMixedState(checkpoint, budget=budget) as state:
        value = state.leg('plus', 0)
        require('entire_computed_zero_matrix_retains_shape_dtype_values',
                value.shape==shape and value.dtype==np.complex128
                and np.array_equal(value,zero) and not value.flags.writeable)
    require('exact_zero_has_no_fictitious_payload_file', not (work/'computed-zero.npy').exists())
    # A tiny nonzero and complex signs must survive, not take a zero shortcut.
    data = np.zeros((97, 61), dtype=np.complex128, order='F')
    data[31,27] = 1.e-300-2.e-300j
    data[90,3] = -.125+.375j
    actual = write_payload(work/'nonzero.npy', data, budget)
    outputs.append(actual)
    value = np.load(authenticate(actual), mmap_mode='r', allow_pickle=False)
    require('tiny_nonzero_and_all_complex_values_are_exact', np.array_equal(value,data)
            and value[31,27] != 0 and value.shape==data.shape and value.dtype==data.dtype)
    value._mmap.close()
    for label, bad in [('nan',np.full((2,3),np.nan+0j)),
                       ('infinite',np.full((2,3),complex(0,np.inf)))]:
        try:
            write_payload(work/(label+'.npy'), bad, budget)
        except ValueError:
            require('reject_'+label+'_before_publication', not (work/(label+'.npy')).exists())
        else:
            raise AssertionError('accepted '+label)
    before = file_binding(work/'nonzero.npy')
    try:
        write_payload(work/'nonzero.npy', np.zeros(data.shape, complex), budget)
    except FileExistsError:
        require('existing_original_never_overwritten', file_binding(work/'nonzero.npy')==before)
    else:
        raise AssertionError('overwrote an existing payload')
    result = dict(schema='P254-exact-zero-computed-payload-oracle-v2',status='PASS',
        scope='actual numeric storage boundaries only; no material propagation',checks=checks,
        outputs=outputs,budget=budget.snapshot(),wall_seconds=perf_counter()-start,
        process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        source_bindings=[file_binding(HERE/'mixed_payload_output.py'),
                         file_binding(HERE/'mixed_response_storage.py'),file_binding(__file__)])
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(status='PASS',checks=len(checks),output=str(args.output))))


if __name__ == '__main__':
    main()
