"""Paid q1 first column, genuine0..30; all128/64/12, no endpoint replay."""
from pathlib import Path
import argparse
import json
import resource
import sys
from time import perf_counter
from column_common import (HERE, BUNDLE, CHANNELS, begin, authenticate,
                           atomic_json, file_binding)
import numpy as np
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors
from general_moving_reference import SavedMovingReference
from streamed_ee import stream_ee, NumericSink, read_blocks


def produce(output):
    started = perf_counter()
    out, spec = begin(output, 6)
    paid_path = HERE.parent/'0064/pilot-001/result.json'
    paid = json.loads(paid_path.read_bytes())
    endpoint = authenticate(spec['saved_endpoints_ee'])
    if spec['saved_endpoints_ee'] not in paid['outputs'] or paid['status'] != 'complete':
        raise ValueError('complete original paid endpoint output required')
    c = load_coulomb_inputs()
    factors = BareCoulombFactors(c)
    reference = SavedMovingReference(c, BUNDLE, c.bindings['manifest'])
    if (c.bindings != paid['inputs']['material'] or
            file_binding(BUNDLE) not in paid['inputs']['reference'] or
            reference.spec['U'] not in paid['inputs']['reference'] or
            not np.array_equal(reference.times, np.arange(31, dtype=float))):
        raise ValueError('same actual material and genuine U31 source required')
    atomic_json(out/'actual-inputs.json', dict(material=c.bindings,
        reference=reference.bindings, paid_endpoint=file_binding(endpoint)))
    arrays = {name:np.empty((31,128,12,12), dtype=complex) for name in CHANNELS}
    seen = np.zeros((31,128), dtype=bool)
    with np.load(endpoint, allow_pickle=False) as archive:
        for name in CHANNELS:
            old = archive[name+'_total']
            if old.shape != (128,2,2,12,12):
                raise ValueError('paid all-slot ordered endpoint axes differ')
            arrays[name][0] = old[:,0,0]
            arrays[name][30] = old[:,1,0]
    seen[[0,30]] = True
    sink = NumericSink(out/'new-blocks.npyblocks')
    blocks = 0
    def consume(block):
        nonlocal blocks
        sink(block)
        for column,(a,b) in enumerate(block.pairs):
            if b != 0 or np.any(seen[a,list(block.external)]):
                raise ValueError('duplicated or unowned physical time/slot row')
            for name in CHANNELS:
                arrays[name][a,list(block.external)] = block.values[name+'_total'][:,column]
            seen[a,list(block.external)] = True
        blocks += 1
        if blocks % 32 == 0:
            receipt = dict(elapsed_seconds=perf_counter()-started,
                completed_time_slot_rows=int(seen.sum()),required_time_slot_rows=3968,
                actual_blocks=blocks,physical_scope='q1 first column only')
            atomic_json(out/'progress.json', receipt)
            with (out/'telemetry.jsonl').open('a') as stream:
                stream.write(json.dumps(receipt)+'\n')
            print(json.dumps(receipt), flush=True)
    tick = perf_counter()
    try:
        manifest = stream_ee(factors, reference.electron, reference.times, consume,
            owned_q=(1,), owned_pairs=tuple((a,0) for a in range(1,30)),
            pair_batch=2, external_batch=2, static_slots=16,
            source_bindings=reference.bindings)
        stream_record = sink.close(manifest)
        seconds = perf_counter()-tick
        if not seen.all() or not all(np.isfinite(value).all() for value in arrays.values()):
            raise ValueError('complete finite owned first column not produced')
        archive_path = out/'column.npz'
        np.savez(archive_path, times=reference.times,
            ordered_pairs=np.array([(a,0) for a in range(31)]), outer_q=np.array(1),
            **{name+'_total':value for name,value in arrays.items()})
        persisted_blocks = sum(1 for _ in read_blocks(out/'new-blocks.npyblocks'))
        if persisted_blocks != manifest['blocks']:
            raise ValueError('paid numeric stream incomplete')
        with np.load(archive_path, allow_pickle=False) as saved:
            consumer = dict(retarded_initial_column_action_norm=float(np.linalg.norm(
                saved['retarded_total']@c.rho0)))
        record = dict(schema='p254-paid-ee-long-history-column-v1',
            status='complete-first-column-q1',source_reference=file_binding(BUNDLE),
            inputs=file_binding(out/'actual-inputs.json'),paid_endpoint=file_binding(endpoint),
            output=file_binding(archive_path),stream=stream_record,
            times_inverse_ry=reference.times.tolist(),outer_q=1,external_slots=128,
            inner_h=64,bands=12,spins=2,completed_time_slot_rows=int(seen.sum()),
            reused_time_pairs=[[0,0],[30,0]],new_time_pairs=[[a,0] for a in range(1,30)],
            new_evaluation_seconds=seconds,elapsed_seconds=perf_counter()-started,
            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            consumer=consumer,full_time_square=False,full_outer_mesh=False,
            source_bindings=[file_binding(module.__file__) for module in list(sys.modules.values())
                if getattr(module,'__file__',None) and str(HERE.parent) in str(module.__file__)
                and str(module.__file__).endswith('.py')],
            limits=['First column/q1 of Gaussian first map, not complete intrinsic evolution.',
                    'No new performance method, approximation, state repair or conservation claim.'])
        atomic_json(out/'result.json', record)
        print('A66_EE_NATIVE_COMPLETE '+json.dumps(file_binding(out/'result.json')), flush=True)
    finally:
        reference.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    produce(parser.parse_args().output)
