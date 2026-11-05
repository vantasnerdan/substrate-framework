"""Pack disjoint paid EE sectors; no physics contraction or normalization."""
from pathlib import Path
import json
import os
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0063'))
from history_storage import StoredNPZ, authenticate, binding, descriptor
from production_transport import atomic_json

ORDERS = ('greater', 'lesser', 'retarded')
RESULTS = (
    Path('/var/tmp/p254-a67-ee-full31-001/result.json'),
    Path('/var/tmp/p254-a70-ee-full31-fused-001/result.json'),
    Path('/var/tmp/p254-a73-ee-q1-complete-001/result.json'),
)
BUNDLE = Path('/var/tmp/p254-a52-requested-reference-001/producer/bundle.json')


def pack(output):
    output = Path(output).resolve()
    output.mkdir(exist_ok=False)
    reference = binding(BUNDLE)
    records, seen, inputs = [], set(), []
    for path in RESULTS:
        record = json.loads(path.read_bytes())
        slots = record['slots']
        actual_inputs = json.loads(authenticate(record['inputs']).read_bytes())
        if (record['status'] != 'complete-at-owned-slots-q1' or
                record['genuine_times_inverse_ry'] != list(range(31)) or
                record['outer_q'] != 1 or record['inner_h'] != 64 or
                record['internal_modes'] != 1728 or record['bands'] != 12 or
                record['actual_owned_cells'] != len(slots)*961 or
                reference not in actual_inputs['reference'] or
                len(set(slots)) != len(slots) or seen.intersection(slots) or
                len(record['overlap_findings']) != 3*len(slots) or
                not all(row['passed'] for row in record['overlap_findings'])):
            raise ValueError('paid full-history sector ownership/source differs')
        authenticate(record['output'])
        seen.update(slots)
        records.append(record)
        inputs.append(dict(result=binding(path), numeric_archive=record['output'],
                           actual_inputs=record['inputs'], slots=slots,
                           historical_implementation=record['source_bindings']))
    if seen != set(range(128)):
        raise ValueError('all128 disjoint paid sectors are required')
    comparisons = 0
    descriptions = {}
    for order in ORDERS:
        path = output/('sigma-'+order+'.npy')
        target = np.lib.format.open_memmap(path, mode='w+', dtype='<c16',
                                           shape=(31,31,128,12,12))
        try:
            for record in records:
                with StoredNPZ(record['output']['path']) as source:
                    slots = record['slots']
                    values = source[order]
                    if (values.shape != (len(slots),31,31,12,12) or
                            values.dtype != np.dtype('<c16') or
                            not np.array_equal(source['slots'], slots) or
                            not np.array_equal(source['times'], np.arange(31)) or
                            int(source['outer_q']) != 1):
                        raise ValueError('paid member layout differs')
                    for index, slot in enumerate(slots):
                        if not np.isfinite(values[index]).all():
                            raise ValueError('nonfinite paid sector')
                        target[:,:,slot] = values[index]
            target.flush()
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
        finally:
            target._mmap.close()
        # Reopen the persisted consumer path; compare every actual paid matrix.
        target = np.load(path, mmap_mode='r', allow_pickle=False)
        try:
            for record in records:
                with StoredNPZ(record['output']['path']) as source:
                    for index, slot in enumerate(record['slots']):
                        if not np.array_equal(target[:,:,slot], source[order][index]):
                            raise ValueError('saved axis packing changed paid values')
                        comparisons += 961
        finally:
            target._mmap.close()
        descriptions[order] = descriptor(path)
        descriptions[order]['dtype'] = 'complex128'
    pairs = [[a,b] for a in range(31) for b in range(31)]
    catalog = dict(schema='p254-additive-history-component-v1', component='ee',
        status='complete-owned-component-q1', input_dataset=reference,
        times_inverse_ry=list(map(float,range(31))), owned_q=[1], completed_q=[1],
        slots=128, bands=12, inner_h=64, internal_modes=1728,
        axes=['time','time','slot','band','band'], sigma=descriptions,
        coverage=dict(owned_pairs=pairs, completed_pairs=pairs,
                      external_targets=list(range(128))),
        normalization=dict(original_N=64,
            factor='Literal already-normalized A64 total=direct+exchange; no further factor during packing'),
        input_bindings=[dict(result=x['result'], numeric_archive=x['numeric_archive'],
                             actual_inputs=x['actual_inputs']) for x in inputs],
        paid_provenance=[dict(archive=x['numeric_archive'],metadata=x) for x in inputs],
        implementation_bindings=[binding(__file__)],
        consumer=dict(all_saved_channel_slot_pairs_byte_identical=True,
                      compared_channel_slot_pairs=comparisons),
        full_time_square=True, full_mesh=False, full_solve_eligible=False,
        limits=['One original outerq1 only; not all64 momentum mesh or trajectory.',
                'Packing preserves literal paid values and adds no physics or accuracy claim.'])
    atomic_json(output/'catalog.json', catalog)
    return output/'catalog.json'
