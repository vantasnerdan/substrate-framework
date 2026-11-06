"""Actual-data regression for the original Gamma integer-representation join.

Uses the real saved Gamma descriptor, coordinate archive and accepted factor
receipt. The coverage object projects only fields consumed by gamma_metadata;
this is not a new full-grid admission, population fill, state or propagation.
Run only under the bound Root018 five-second input-oracle allowance.
"""
from dataclasses import replace
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT / 'src', HERE.parent / '0025', HERE.parent / '0024'):
    sys.path.insert(0, str(directory))

from material_pilot_io import array_record
from material_reference_output import load_arrays
from material_run_contract import binding, bound_json
from substrate_framework.material_bloch_records import decode_records
from substrate_framework.material_preparation_gamma import gamma_metadata


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    stop = json.loads((HERE / 'root-preparation-join-stop-017.json').read_bytes())
    failed = bound_json(stop['original_failed_terminal'])
    inputs, settings = failed['recipe']['inputs'], failed['recipe']['settings']
    bloch, factor = (bound_json(inputs[key]) for key in ('bloch_terminal', 'factor_terminal'))
    source = next(entry for entry in bloch['outputs'] if entry['key'] == 'source-1-1')
    source_arrays = load_arrays(source)
    records = decode_records(source_arrays['slot_records_json'])
    require(len(records) == 1 and records[0].mesh_index == (0, 0, 0),
            'original Gamma descriptor required')
    record = records[0]
    entry = next(item for item in factor['outputs'] if item['key'] == 'spin1-functional')
    coordinates = entry['receipt']['metadata']['provenance']['original_coordinates']
    legacy = load_arrays(coordinates)
    native = legacy['pw_miller']
    canonical = native.astype(np.int64)
    require(array_record(native) == record.source_miller_record,
            'actual original Native integer32 source identity')
    require(array_record(canonical) == record.target_miller_record,
            'actual stored canonical target matches exact ordered integer embedding')
    require(np.array_equal(native, canonical) and native.dtype == np.dtype('int32'),
            'lossless full original integer coordinates')
    require(record.source_miller_record != record.target_miller_record,
            'the real regression exercises distinct typed byte records')
    reference = SimpleNamespace(cell_bohr=source_arrays['cell_bohr'],
        reciprocal_bohr_inv=source_arrays['reference_reciprocal_bohr_inv'])

    def coverage(selected):
        return SimpleNamespace(slot_records=(selected,), source_receipts={(1, 1): source},
            grid=SimpleNamespace(reference=reference), coverage_identity=bloch['recipe_sha256'])

    gamma_metadata(inputs, settings, coverage(record), factor)
    changed_value = canonical.copy()
    changed_value[0, 0] += 1
    changed_order = canonical[:, ::-1].copy()
    require(not np.array_equal(changed_order, canonical), 'row-order mutation must change rows')
    candidates = {
        'native_storage_mislabelled_as_target': record.source_miller_record,
        'changed_integer_coordinate': array_record(changed_value),
        'permuted_coordinate_order': array_record(changed_order),
    }
    exposed = {}
    for name, candidate in candidates.items():
        try:
            gamma_metadata(inputs, settings, coverage(replace(record,
                target_miller_record=candidate)), factor)
        except ValueError as error:
            require(str(error) == 'original coordinate order/cell/reciprocal identity does not join',
                    name + ': mutation failed at an unrelated gate')
            exposed[name] = str(error)
        else:
            raise AssertionError(name + ': wrong target identity accepted')
    files = (Path(__file__), ROOT / 'src/substrate_framework/material_preparation_gamma.py',
             ROOT / 'src/substrate_framework/material_bloch_records.py')
    source_bindings = [binding(path) for path in files]
    print(json.dumps(dict(scope=__doc__, status='passed',
        authority=binding(HERE / 'root-preparation-join-repair-018.json'),
        original_failed_terminal=stop['original_failed_terminal'],
        original_bloch_terminal=inputs['bloch_terminal'],
        original_factor_terminal=inputs['factor_terminal'],
        original_coordinate_archive=coordinates['archive'], original_source_archive=source['archive'],
        source_binding=record.source_binding, native_record=record.source_miller_record,
        canonical_target_record=record.target_miller_record, exact_integer_order_preserved=True,
        actual_gamma_join_accepted=True, mutation_exposures=exposed,
        source_bindings=source_bindings, numpy=np.__version__), indent=2))


if __name__ == '__main__':
    main()
