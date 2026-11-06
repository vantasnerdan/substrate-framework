"""Reuse a complete paid A56 map in the canonical A62 NPY representation."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0062'))
from history_banks import FinalBanks, OwnedKernels, ORDERS, authenticate, binding, REFERENCE, PI_CONVENTION


def transcode(catalog, output):
    catalog = Path(catalog).resolve()
    old = json.loads(catalog.read_bytes())
    if (old['schema'] != 'p254-dressed-functional-components-v1' or old['status'] != 'complete'
            or old['reference_source_kind'] != REFERENCE or old['pi_convention'] != PI_CONVENTION):
        raise ValueError('complete original common functional map required')
    times = np.asarray(old['times_inverse_ry'], dtype=float)
    if times.shape != (3,):
        raise ValueError('this historical producer has exactly three genuine nodes')
    ee_path = authenticate(old['complete_ee_direct_exchange_total'])
    mixed_binding = old['component_bindings']['mixed']
    cubic_binding = old['component_bindings']['cubic']
    mixed = json.loads(authenticate(mixed_binding).read_bytes())
    cubic = json.loads(authenticate(cubic_binding).read_bytes())
    if mixed['status'] != 'complete' or cubic['status'] != 'complete':
        raise ValueError('partial owning components')
    ei_q = [item['q'] for item in mixed['timings']]
    with np.load(ee_path, allow_pickle=False) as ee, \
            np.load(authenticate(cubic['archive']), allow_pickle=False) as pi:
        ee_q = ee['completed_outer_q'].tolist()
        cubic_q = pi['completed_internal_q'].tolist()
        if ee_q != list(range(64)) or ei_q != list(range(64)) or cubic_q != list(range(64)):
            raise ValueError('complete disjoint all64 ownership required')
        if set(map(tuple, pi['updated_pairs'].tolist())) != {(a, b) for a in range(3) for b in range(3)}:
            raise ValueError('all genuine cubic time pairs required')
    for item in (old['actual_saved_numeric_consumer'], old['input_dataset'],
                 mixed['archive'], mixed['joint_archive'], mixed['ei_archive']):
        authenticate(item)
    archive = authenticate(old['archive'])
    coverage = {name: dict(owned_pairs='all', completed_pairs='all')
                for name in ('ee', 'mixed_ei', 'joint', 'cubic')}
    for name, owned in (('ee', ee_q), ('mixed_ei', ei_q), ('cubic', cubic_q)):
        coverage[name].update(owned_q=owned, completed_q=owned)
    with FinalBanks(output, times) as banks, np.load(archive, allow_pickle=False) as paid:
        for kind in ('sigma', 'pi'):
            for order in ORDERS:
                array = paid[kind+'_'+order]
                target = banks.banks[kind][order]
                if array.shape != target.shape or array.dtype != np.complex128 or not np.isfinite(array).all():
                    raise ValueError('complete paid axes/precision/values required')
                np.copyto(target, array)
        banks.publish(dict(times_inverse_ry=times.tolist(), input_dataset=old['input_dataset'],
            coverage=coverage, source_bindings=[binding(catalog), old['archive'],
                old['complete_ee_direct_exchange_total'], mixed_binding, cubic_binding, binding(__file__)],
            approximation=old['exact_physics_scope'], paid_functional_map=binding(catalog),
            computation='Representation-only reuse of all six paid kernels; zero new functional evaluations.'))
    result = Path(output)/'catalog.json'
    with OwnedKernels(result) as saved, np.load(archive, allow_pickle=False) as paid:
        for kind in ('sigma', 'pi'):
            for order in ORDERS:
                if not np.array_equal(saved.banks[kind][order], paid[kind+'_'+order]):
                    raise ValueError('saved kernel representation changed numeric values')
    return binding(result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(transcode(args.catalog, args.output)), flush=True)
