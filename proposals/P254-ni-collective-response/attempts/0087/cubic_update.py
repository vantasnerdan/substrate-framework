"""Complete retained all64 cubic - original q0 + actual-provider q0.

No nonzero-q contraction is rerun. The original nine physical 0/15/30 pairs
are copied literally; only the missing 952 original pairs are computed. Three
independent durable stores retain old0/new0/family numerical pairs, so resume
skips each committed stage. The old-q0.json cache may be supplied as old_q0_paid
on later iterations, avoiding every original-q0 contraction thereafter.

API: produce_cubic_update(provider, base_catalog, old_q0_paid, output, *, resume=False)
CLI: python cubic_update.py --provider G1 --base-catalog OLD --old-q0-paid PAID --output NEW [--resume]
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from produce_joint import (HERE, PAIRS, PairStore, atomic_save, binding, authenticate,
                           contract, material_reference, output_lock, validate_pair,
                           retain_dependencies)
import numpy as np
from history_banks import ORDERS
from history_sources import open_reference
from cubic_contour import CubicContour, ordered_matrix


def numeric_identity(value):
    value = np.ascontiguousarray(value)
    return dict(shape=list(value.shape), dtype=str(value.dtype), bytes=value.nbytes,
        sha256=hashlib.sha256(memoryview(value).cast('B')).hexdigest(),
        hash_scope='C-order numeric payload')


def mapped_orders(descriptors, shape):
    arrays = {}
    try:
        if set(descriptors) != set(ORDERS):
            raise ValueError('all three original contour orders required')
        for order, item in descriptors.items():
            arr = np.load(authenticate(item), mmap_mode='r', allow_pickle=False)
            arrays[order] = arr
            if (arr.shape != shape or arr.dtype != np.complex128 or
                    item['shape'] != list(shape) or item['dtype'] != 'complex128'):
                raise ValueError('authenticated cubic dimensions/dtype differ')
        return arrays
    except BaseException:
        close_orders(arrays)
        raise


def close_orders(arrays):
    for arr in arrays.values():
        arr._mmap.close()
    arrays.clear()


def full_pairs(value):
    if value == 'all':
        return True
    return (isinstance(value, list) and len(value) == 961 and
            {tuple(p) for p in value} == set(PAIRS))


def read_base(base_catalog, prereg):
    if binding(base_catalog) != prereg['complete_old_bare_catalog']:
        raise ValueError('authenticated retained original all64 catalog required')
    base = json.loads(Path(base_catalog).read_bytes())
    cov = base['coverage']['cubic']
    descriptors = base['families']['cubic']['arrays']['cubic']
    if (base['status'] != 'complete' or
            base['input_dataset'] != prereg['original_numerical_reference'] or
            base['times_inverse_ry'] != list(range(31)) or
            cov['completed_q'] != list(range(64)) or cov['owned_q'] != list(range(64)) or
            not full_pairs(cov['completed_pairs']) or not full_pairs(cov['owned_pairs']) or
            descriptors != prereg['actual_retained_all64_cubic']):
        raise ValueError('retained original family lineage/coverage differs')
    return descriptors, mapped_orders(descriptors, (31,31,111,111))


def read_old_paid(path, prereg):
    identity = binding(path)
    record = json.loads(Path(path).read_bytes())
    if record['schema'] == 'p254-additive-history-component-v1':
        if (identity != prereg['actual_paid_old_q0_nine'] or
                record['component'] != 'cubic' or record['owned_q'] != [0] or
                record['completed_q'] != [0] or record['status'] != 'complete-owned-component-q0' or
                record['input_dataset'] != prereg['original_numerical_reference'] or
                record['times_inverse_ry'] != [0,15,30] or
                record['selected_original_time_indices'] != [0,15,30] or
                record['external_q'] != 0 or record['internal_q'] != 0 or record['partner_q'] != 0 or
                set(map(tuple,record['coverage']['completed_pairs'])) !=
                {(a,b) for a in range(3) for b in range(3)}):
            raise ValueError('literal authenticated original nine q0 pairs required')
        inputs = json.loads(authenticate(record['input_bindings']).read_bytes())
        coverage = np.load(authenticate(record['coverage']['array']), mmap_mode='r', allow_pickle=False)
        try:
            if coverage.shape != (3,3,1) or coverage.dtype != np.bool_ or not coverage.all():
                raise ValueError('paid original nine bitmap incomplete')
        finally:
            coverage._mmap.close()
        indices = (0,15,30)
        pairs = {(a,b):(i,j) for i,a in enumerate(indices) for j,b in enumerate(indices)}
        arrays = mapped_orders(record['arrays'], (3,3,111,111))
    elif record['schema'] == 'p254-cubic-q0-cache-v1':
        cov = record['coverage']
        if (record['status'] != 'complete' or
                record['input_dataset'] != prereg['original_numerical_reference'] or
                record['preparation_origin'] != prereg['original_preparation'] or
                record['numerical_reference'] != prereg['original_numerical_reference'] or
                record['original_paid_q0'] != prereg['actual_paid_old_q0_nine'] or
                record['owned_q'] != [0] or record['completed_q'] != [0] or
                record['times_inverse_ry'] != list(range(31)) or
                cov['owned_q'] != [0] or cov['completed_q'] != [0] or
                cov['all111_coordinates'] is not True or
                not full_pairs(cov['owned_pairs']) or not full_pairs(cov['completed_pairs']) or
                record['route'] != dict(external_q=0,internal_q=0,partner_q=0)):
            raise ValueError('complete restored original-q0 cache required')
        coverage = np.load(authenticate(cov['array']), mmap_mode='r', allow_pickle=False)
        try:
            if coverage.shape != (31,31) or coverage.dtype != np.bool_ or not coverage.all():
                raise ValueError('restored original q0 bitmap incomplete')
        finally:
            coverage._mmap.close()
        inputs = record['cubic_inputs']
        pairs = {p:p for p in PAIRS}
        arrays = mapped_orders(record['arrays'], (31,31,111,111))
    else:
        raise ValueError('original paid q0 component or complete restored cache required')
    try:
        for physical, local in pairs.items():
            validate_pair({o:arr[local] for o,arr in arrays.items()}, *physical)
    except BaseException:
        close_orders(arrays)
        raise
    return identity, pairs, arrays, inputs


def q0_values(reference, a, b, partner, left, right):
    # Raw ordered legs are requested independently, exactly as A76. Neither
    # the reverse vertex nor any D leg is manufactured by a dagger/refill.
    first = reference.boson(0,a,b)
    second = reference.boson(partner,a,b)
    greater = ordered_matrix(left,right,first[0],second[0])
    lesser = ordered_matrix(left,right,first[1],second[1])
    retarded = greater-lesser if reference.times[a] >= reference.times[b] else np.zeros_like(greater)
    return dict(greater=greater,lesser=lesser,retarded=retarded)


def cache_record(store, source, prereg, cubic_inputs, original_paid, paid_source):
    completed = store.completed()
    complete = len(completed) == 961
    return dict(schema='p254-cubic-q0-cache-v1',status='complete' if complete else 'partial',
        input_dataset=source, numerical_reference=prereg['original_numerical_reference'],
        preparation_origin=prereg['original_preparation'], times_inverse_ry=list(range(31)),
        arrays=store.descriptors('cubic'), owned_q=[0], completed_q=[0] if complete else [],
        original_paid_q0=original_paid, reused_source=paid_source,
        cubic_inputs=cubic_inputs, route=dict(external_q=0,internal_q=0,partner_q=0),
        coverage=dict(owned_q=[0],completed_q=[0] if complete else [],
            owned_pairs=[list(p) for p in PAIRS],completed_pairs=completed,
            all111_coordinates=True,array=store.bitmap()),
        full_time_square=complete, full_solve_eligible=False,
        scope='One genuine internal q0 contribution, never the all64 cubic family')


def validate_family(catalog):
    """Read/authenticate actual saved outputs and exact pairwise replacement."""
    record = json.loads(Path(catalog).read_bytes())
    if record['status'] != 'complete':
        raise ValueError('complete persisted cubic family required')
    manifests = {name:json.loads(authenticate(item).read_bytes())
                 for name,item in record['provenance'].items()}
    banks = {}
    try:
        banks['family'] = mapped_orders(record['arrays'],(31,31,111,111))
        for name in ('base','old_q0','new_q0'):
            banks[name] = mapped_orders(manifests[name]['arrays'],(31,31,111,111))
        for a,b in PAIRS:
            for bank in banks.values():
                validate_pair({o:arr[a,b] for o,arr in bank.items()},a,b)
            for o in ORDERS:
                expected = banks['base'][o][a,b]-banks['old_q0'][o][a,b]+banks['new_q0'][o][a,b]
                if not np.array_equal(banks['family'][o][a,b],expected):
                    raise ValueError('persisted all64-old0+new0 identity differs')
        return dict(schema='p254-a87-persisted-cubic-validation-v1',status='passed',
            catalog=binding(catalog),pairs_read=961,all111_coordinates=True,
            inherited_q=list(range(1,64)),replaced_q=[0],
            exact_saved_replacement=True, nonzero_q_contractions=0)
    finally:
        for bank in banks.values():
            close_orders(bank)


def produce_cubic_update(provider, base_catalog, old_q0_paid, output, *, resume=False):
    prereg = contract()
    source = binding(provider)
    request = dict(input_dataset=source, base_catalog=binding(base_catalog),
        old_q0_paid=binding(old_q0_paid), preparation_origin=prereg['original_preparation'],
        numerical_reference=prereg['original_numerical_reference'],
        implementation=[binding(HERE/n) for n in ('cubic_update.py','produce_joint.py')])
    for key in ('original_preparation','original_numerical_reference','actual_paid_old_q0_nine'):
        authenticate(prereg[key])
    began = time.monotonic()
    with output_lock(output,resume) as out:
        if resume:
            record = json.loads((out/'checkpoint.json').read_bytes())
            if record.get('request') != request or record.get('schema') != 'p254-dressed-cubic-family-v1':
                raise ValueError('resume source/request/implementation differs')
            for item in record.get('implementation_bindings', ()):
                authenticate(item)
        else:
            record = dict(schema='p254-dressed-cubic-family-v1',status='partial',request=request,
                input_dataset=source, numerical_reference=prereg['original_numerical_reference'],
                preparation_origin=prereg['original_preparation'], times_inverse_ry=list(range(31)),
                axes=['time','time','ionic_coordinate','ionic_coordinate'],
                inherited_q=list(range(1,64)),replaced_q=[0], full_solve_eligible=False,
                full_time_square=False, full_mesh=False,
                combination='Cbare_all64 - original_q0 + actual_provider_q0',
                coverage=dict(owned_q=list(range(64)),completed_q=[],
                    owned_pairs=[list(p) for p in PAIRS],completed_pairs=[],
                    all111_coordinates=True,
                    checkpoint='old0/new0/family each numeric flush/fsync then bitmap flush/fsync; atomic JSON+directory fsync'),
                stage_completed_pairs=dict(old_q0=[],new_q0=[],family=[]),timings=[],
                boundary_loss=dict(incomplete_pair=None,incomplete_stage=None,
                    static_rebuild='original material, bare transitions and q0 vertices rebuilt; bitmap-committed contractions skipped'),
                scope='Complete all64 ionic cubic family by q0 replacement only, not a single-q catalog or complete collision map')
        stores, base_arrays, paid_arrays = {}, {}, {}
        ref = bare = None
        cubic_inputs = None
        active = stage = None
        try:
            for name in ('old_q0','new_q0','family'):
                stores[name] = PairStore(out/name, {'cubic-'+o:(111,111) for o in ORDERS},resume=resume)
            if np.any(stores['family'].seen & ~(stores['old_q0'].seen & stores['new_q0'].seen)):
                raise ValueError('persisted family coverage precedes its old0/new0 dependencies')
            record['stage_completed_pairs'] = {k:s.completed() for k,s in stores.items()}
            record['coverage']['completed_pairs'] = stores['family'].completed()
            atomic_save(out/'checkpoint.json',record)
            descriptors, base_arrays = read_base(base_catalog,prereg)
            paid_identity, paid_pairs, paid_arrays, paid_inputs = read_old_paid(old_q0_paid,prereg)
            c, ref, coefficients = material_reference(provider)
            bare = open_reference(c,coefficients['manifest'],bundle=authenticate(prereg['original_numerical_reference']),times=np.arange(31))
            cubic_source = authenticate(paid_inputs['cubic_source'])
            cubic = CubicContour.from_source(cubic_source)
            partner,left,right = cubic.pair(0,0)
            routed = [int(cubic.source.qindex(cubic.source.indices[0]-cubic.source.indices[k])) for k in range(64)]
            if (partner != 0 or left.shape != (111,111,111) or right.shape != left.shape or
                    [k for k,l in enumerate(routed) if k == 0 or l == 0] != [0]):
                raise ValueError('q0 replacement is not a disjoint original all64 routing')
            if numeric_identity(left) != paid_inputs['left'] or numeric_identity(right) != paid_inputs['right']:
                raise ValueError('independently rebuilt raw q0 vertices differ from original paid input')
            cubic_inputs = dict(cubic_source=paid_inputs['cubic_source'],
                cubic_source_arrays=cubic.source.source_record['arrays'],
                cubic_source_bindings=cubic.source.source_record['source_bindings'],
                left=numeric_identity(left),right=numeric_identity(right),
                coefficient_manifest=coefficients['manifest'])
            dependencies = [binding(HERE/'cubic_update.py'),
                binding(HERE/'produce_joint.py'), binding(Path(ordered_matrix.__code__.co_filename)),
                binding(Path(open_reference.__code__.co_filename)),
                binding(Path(ref.electron.__func__.__code__.co_filename))]
            retain_dependencies(record, dependencies, stores.values())
            record.update(cubic_inputs=cubic_inputs,source_bindings=ref.bindings,
                original_source_bindings=bare.bindings,
                routed_partner_q=routed,route=dict(external_q=0,internal_q=0,partner_q=0),
                paid_original_reuse=dict(source=paid_identity,
                    physical_pairs=[list(p) for p in sorted(paid_pairs)],literal_copy=True),
                original_q0_missing_count=961-len(paid_pairs),nonzero_q_contractions=0)
            atomic_save(out/'checkpoint.json',record)
            atomic_save(out/'retained-base.json',dict(schema='p254-retained-cubic-base-v1',status='complete',
                arrays=descriptors,base_catalog=request['base_catalog'],
                original_input_dataset=prereg['original_numerical_reference'],
                preparation_origin=prereg['original_preparation'],times_inverse_ry=list(range(31)),
                coverage=dict(owned_q=list(range(64)),completed_q=list(range(64)),
                              owned_pairs='all',completed_pairs='all',all111_coordinates=True)))
            for a,b in PAIRS:
                if all(store.seen[a,b] for store in stores.values()):
                    continue
                active = [a,b]
                tick = time.monotonic()
                for stage in ('old_q0','new_q0','family'):
                    store = stores[stage]
                    if store.seen[a,b]:
                        continue
                    if stage == 'old_q0':
                        if (a,b) in paid_pairs:
                            local = paid_pairs[a,b]
                            values = {o:paid_arrays[o][local] for o in ORDERS}
                        else:
                            values = q0_values(bare,a,b,partner,left,right)
                    elif stage == 'new_q0':
                        values = q0_values(ref,a,b,partner,left,right)
                    else:
                        if not stores['old_q0'].seen[a,b] or not stores['new_q0'].seen[a,b]:
                            raise ValueError('family cannot precede durable old0/new0 pair')
                        values = {o:base_arrays[o][a,b]-stores['old_q0'].arrays['cubic-'+o][a,b]
                                    +stores['new_q0'].arrays['cubic-'+o][a,b] for o in ORDERS}
                    store.commit(a,b,{'cubic-'+o:values[o] for o in ORDERS})
                    record['stage_completed_pairs'] = {k:s.completed() for k,s in stores.items()}
                    record['coverage']['completed_pairs'] = stores['family'].completed()
                    record['elapsed_seconds_this_invocation'] = time.monotonic()-began
                    atomic_save(out/'checkpoint.json',record)
                record['timings'].append(dict(pair=[a,b],seconds=time.monotonic()-tick))
                active = stage = None
                print(json.dumps(dict(event='actual-dressed-cubic-family-pair-committed',pair=[a,b],
                    completed_pairs=len(stores['family'].completed()))),flush=True)
            record['failure'] = None
        except BaseException as exc:
            record['failure'] = dict(type=type(exc).__name__,message=str(exc))
            record['boundary_loss'].update(incomplete_pair=active,incomplete_stage=stage)
            raise
        finally:
            if bare is not None:
                bare.close()
            if ref is not None:
                ref.close()
            close_orders(base_arrays)
            close_orders(paid_arrays)
            try:
                if len(stores) == 3:
                    completed = stores['family'].completed()
                    complete = (len(completed) == 961 and all(len(s.completed()) == 961 for s in stores.values()))
                    record.update(status='complete' if complete else 'partial',
                        full_time_square=complete,full_mesh=complete,
                        elapsed_seconds_this_invocation=time.monotonic()-began,
                        arrays=stores['family'].descriptors('cubic'))
                    record['stage_completed_pairs'] = {k:s.completed() for k,s in stores.items()}
                    record['coverage'].update(completed_pairs=completed,
                        completed_q=list(range(64)) if complete else [],
                        array=stores['family'].bitmap())
                    if cubic_inputs is not None:
                        atomic_save(out/'old-q0.json',cache_record(stores['old_q0'],
                            prereg['original_numerical_reference'],prereg,cubic_inputs,
                            prereg['actual_paid_old_q0_nine'],request['old_q0_paid']))
                        atomic_save(out/'new-q0.json',cache_record(stores['new_q0'],source,
                            prereg,cubic_inputs,prereg['actual_paid_old_q0_nine'],None))
                        record['provenance'] = dict(base=binding(out/'retained-base.json'),
                            old_q0=binding(out/'old-q0.json'),new_q0=binding(out/'new-q0.json'))
                    atomic_save(out/'checkpoint.json',record)
                    atomic_save(out/'catalog.json',record)
            finally:
                for store in stores.values():
                    store.close()
        atomic_save(out/'validation.json',validate_family(out/'catalog.json'))
        return out/'catalog.json'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--provider',type=Path,required=True)
    parser.add_argument('--base-catalog',type=Path,required=True)
    parser.add_argument('--old-q0-paid',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args()
    print(produce_cubic_update(args.provider,args.base_catalog,args.old_q0_paid,args.output,resume=args.resume))
