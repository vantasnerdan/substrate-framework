"""Full-Gamma prepared reuse and deterministic, fail-closed full-q manifest union.

This only reuses the approved maximum12-band/27-family physical projection:
all128 explicit spin/k slots, actual prepared populations/raw-D bindings, all64 q,
full original paired ionic support, independent complex reverse coefficients,
and BOTH frozen passive splits. No source product is rebound to a new producer.
Numerical arrays remain in their authenticated absolute immutable archives.

The producer-history directory must contain exact snapshots named
``<original stem>-<sha256><original suffix>`` for every implementation AND
physical_implementation source in the settled Gamma specification. The original
Gamma manifest, its pilot manifest, snapshots and numerical archives must remain
available and immutable. Missing/foreign/colliding inputs are errors, not zeros.
"""
from copy import deepcopy
from contextlib import ExitStack
import fcntl
import json
from pathlib import Path

import numpy as np
from substrate_framework.streamed_material_state import authenticate, file_binding
from reduced_active_basis import occupation_contrast_priority
from reduced_projection_cache import SCHEMA, atomic_json
from reduced_projection_storage import ImmutableItems

FULLQ_SCHEMA = 'p254-explicit-fullq-reuse-v1'


def shared_keys():
    keys = {'initial', 'mechanical', 'kernel/00', 'passive/08', 'passive/12'}
    for slot in range(128):
        keys.update(('slot/%03d' % slot, 'preparation/%03d' % slot,
                     'edge/00/%03d' % slot))
    keys.update('passive-chunk/%02d/%02d' % (m, c)
                for m in (8, 12) for c in range(16))
    return keys


def complete_keys():
    return shared_keys() | {'kernel/%02d' % q for q in range(64)} | {
        'edge/%02d/%03d' % (q, s) for q in range(64) for s in range(128)}


def _read_manifest(path):
    binding = file_binding(path)
    if not Path(binding['path']).is_absolute():
        raise ValueError('absolute immutable manifest binding required')
    manifest = json.loads(authenticate(binding).read_text())
    if manifest.get('schema') != SCHEMA:
        raise ValueError('foreign projection manifest schema')
    return binding, manifest


def _authenticate_items(items, reader=None):
    reader = reader if reader is not None else ImmutableItems(items)
    for key, record in items.items():
        binding = record['archive']
        if not Path(binding['path']).is_absolute():
            raise ValueError('relative immutable archive binding: '+key)
        reader.authenticate(binding)
    return reader


def _history(spec, directory):
    history = []
    for family in ('implementation', 'physical_implementation'):
        for original in spec[family]:
            source = Path(original['path'])
            snapshot = file_binding(Path(directory) / (
                source.stem+'-'+original['sha256']+source.suffix))
            if (snapshot['sha256'] != original['sha256'] or
                    int(snapshot['bytes']) != int(original['bytes'])):
                raise ValueError('historical producer snapshot bytes differ: '+str(source))
            history.append(dict(family=family, original_binding=original,
                                archived_snapshot=snapshot))
    return history


def _authenticate_prior_lineage(items):
    reader = ImmutableItems({})
    for record in items.values():
        reuse = record.get('metadata', {}).get('explicit_reuse')
        if reuse is None:
            continue
        reader.authenticate(reuse['manifest'])
        for snapshot in reuse.get('original_producer_snapshots', []):
            original = snapshot['original_binding']
            archived = snapshot['archived_snapshot']
            if (original['sha256'] != archived['sha256'] or
                    int(original['bytes']) != int(archived['bytes'])):
                raise ValueError('inherited pilot producer history bytes differ')
            reader.authenticate(archived)


def _compare_specifications(old, current, gamma_binding):
    # Source-code changes are authorized ONLY through the archived old producer
    # history and explicit item lineage below. Compare scientific fields,
    # including population context, boundary, selections and feature partners.
    # The destination's explicit seed is this settled Gamma bank, irrespective
    # of whether that bank was produced fresh or reused an older pilot.
    excluded = {'implementation', 'explicit_gamma_pilot_reuse'}
    if {k: v for k, v in old.items() if k not in excluded} != {
            k: v for k, v in current.items() if k not in excluded}:
        names = sorted(k for k in set(old) | set(current)
                       if k not in excluded and old.get(k) != current.get(k))
        raise ValueError('Gamma scientific/input specification differs: '+', '.join(names))
    if old.get('bands') != 12 or old.get('features') != 27:
        raise ValueError('complete maximum12-band/27-family Gamma bank required')
    if current.get('explicit_gamma_pilot_reuse') != gamma_binding:
        raise ValueError('destination must explicitly seed from the admitted Gamma bank')
    for family in ('implementation', 'physical_implementation'):
        for binding in current[family]:
            authenticate(binding)
    authenticate(current['original_boundary'])
    authenticate(current['neutral_labels'])
    authenticate(current['explicit_gamma_pilot_reuse'])


def _validate_prepared_inputs(exporter, manifest, reader):
    """Check saved selections against actual current populations without Gram/H0.

    Only small selected arrays are decompressed; active columns, raw WFC, Gram
    and full H0 are never rebuilt. Authenticated context and boundary receipts
    bind the complete original physical inputs, not just their label strings.
    """
    source_reader = ImmutableItems({})
    for index, slot in enumerate(exporter.slots):
        key = 'slot/%03d' % index
        prep_key = 'preparation/%03d' % index
        metadata = manifest['items'][key]['metadata']
        prepared = manifest['items'][prep_key]['metadata']
        identity = str(slot.ispin)+':'+','.join(map(str, slot.mesh_index))
        original = exporter.records[identity]
        state_binding = exporter.case.source_states[(slot.ispin, slot.source_ik)]
        state = json.loads(source_reader.authenticate(state_binding).read_text())
        expected = dict(slot=index, identity=identity, frame=original.frame,
                        raw_d=original.coefficient, source=dict(slot.source_binding),
                        source_state=state_binding, source_archive=state['archive'],
                        basis_identity=original.basis_identity,
                        frame_identity=original.gauge_identity,
                        orbital_conjugation=slot.orbital_conjugation,
                        operation=slot.operation.index,
                        reciprocal_shift=slot.reciprocal_shift.tolist())
        for name, value in expected.items():
            if metadata.get(name) != value or prepared.get(name) != value:
                raise ValueError('foreign saved preparation %03d: %s' % (index, name))
        if (original.component != slot.ispin-1 or original.weight != slot.weight or
                original.source_identity != slot.source_binding['sha256']):
            raise ValueError('original sector population/measure/source differs')
        with np.load(source_reader.authenticate(state['archive']), allow_pickle=False) as archive:
            populations = archive['populations']
        labels = np.asarray(slot.eigenvalues_ry)
        priority, groups = occupation_contrast_priority(labels, exporter.interval, populations)
        selected = priority[:12]
        for item_metadata in (metadata, prepared):
            diagnostic = item_metadata['diagnostics']
            if (item_metadata.get('selected') != selected.tolist() or
                    diagnostic.get('selected_raw_labels') != selected.tolist() or
                    diagnostic.get('selected_population_reference') != populations[selected].tolist() or
                    diagnostic.get('selected_labels_ry') != labels[selected].tolist() or
                    diagnostic.get('population_groups_fractional_filled_empty') != groups or
                    diagnostic.get('neutral_interval_ry') != list(exporter.interval) or
                    diagnostic.get('no_population_refill') is not True or
                    diagnostic.get('selection_recipe') != manifest['specification']['selection_recipe']):
                raise ValueError('actual populations/selection differ for slot %03d' % index)
        fixed = reader.get(key, ('selected', 'selected_labels', 'rho0'))
        with np.load(reader.authenticate(manifest['items'][prep_key]['archive']),
                     allow_pickle=False) as archive:
            required = {'miller', 'k_cart_bohr_inv', 'reciprocal_bohr_inv', 'row_map',
                        'source_miller', 'active_columns', 'rho0', 'labels', 'selected'}
            for number, _ in enumerate(prepared['overlap_names']):
                required.update(('overlap_base_%d' % number, 'overlap_derivative_%d' % number))
            if not required.issubset(archive.files):
                raise ValueError('incomplete immutable prepared arrays: '+prep_key)
            if (not np.array_equal(fixed['selected'], selected) or
                    not np.array_equal(archive['selected'], selected) or
                    not np.array_equal(fixed['selected_labels'], labels[selected]) or
                    not np.array_equal(archive['labels'], labels[selected]) or
                    not np.array_equal(fixed['rho0'], archive['rho0'])):
                raise ValueError('slot/preparation arrays disagree: '+prep_key)
        edge = manifest['items']['edge/00/%03d' % index]['metadata']
        edge_expected = dict(source=index, target=index, q_index=0,
                             counterpart_q=0, counterpart_source=index,
                             source_binding=original.frame, target_binding=original.frame,
                             source_d=original.coefficient, target_d=original.coefficient,
                             source_basis=original.basis_identity, target_basis=original.basis_identity,
                             feature_g=manifest['specification']['feature_g'][0],
                             feature_partner=manifest['specification']['feature_partner'][0])
        if any(edge.get(name) != value for name, value in edge_expected.items()):
            raise ValueError('foreign Gamma edge metadata: %03d' % index)
        exporter.check_budget()
    for m in (8, 12):
        split = manifest['items']['passive/%02d' % m]['metadata']
        if split.get('bands') != m or split.get('frozen') is not True:
            raise ValueError('foreign passive split')
        for chunk in range(16):
            metadata = manifest['items']['passive-chunk/%02d/%02d' % (m, chunk)]['metadata']
            if metadata.get('bands') != m or metadata.get('slots') != list(range(8*chunk, 8*(chunk+1))):
                raise ValueError('foreign passive chunk')


def adopt_complete_gamma(exporter, manifest_path, producer_history):
    """Explicitly reference complete Gamma inputs; never modify the source bank.

    Call AFTER ordinary MaterialExporter construction with the ORIGINAL explicit
    reuse_gamma_pilot argument and BEFORE any source/edge/kernel calls. Existing
    MappedPreparations and ImmutableItems see the same dictionary, updated in place.
    """
    binding, old = _read_manifest(manifest_path)
    if exporter.cache.path.resolve() == Path(binding['path']):
        raise ValueError('source Gamma bank and destination worker must be disjoint')
    if old.get('status') not in ('gamma-complete', 'complete'):
        raise ValueError('settled complete128-sector Gamma bank required')
    missing = sorted(shared_keys()-old['items'].keys())
    if missing:
        raise ValueError('Gamma bank lacks required inputs: '+', '.join(missing))
    _compare_specifications(old['specification'], exporter.cache.manifest['specification'], binding)
    history = _history(old['specification'], producer_history)
    items = {key: old['items'][key] for key in sorted(shared_keys())}
    reader = _authenticate_items(items)
    _authenticate_prior_lineage(items)
    _validate_prepared_inputs(exporter, old, reader)
    authenticate(binding)
    mapped = {}
    for slot in range(128):
        receipt_path = Path(binding['path']).parent/'prepared-maps'/('%03d' % slot)/'completed.json'
        receipt_binding = file_binding(receipt_path)
        receipt = json.loads(authenticate(receipt_binding).read_bytes())
        if (receipt.get('schema') != 'p254-mapped-active-preparation-v1' or
                receipt['source_archive'] != old['items']['preparation/%03d' % slot]['archive']):
            raise ValueError('original prepared-map receipt differs: %03d' % slot)
        mapped['%03d' % slot] = receipt_binding
    descriptor = dict(schema=FULLQ_SCHEMA, gamma_manifest=binding,
                      original_producer_history=history,
                      prepared_map_receipts=mapped,
                      helper_implementation=[file_binding(Path(__file__)),
                          file_binding(Path(__file__).with_name('run_reduced_fullq_projection.py'))])
    cache = exporter.cache
    existing = cache.manifest.get('fullq')
    if existing is not None and existing != descriptor:
        raise ValueError('worker full-q adoption/helper source binding differs')
    adopted = {}
    for key, record in items.items():
        value = deepcopy(record)
        value['metadata']['fullq_item_lineage'] = dict(
            gamma_manifest=binding, original_item=key,
            producer_history='fullq.original_producer_history',
            scope='unchanged admitted Gamma product; NOT generated by current implementation')
        previous = cache.manifest['items'].get(key)
        if previous is not None and previous != value:
            # The ordinary constructor has just explicitly adopted these three
            # pilot products. Only an identical archive can be promoted into its
            # complete-Gamma lineage, never an arbitrary old worker collision.
            if (existing is not None or key not in ('initial', 'mechanical', 'kernel/00') or
                    previous['archive'] != value['archive'] or
                    'explicit_reuse' not in previous['metadata']):
                raise ValueError('Gamma adoption item collision: '+key)
        adopted[key] = value
    if existing is None and set(cache.manifest['items'])-{'initial', 'mechanical', 'kernel/00'}:
        raise ValueError('new worker already contains non-pilot products; explicit adoption refused')
    cache.manifest['items'].update(adopted)
    cache.manifest['fullq'] = descriptor
    for key in adopted:
        cache.manifest['timings'][key] = deepcopy(old.get('timings', {}).get(key, {}))
    cache.manifest['status'] = 'gamma-complete'
    cache.publish()
    return descriptor


class FullQCache:
    """Delegate immutable storage while restricting one writer to its assignment."""
    def __init__(self, cache, q_indices, sources, kernels):
        self._cache = cache
        assignment = dict(q=list(q_indices), sources=list(sources), kernels=bool(kernels))
        previous = cache.manifest.get('fullq_shard')
        if previous is not None and previous != assignment:
            raise ValueError('explicit shard assignment differs on resume')
        if cache.manifest.get('fullq_union') is not None:
            raise ValueError('a settled union cannot be resumed as a worker')
        cache.manifest['fullq_shard'] = assignment
        self.allowed = {'edge/%02d/%03d' % (q, s) for q in q_indices for s in sources}
        if kernels:
            self.allowed.update('kernel/%02d' % q for q in q_indices)
        unexpected = set(cache.manifest['items'])-shared_keys()-self.allowed
        if unexpected:
            raise ValueError('worker contains unassigned products: '+', '.join(sorted(unexpected)))
        cache.publish()

    def __getattr__(self, name):
        return getattr(self._cache, name)

    def put(self, key, arrays, *, metadata=None, timings=None):
        if key not in self.allowed or key in shared_keys():
            raise ValueError('worker may not rebuild admitted/unassigned products: '+key)
        value = dict(metadata or {})
        value['fullq_item_lineage'] = dict(
            producer='specification.implementation + specification.physical_implementation',
            helpers='fullq.helper_implementation', prepared_inputs='fullq.gamma_manifest',
            actual_prepared_maps='fullq.prepared_map_receipts',
            scope='new finite-q product of current explicitly bound producer')
        self._cache.put(key, arrays, metadata=value, timings=timings)


def _validate_worker(manifest):
    if manifest.get('fullq', {}).get('schema') != FULLQ_SCHEMA:
        raise ValueError('not an explicitly adopted full-q worker')
    if manifest.get('fullq_union') is not None:
        raise ValueError('union inputs must be independent workers, not nested unions')
    if manifest.get('fullq_worker_status') != 'settled':
        raise ValueError('worker is not settled successfully')
    if any(run.get('status') == 'running' for run in manifest.get('fullq_commands', [])):
        raise ValueError('worker still has an active/unsettled command')
    assignment = manifest['fullq_shard']
    required = {'edge/%02d/%03d' % (q, s)
                for q in assignment['q'] for s in assignment['sources']}
    if assignment['kernels']:
        required.update('kernel/%02d' % q for q in assignment['q'])
    expected = shared_keys() | required
    if set(manifest['items']) != expected:
        raise ValueError('worker missing assigned blocks or has foreign/unassigned products')
    return required-shared_keys()


def _union_manifests(manifest_paths, destination):
    """Publish ONLY a complete deterministic union of settled independent shards.

    Each non-Gamma edge/kernel has exactly one owner; even byte-identical duplicate
    nonshared records are collisions. Shared records AND timings must be identical.
    No new arrays are generated and no numerical block is filled or repaired.
    """
    inputs = sorted({str(Path(path).resolve()) for path in manifest_paths})
    if len(inputs) != len(manifest_paths) or not inputs:
        raise ValueError('distinct explicit shard manifests required')
    output = Path(destination).resolve()
    if output.suffix != '.json':
        raise ValueError('union destination must be an explicit .json manifest path')
    if str(output) in inputs or output.exists():
        raise ValueError('union destination must be new and disjoint from all source banks')
    base = None
    records, timings, owners, bindings = {}, {}, {}, []
    shared_reader = ImmutableItems({})
    for path in inputs:
        binding, manifest = _read_manifest(path)
        owned = _validate_worker(manifest)
        _authenticate_items(manifest['items'], shared_reader)
        if base is None:
            base = manifest
        elif (manifest['specification'] != base['specification'] or
              manifest['fullq'] != base['fullq']):
            raise ValueError('shards do not share the exact scientific/source/reuse specification')
        for key, item in sorted(manifest['items'].items()):
            timing = manifest.get('timings', {}).get(key)
            if key in records:
                if key not in shared_keys():
                    raise ValueError('independent nonshared ownership collision: '+key)
                if records[key] != item or timings.get(key) != timing:
                    raise ValueError('shared item/timing records differ: '+key)
            else:
                records[key] = item
                if timing is not None:
                    timings[key] = timing
            if key in owned:
                owners[key] = binding
        bindings.append(binding)
    missing = sorted(complete_keys()-records.keys())
    foreign = sorted(records.keys()-complete_keys())
    if missing or foreign:
        raise ValueError('union refuses incomplete/foreign coverage: '+json.dumps(
            dict(missing=missing, foreign=foreign), sort_keys=True))
    for binding in base['specification']['implementation']+base['specification']['physical_implementation']:
        authenticate(binding)
    authenticate(base['fullq']['gamma_manifest'])
    for record in base['fullq']['original_producer_history']:
        authenticate(record['archived_snapshot'])
    for binding in base['fullq']['helper_implementation']:
        authenticate(binding)
    result = dict(schema=SCHEMA, status='complete', specification=base['specification'],
                  fullq=base['fullq'], items=dict(sorted(records.items())),
                  timings=dict(sorted(timings.items())), failures=[],
                  fullq_union=dict(schema=FULLQ_SCHEMA, source_manifests=bindings,
                                   item_owners=dict(sorted(owners.items())),
                                   coverage=dict(q=64, slots=128, passive_bands=[8, 12]),
                                   source_commands_preserved_in_original_manifests=True))
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output, result)
    return output


def union_manifests(manifest_paths, destination):
    """Hold read ownership of every settled worker until union publication."""
    paths = sorted({Path(path).resolve() for path in manifest_paths})
    if len(paths) != len(manifest_paths) or not paths:
        raise ValueError('distinct explicit shard manifests required')
    output = Path(destination).resolve()
    if any(output.parent == path.parent for path in paths):
        raise ValueError('union output must be in a disjoint destination directory')
    with ExitStack() as stack:
        for directory in sorted({path.parent for path in paths}):
            lock = stack.enter_context((directory/'.fullq-writer.lock').open('r'))
            try:
                fcntl.flock(lock.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError('union input still has an active writer: '+str(directory)) from exc
        return _union_manifests(manifest_paths, destination)
