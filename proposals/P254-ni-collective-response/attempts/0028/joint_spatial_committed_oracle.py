"""Temporary native-schema fixtures; no material trajectories or results."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path

from substrate_framework.joint_response_observables import spatial_pair
from substrate_framework.joint_spatial_readout import coherent_local
from substrate_framework.material_control_reference import bound_json
from joint_spatial_provenance import COMMITTED_SCHEMA, authenticate, canonical

HERE = Path(__file__).resolve().parent


def _native_fixture(save, forward, request_binding, worker, index, *, physical=None, scalar=None):
    request = bound_json(request_binding)
    sources = forward['source_bindings']+[worker]
    consumer = sha256(canonical(dict(model=forward['model_identity'], q=forward['q'],
        dynamics=forward['dynamics'], sources=sources)).encode()).hexdigest()
    provenance = dict(baseline=forward['baseline'], history=[])
    checkpoint = save('native-checkpoint-%d.json' % index,
        dict(schema='streamed-joint-response-v1', complete=True, kind='forward',
             identity=request['identity'], q=forward['q'], time=request['times'][-1],
             provenance=provenance))
    if physical is None:
        physical = save('native-physical-proof-%d.json' % index,
            dict(schema='P254-joint-partitioned-q-oracle-v1', status='passed', finite_fixture_only=True))
    if scalar is None:
        scalar = save('native-readout-proof-%d.json' % index,
            dict(schema='P254-joint-q-readout-oracle-v1', status='passed', finite_fixture_only=True))
    operator = dict(case=forward['case'], initialization=forward['initialization'],
        model_identity=forward['model_identity'], specification=forward['finite_specification'],
        history=forward['baseline'], q=forward['q'], time=request['times'][-1],
        atom=request['atom'], axis=request['axis'], observable='total-force', field_axis=None,
        probe_identity=forward['observable'], observable_metadata=None,
        block_size=request['block_size'], dynamics=forward['dynamics'], sources=sources)
    consumer_request = dict(schema=COMMITTED_SCHEMA, checkpoint=checkpoint,
        historical_producer_identity=request['identity'], historical_producer_request=request_binding,
        consumer_identity=consumer, physical_operator_oracle=physical, readout_oracle=scalar,
        operator=operator)
    return dict(schema=COMMITTED_SCHEMA, status='complete', complete=True,
        historical_producer_identity=request['identity'], consumer_identity=consumer,
        historical_producer_request=dict(binding=request_binding, metadata=request),
        committed_checkpoint=checkpoint, checkpoint_provenance=provenance,
        physical_operator_oracle=physical, readout_oracle=scalar, request=consumer_request,
        case=forward['case'], initialization=forward['initialization'], baseline=forward['baseline'],
        model_identity=forward['model_identity'], scalar=forward['scalar'],
        matched_aggregate=forward['matched_aggregate'], q=forward['q'],
        atom=request['atom'], axis=request['axis'], dynamics=forward['dynamics'],
        observable=forward['observable'], observable_selection='total-force',
        units='Ry/bohr', observable_metadata=None, time_ry_inverse=request['times'][-1],
        all_explicit_edges=forward['all_explicit_edges'], fixed_original_D=True,
        source_bindings=sources, readout=forward['coherent']['readout'])


def exercise_committed(*, save, split_entries, specification, consume, q, cells,
                       local_force, expected_ab, expected_ba, close, reject):
    snapshot = save('finite-old-readout-source.json', {'finite_original_worker': True})
    worker = dict(snapshot, path=str(HERE/'partitioned_material_readout.py'))
    snapshots = [dict(original=worker, snapshot=snapshot)]
    reject('native_old_source_without_snapshot', lambda: authenticate(worker))
    authenticate(worker, snapshots)
    wrong_snapshot = save('wrong-readout-source.json', {'different_source_bytes': True})
    reject('native_wrong_snapshot_bytes', lambda: authenticate(worker,
        [dict(original=worker, snapshot=wrong_snapshot)]))
    entries, natives = deepcopy(split_entries), []
    for index, entry in enumerate(entries):
        forward = bound_json(entry['forward']['result'])
        forward['finite_specification'] = specification
        native = _native_fixture(save, forward, entry['forward']['request'], worker, index)
        if native['consumer_identity'] == native['historical_producer_identity']:
            raise AssertionError('native fixture does not distinguish consumer from producer')
        entry['forward']['result'] = save('native-result-%d.json' % index, native)
        natives.append(native)
    for label, selected in (('native_manifest', entries),
                            ('native_legacy_mixed_manifest', entries[:2]+split_entries[2:])):
        _, force, ab, ba = consume(selected)
        close(label+'_coherent_dense', coherent_local(q, force, cells), local_force)
        for ordering, values, expected in (('AB', ab, expected_ab), ('BA', ba, expected_ba)):
            close(label+'_'+ordering+'_dense',
                [spatial_pair(q, values, cell, [1, 0, 0], ncells=len(q)) for cell in cells], expected)
    for name in ('consumer_as_producer', 'wrong_consumer', 'wrong_checkpoint_time',
                 'dual_checkpoint', 'changed_physical_source', 'changed_original_request'):
        changed, native = deepcopy(entries), deepcopy(natives[-1])
        if name == 'consumer_as_producer':
            native['historical_producer_identity'] = native['consumer_identity']
        elif name == 'wrong_consumer':
            native['consumer_identity'] = 'different-consumer'
        elif name == 'wrong_checkpoint_time':
            native['time_ry_inverse'] += .01
            native['request']['operator']['time'] = native['time_ry_inverse']
        elif name == 'dual_checkpoint':
            checkpoint = bound_json(native['committed_checkpoint'])
            checkpoint['kind'] = 'dual'
            binding = save('wrong-dual-checkpoint.json', checkpoint)
            native['committed_checkpoint'] = native['request']['checkpoint'] = binding
        elif name == 'changed_physical_source':
            sources = [wrong_snapshot, worker]
            native['source_bindings'] = native['request']['operator']['sources'] = sources
            identity = sha256(canonical(dict(model=native['model_identity'], q=native['q'],
                dynamics=native['dynamics'], sources=sources)).encode()).hexdigest()
            native['consumer_identity'] = native['request']['consumer_identity'] = identity
        elif name == 'changed_original_request':
            original = bound_json(changed[-1]['forward']['request'])
            changed[-1]['forward']['request'] = save('unbound-original-request-copy.json', original)
        changed[-1]['forward']['result'] = save('native-mismatch-'+name+'.json', native)
        reject('native_'+name, lambda: consume(changed))
