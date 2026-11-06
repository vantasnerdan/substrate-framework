"""Finite receipt joins, separate from independent real-space numerical oracles."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from substrate_framework.joint_spatial_readout import coherent_local
from substrate_framework.joint_response_observables import spatial_pair
from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import file_binding
from joint_spatial_provenance import authenticate, canonical
from joint_spatial_receipts import read_entries, validate_records
from joint_spatial_committed_oracle import exercise_committed


def exercise_receipts(*, q, q_force, ab, ba, local_force, cells, t, s, mesh,
                      expected_ab, expected_ba, encode, close, reject):
    n = len(q)
    # Temporary schema fixtures exercise actual bound-json and scientific-join
    # consumer paths. They are NOT production outputs or material trajectories.
    with TemporaryDirectory(prefix='p254-spatial-oracle-') as temporary:
        directory = Path(temporary)

        def save(name, value):
            path = directory/name
            path.write_text(json.dumps(encode(value), allow_nan=False)+'\n')
            return file_binding(path)

        source = save('source.json', {'finite_oracle_only': True})
        initial = save('initial.json', {'finite_oracle_only': True})
        aggregate = save('aggregate.json', {'finite_oracle_only': True})
        baseline = save('baseline.json', {'finite_oracle_only': True})
        specification = dict(scalar='finite-independent-real-space', aggregate=aggregate)
        saved = dict(case='finite-oracle-only',
                     model_identity=sha256(canonical(specification).encode()).hexdigest(),
                     specification=specification)
        records, requests, all_entries = [], [], []
        for i in range(n):
            identity = sha256(canonical(dict(model=saved['model_identity'], q=q[i].tolist(),
                dynamics='screened-joint', sources=[source])).encode()).hexdigest()
            quantum = dict(AB=dict(zip(('electronic', 'ionic', 'field', 'total'), ab[i])),
                BA=dict(zip(('electronic', 'ionic', 'field', 'total'), ba[i])),
                symmetrized=(ab[i, 3]+ba[i, 3])/2, commutator=ab[i, 3]-ba[i, 3],
                time=t, second_time=s, dynamics='screened-joint', propagation_identity=identity,
                preparation_identity='original-matched-product-'+saved['model_identity'])
            record = dict(schema='P254-joint-q-material-v2', status='complete', q=q[i].tolist(),
                required_q_for_whole_local=n, fixed_original_D=True, all_explicit_edges=2*n,
                case=saved['case'], model_identity=saved['model_identity'], scalar=saved['specification']['scalar'],
                matched_aggregate=aggregate, initialization=initial, baseline=baseline,
                dynamics='screened-joint', source_bindings=[source], time_grid=[0., t],
                second_time_grid=[0., s], observable='original-H1x-total-force', quantum=quantum,
                coherent=dict(readout=q_force[i], amplitude=1.7, cell=[1, 0, 0]))
            request = dict(identity=identity, mode='all', case=saved['case'], initialization=initial,
                baseline=baseline, q=q[i].tolist(), times=[0., t], second_times=[0., s],
                atom=1, axis='x', amplitude=1.7, cell=[1, 0, 0], subdivisions=1, block_size=8,
                checkpoint='None', checkpoint_sha256=None)
            entry = dict(result=save('result-%d.json' % i, record),
                         request=save('request-%d.json' % i, request))
            all_entries.append(entry)
            records.append(bound_json(entry['result']))
            requests.append(bound_json(entry['request']))
        reject('missing_q', lambda: validate_records(records[:-1], requests[:-1], saved, initial, mesh))
        duplicate = deepcopy(records)
        duplicate[-1]['q'] = duplicate[0]['q']
        reject('duplicate_q', lambda: validate_records(duplicate, requests, saved, initial, mesh))
        offgrid = deepcopy(records)
        offgrid[-1]['q'][0] = .17
        reject('off_grid_q', lambda: validate_records(offgrid, requests, saved, initial, mesh))
        changes = {
            'mixed_baseline': lambda r: r.update(baseline=source),
            'mixed_case': lambda r: r.update(case='other'),
            'mixed_observable': lambda r: r.update(observable='original-H1y-total-force'),
            'mixed_preparation': lambda r: r['quantum'].update(preparation_identity='other'),
            'mixed_dynamics': lambda r: r.update(dynamics='bare-electronic'),
            'mixed_time': lambda r: r['quantum'].update(time=t+.01),
            'mixed_sources': lambda r: r.update(source_bindings=[initial]),
            'missing_origin': lambda r: r['quantum']['AB'].pop('field'),
            'inconsistent_origin_sum': lambda r: r['quantum']['AB'].update(total=99.),
            'incomplete_status': lambda r: r.update(status='partial')}
        for name, change in changes.items():
            changed = deepcopy(records)
            change(changed[-1])
            reject(name, lambda: validate_records(changed, requests, saved, initial, mesh))

        def consume(entries):
            forward, forward_requests, pairs, pair_requests = read_entries(entries)
            return validate_records(forward, forward_requests, saved, initial, mesh,
                                    pair_records=pairs, pair_requests=pair_requests)

        split_entries = []
        for i, (record, request) in enumerate(zip(records, requests)):
            forward, pair = deepcopy(record), deepcopy(record)
            del forward['quantum']
            del pair['coherent']
            forward_request, pair_request = dict(request, mode='forward'), dict(request, mode='pair')
            split_entries.append(dict(
                forward=dict(result=save('forward-result-%d.json' % i, forward),
                             request=save('forward-request-%d.json' % i, forward_request)),
                pair=dict(result=save('pair-result-%d.json' % i, pair),
                          request=save('pair-request-%d.json' % i, pair_request))))
        for label, entries in (('all_manifest', all_entries), ('split_manifest', split_entries),
                               ('mixed_forms_manifest', all_entries[:2]+split_entries[2:])):
            selected, f_read, a_read, b_read = consume(entries)
            close(label+'_coherent_dense', coherent_local(q, f_read, cells), local_force)
            close(label+'_AB_dense', [spatial_pair(q, a_read, cell, [1, 0, 0], ncells=n)
                                     for cell in cells], expected_ab)
            close(label+'_BA_dense', [spatial_pair(q, b_read, cell, [1, 0, 0], ncells=n)
                                     for cell in cells], expected_ba)
            if selected != validate_records(records, requests, saved, initial, mesh)[0]:
                raise AssertionError(label+' changed scientific comparison identity')
        exercise_committed(save=save, split_entries=split_entries, specification=specification,
            consume=consume, q=q, cells=cells, local_force=local_force,
            expected_ab=expected_ab, expected_ba=expected_ba, close=close, reject=reject)
        split_changes = (
            'source', 'baseline', 'time', 'second_time', 'seed_amplitude', 'seed_cell',
            'seed_atom', 'seed_axis', 'preparation', 'q', 'scalar', 'observable', 'mode')
        for name in split_changes:
            entries = deepcopy(split_entries)
            pair = bound_json(entries[-1]['pair']['result'])
            request = bound_json(entries[-1]['pair']['request'])
            if name == 'source':
                pair['source_bindings'] = [initial]
            elif name == 'baseline':
                pair['baseline'] = request['baseline'] = source
            elif name in ('time', 'second_time'):
                record_key, request_key = ('time_grid', 'times') if name == 'time' else ('second_time_grid', 'second_times')
                new_time = pair['quantum'][name]+.01
                pair[record_key] = request[request_key] = [0., new_time]
                pair['quantum'][name] = new_time
            elif name.startswith('seed_'):
                key = name.removeprefix('seed_')
                request[key] = dict(amplitude=1.9, cell=[2, 0, 0], atom=2, axis='y')[key]
            elif name == 'preparation':
                pair['quantum']['preparation_identity'] = 'different-original-preparation'
            elif name == 'q':
                pair['q'] = request['q'] = [0., 0., 0.]
            elif name == 'scalar':
                pair['scalar'] = 'different-scalar'
            elif name == 'observable':
                pair['observable'], request['axis'] = 'original-H1y-total-force', 'y'
            elif name == 'mode':
                request['mode'] = 'forward'
            identity = sha256(canonical(dict(model=pair['model_identity'], q=pair['q'],
                dynamics=pair['dynamics'], sources=pair['source_bindings'])).encode()).hexdigest()
            request['identity'] = pair['quantum']['propagation_identity'] = identity
            entries[-1]['pair'] = dict(result=save('mismatch-'+name+'-result.json', pair),
                                       request=save('mismatch-'+name+'-request.json', request))
            reject('split_mismatched_'+name, lambda: consume(entries))
        reject('split_forward_alone', lambda: consume([split_entries[0]['forward']]+split_entries[1:]))
        reject('split_missing_pair_binding', lambda: consume([dict(forward=split_entries[0]['forward'])]+split_entries[1:]))
        Path(split_entries[0]['pair']['result']['path']).write_text('changed split bytes\n')
        reject('split_changed_bound_result', lambda: consume(split_entries))
        Path(source['path']).write_text('changed bytes\n')
        reject('changed_bound_source', lambda: authenticate(source))
        reject('changed_bound_json', lambda: bound_json(source))
