#!/usr/bin/env python3
"""Authenticate committed A35 sectors and sum raw per-cell current covariance.

All nine directed complex pairs and RR/RE/ER/EE blocks survive unchanged except
for the original 1/64 quadrature applied once. Partial sums are not normalized.
No physical actions or work launch on import. Exit: full128=0, invalid=1, partial=2.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

A34 = Path(__file__).resolve().parent.parent / '0034'
sys.path.insert(0, str(A34))
from aggregate_mode_noise import (authenticate, binding, complex_record, complex_sum,
                                 complex_value, load_json, number, same_binding, save_json)
from export_uniform_current_noise import scientific_identity

PROBES = ('current_x', 'current_y', 'current_z')
PARTS = ('retained', 'omitted')
UNITS = 'electronic_current_hartree'
RESULT_SCHEMA = 'P254-A35-original-sector-uniform-current-noise-v1'
AGGREGATE_SCHEMA = 'P254-A35-initial-Gamma-uniform-current-noise-aggregate-v1'
BATCH_SCHEMA = 'P254-A35-uniform-current-noise-batch-summary-v1'
SCOPE = ('Initial externalGamma raw physical per-cell electronic-current covariance, J=-L/2; '
         'not canonical photon generalized-force noise, finite-q density, stationarity, '
         'temperature/FDT, loss, linewidth, heat, trajectory or full parent closure.')


def parse_covariance(rows):
    if not isinstance(rows, list) or len(rows) != 9:
        raise ValueError('all nine directed current covariance rows required')
    full, blocks, defects = {}, {}, {}
    for row in rows:
        pair = (row['left'], row['right'])
        if pair[0] not in PROBES or pair[1] not in PROBES or pair in full:
            raise ValueError('duplicate or unknown directed current pair')
        full[pair] = complex_value(row['full'])
        defects[pair] = complex_value(row['reconstruction_defect'])
        if not isinstance(row['blocks'], list) or len(row['blocks']) != 4:
            raise ValueError('all four RR/RE/ER/EE blocks required')
        for block in row['blocks']:
            key = pair + (block['part_left'], block['part_right'])
            if key[2] not in PARTS or key[3] not in PARTS or key in blocks:
                raise ValueError('duplicate or unknown retained/omitted block')
            blocks[key] = complex_value(block['value'])
    return full, blocks, defects


def slot_metadata(slot):
    return {key: slot[key] for key in ('identity', 'component', 'weight', 'shape',
                                      'basis_identity', 'gauge_identity', 'source_identity')}


def load_context(initial_result, coefficient_manifest, contract_path):
    fixed = Path(__file__).with_name('implementation-contract-001.json')
    contract_binding = binding(contract_path)
    if not same_binding(contract_binding, binding(fixed)):
        raise ValueError('fixed adjacent A35 implementation contract required')
    identity = scientific_identity()
    contract = load_json(contract_path)
    if contract['schema'] != 'P254-A35-original-full-PW-uniform-current-noise-plan-v1':
        raise ValueError('wrong A35 implementation contract')
    inputs = dict(initial_result=binding(initial_result),
                  coefficient_manifest=binding(coefficient_manifest))
    for key, authority in (('initial_result', 'original_initial'),
                           ('coefficient_manifest', 'original_coefficients')):
        if not same_binding(inputs[key], contract[authority]):
            raise ValueError('input differs from original A35 authority: ' + key)
    initial = load_json(initial_result)
    if (initial['schema'] != 'P254-matched-joint-case-v1' or initial['status'] != 'complete'
            or initial['mode'] != 'initialize' or initial['material_trajectory'] is not False):
        raise ValueError('complete original matched initialization required')
    saved = initial['saved']
    original = authenticate(saved['original_boundary'])
    preparation = authenticate(saved['preparation'])
    if (not same_binding(original, contract['original_boundary'])
            or not same_binding(preparation, contract['preparation'])):
        raise ValueError('original boundary/preparation differs from fixed A35 authority')
    boundary = load_json(original['path'])
    slots = boundary['slots']
    if (boundary['schema'] != 'streamed-material-mean-v1' or len(slots) != 128
            or boundary['time_ry_inverse'] != 0 or any(boundary['u_bohr'])
            or any(boundary['p_canonical'])
            or boundary['model_identity'] != saved['original_model_identity']
            or boundary['evidence']['context_identity'] != saved['model_identity']
            or boundary['evidence']['case'] != initial['case']
            or len({s['identity'] for s in slots}) != 128):
        raise ValueError('literal original complete128 t=0 boundary required')
    manifest = load_json(coefficient_manifest)
    specification = manifest['specification']
    if (manifest['schema'] != 'p254-physical-activeband-cache-v1'
            or manifest['status'] != 'complete' or specification['bands'] != 12
            or not same_binding(specification['original_boundary'], original)):
        raise ValueError('matching original saved12 coefficient manifest required')
    physical_spec = saved['specification']
    for name in ('case', 'context', 'aggregate', 'source_states', 'scalar', 'preparation'):
        if specification['context'][name] != physical_spec[name]:
            raise ValueError('saved coefficient/preparation context differs: ' + name)
    selected = {}
    for key, item in manifest['items'].items():
        if not key.startswith('slot/'):
            continue
        name = item['metadata']['identity']
        if name in selected:
            raise ValueError('duplicate coefficient sector identity')
        selected[name] = (key, item)
    if set(selected) != {s['identity'] for s in slots}:
        raise ValueError('coefficient manifest must cover original128 exactly')
    sectors = set()
    for slot in slots:
        if number(slot['weight'], 'quadrature') != 1 / 64 or slot['shape'][1] != 346:
            raise ValueError('literal 1/64 weight and all346 raw columns required')
        spin_text, mesh_text = slot['identity'].split(':')
        spin, mesh = int(spin_text), tuple(int(i) for i in mesh_text.split(','))
        if (spin not in (1, 2) or len(mesh) != 3 or any(i not in range(4) for i in mesh)
                or slot['component'] != spin - 1):
            raise ValueError('literal explicit spin/internal4x4x4 identity required')
        sectors.add((spin, mesh))
        _, item = selected[slot['identity']]
        metadata = item['metadata']
        indices = metadata['selected']
        if (metadata['basis_identity'] != slot['basis_identity']
                or metadata['frame_identity'] != slot['gauge_identity']
                or not same_binding(metadata['frame'], slot['frame'])
                or not same_binding(metadata['raw_d'], slot['coefficient'])
                or len(indices) != 12 or any(type(i) is not int for i in indices)
                or len(set(indices)) != 12 or min(indices) < 0 or max(indices) >= 346):
            raise ValueError('saved12 selection differs from original frame/D/basis')
    if len(sectors) != 128:
        raise ValueError('literal original128 explicit spin/internal-k coverage required')
    return dict(contract=contract, contract_binding=contract_binding, inputs=inputs,
                initial=initial, boundary=boundary, slots=slots, selected=selected,
                original_boundary=original, preparation=preparation,
                coefficient_case_identity=specification['case_identity'],
                selection_recipe=specification['selection_recipe'],
                scientific_compatibility_identity=identity)


def load_case(done_path, context, expected_record=None):
    """Admit only exact producer-identity, original-source committed cases."""
    done_path = Path(done_path).resolve(strict=True)
    if done_path.name != 'done.json':
        raise ValueError('case completion must be named done.json')
    done_binding = binding(done_path)
    done = load_json(done_path)
    result_binding = authenticate(done['result'], done_path.parent / 'result.json')
    result = load_json(result_binding['path'])
    if result['schema'] != RESULT_SCHEMA or result['status'] != 'complete':
        raise ValueError('complete A35 original-sector result required')
    record = result['original_record']
    if type(record) is not int or record not in range(128):
        raise ValueError('exact original record in 0..127 required')
    if expected_record is not None and record != expected_record:
        raise ValueError('completion belongs to another original record')
    if result['scientific_compatibility_identity'] != context['scientific_compatibility_identity']:
        raise ValueError('current A35 scientific source identity required; no legacy admission')
    if (result['quadrature_applied'] is not False or result['probe_names'] != list(PROBES)
            or result['units'] != UNITS):
        raise ValueError('unweighted exact Cartesian electronic-current probes required')
    slot = context['slots'][record]
    if number(result['quadrature_weight'], 'quadrature') != slot['weight']:
        raise ValueError('quadrature differs from original sector')
    source = result['source']
    if type(source['original_record']) is not int or source['original_record'] != record:
        raise ValueError('source original record mismatch')
    for key, expected in (('initial_result', context['inputs']['initial_result']),
                          ('coefficient_manifest', context['inputs']['coefficient_manifest']),
                          ('original_boundary', context['original_boundary']),
                          ('preparation', context['preparation']),
                          ('original_raw_d', slot['coefficient'])):
        if not same_binding(source[key], expected):
            raise ValueError('source authority differs: ' + key)
    if (not same_binding(result['frame'], slot['frame'])
            or result['pw_basis_identity'] != slot['basis_identity']
            or source['original_metadata'] != slot_metadata(slot)):
        raise ValueError('original full frame/basis/metadata differs')
    key, item = context['selected'][slot['identity']]
    metadata = item['metadata']
    if (source['coefficient_item'] != key
            or not same_binding(source['coefficient_archive'], item['archive'])
            or source['selected'] != metadata['selected']
            or source['coefficient_case_identity'] != context['coefficient_case_identity']
            or source['selection_recipe'] != context['selection_recipe']):
        raise ValueError('source saved12 selection lineage differs')
    if (source['case'] != context['initial']['case']
            or source['model_identity'] != context['initial']['saved']['model_identity']
            or source['original_model_identity'] != context['boundary']['model_identity']):
        raise ValueError('source model/preparation mismatch')
    if type(source['spin']) is not int or source['spin'] not in (1, 2):
        raise ValueError('literal explicit spin1/2 required')
    mesh = source['mesh_index']
    if not isinstance(mesh, list) or len(mesh) != 3 or any(type(i) is not int for i in mesh):
        raise ValueError('exact three-index internal mesh required')
    if (str(source['spin']) + ':' + ','.join(str(i) for i in mesh) != slot['identity']
            or source['spin'] - 1 != slot['component']):
        raise ValueError('spin/internal mesh differs from original record identity')
    for field in ('source', 'source_state', 'source_archive'):
        if not same_binding(source[field], metadata[field]):
            raise ValueError('saved source lineage differs: ' + field)
    for field in ('orbital_conjugation', 'operation', 'reciprocal_shift', 'source_admission'):
        if source[field] != metadata[field]:
            raise ValueError('saved source transport differs: ' + field)
    physical_spec = context['initial']['saved']['specification']
    expected_state = physical_spec['source_states'][str((source['spin'], source['source_ik']))]
    if not same_binding(source['source_state'], expected_state):
        raise ValueError('source state differs from original physical preparation')
    authenticate(source['source_state'])
    state = load_json(source['source_state']['path'])
    if not same_binding(state['archive'], source['source_archive']):
        raise ValueError('source-state archive differs')
    authority = dict(initial=context['inputs']['initial_result'], preparation=source['preparation'],
                     boundary=source['original_boundary'], record=slot_metadata(slot),
                     frame=slot['frame'], raw_d=slot['coefficient'],
                     state=source['source_state'], archive=state['archive'])
    digest = hashlib.sha256(json.dumps(authority, sort_keys=True, allow_nan=False).encode()).hexdigest()
    if result['state_authority'] != digest:
        raise ValueError('original full-state authority differs')
    work = result['work']
    for field in ('elapsed_seconds', 'contraction_seconds', 'peak_rss_bytes'):
        number(work[field], field, nonnegative=True)
    expected_work = dict(current_action_evaluations=3, adjoint_reuses_by_hermiticity=3,
                         force_forward_actions=0, force_adjoint_actions=0, h0_actions=0,
                         trajectory_actions=0, native_calls=0, gram_recomputations=0,
                         population_refills=0, preparation_replays=0)
    for field, expected in expected_work.items():
        if type(work['operator_work'][field]) is not int or work['operator_work'][field] != expected:
            raise ValueError('current-only work contract differs: ' + field)
    if work['raw_action_payloads_written'] is not False:
        raise ValueError('raw action arrays must not be checkpointed')
    return dict(original_record=record, done=done_binding, result=result_binding,
                payload=result, parsed=parse_covariance(result['covariance']),
                scientific_admission='producer-scientific-identity')


def aggregate(cases, context):
    by_record, identities = {}, set()
    for case in cases:
        record = case['original_record']
        source = case['payload']['source']
        identity = (source['spin'], tuple(source['mesh_index']))
        if record in by_record or identity in identities:
            raise ValueError('duplicate original record or explicit spin/internal-k sector')
        by_record[record] = case
        identities.add(identity)
    present = sorted(by_record)
    missing = sorted(set(range(128)) - set(present))
    rows = []
    for a in PROBES:
        for b in PROBES:
            pair = (a, b)
            full = complex_sum(c['payload']['quadrature_weight'] * c['parsed'][0][pair]
                               for c in by_record.values())
            blocks = []
            for pa in PARTS:
                for pb in PARTS:
                    value = complex_sum(c['payload']['quadrature_weight']
                                        * c['parsed'][1][pair + (pa, pb)]
                                        for c in by_record.values())
                    blocks.append(dict(part_left=pa, part_right=pb, value=complex_record(value)))
            reconstructed = complex_sum(complex_value(block['value']) for block in blocks)
            defect = complex_sum(c['payload']['quadrature_weight'] * c['parsed'][2][pair]
                                 for c in by_record.values())
            rows.append(dict(left=a, right=b, full=complex_record(full), blocks=blocks,
                             block_reconstruction=complex_record(reconstructed),
                             reconstruction_defect=complex_record(full - reconstructed),
                             weighted_raw_reconstruction_defect=complex_record(defect)))
    return dict(schema=AGGREGATE_SCHEMA, status='partial' if missing else 'complete',
                original_records=present, missing_original_records=missing,
                required_original_records=list(range(128)), sector_count=len(present),
                contract=context['contract_binding'], inputs=context['inputs'],
                scientific_compatibility_identity=context['scientific_compatibility_identity'],
                original_boundary=context['original_boundary'], preparation=context['preparation'],
                probe_names=list(PROBES), units=UNITS, covariance_units=UNITS + '^2',
                normalization=dict(rule='sum_k actual_quadrature_weight_k * unweighted_C_k',
                                   quadrature_applied_once=True, quadrature_weight=1 / 64,
                                   required_explicit_spin_sectors=128, extra_spin_factor=1,
                                   final_spatial_factor=1, partial_renormalization=False,
                                   canonical_field_coupling_applied=False),
                covariance=rows, numerical_repairs_applied=False,
                sectors=[dict(original_record=r, done=by_record[r]['done'],
                              result=by_record[r]['result'],
                              scientific_admission=by_record[r]['scientific_admission'],
                              quadrature_weight=by_record[r]['payload']['quadrature_weight'],
                              source=by_record[r]['payload']['source'],
                              work=by_record[r]['payload']['work'],
                              raw_unweighted_covariance=by_record[r]['payload']['covariance'])
                         for r in present], scope=SCOPE)


def add_input_arguments(parser):
    parser.add_argument('--initial-result', type=Path, required=True)
    parser.add_argument('--coefficient-manifest', type=Path, required=True)
    parser.add_argument('--contract', type=Path,
                        default=Path(__file__).with_name('implementation-contract-001.json'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_input_arguments(parser)
    parser.add_argument('--done', type=Path, action='append', default=[])
    parser.add_argument('--batch-summary', type=Path, action='append', default=[])
    parser.add_argument('--output', type=Path, required=True,
                        help='new exclusive JSON path; parent directory must exist')
    args = parser.parse_args()
    try:
        context = load_context(args.initial_result, args.coefficient_manifest, args.contract)
        paths = list(args.done)
        for path in args.batch_summary:
            summary = load_json(path)
            if (summary['schema'] != BATCH_SCHEMA or summary['scientific_compatibility_identity']
                    != context['scientific_compatibility_identity']):
                raise ValueError('wrong batch summary schema or scientific identity')
            for item in summary['completed']:
                authenticate(item['done'])
                paths.append(Path(item['done']['path']))
        if not paths:
            raise ValueError('supply explicit committed --done or --batch-summary cases')
        result = aggregate([load_case(path, context) for path in paths], context)
        code = 0 if result['status'] == 'complete' else 2
    except Exception as error:
        result = dict(schema=AGGREGATE_SCHEMA, status='invalid', error=repr(error), scope=SCOPE)
        code = 1
    receipt = save_json(args.output, result)
    print(json.dumps(dict(event='aggregate-exit', status=result['status'], output=receipt,
                          exit_code=code), allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
