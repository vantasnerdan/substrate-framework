#!/usr/bin/env python3
"""Authenticate and sum A34 unweighted initial-Gamma force covariances.

A complete result requires the literal original records 0..127. Partial sums are
labelled partial, never renormalized. This is not late Ni noise, a bath bound,
64-external-q coverage, damping, heat, or irreversible loss. No PSD repair,
symmetrization, spin multiplier, or final spatial normalization is performed.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from mode_noise_compatibility import admit_result, load_compatibility

RESULT_SCHEMA = 'P254-A34-original-sector-two-mode-covariance-v1'
MODES = ('mode0', 'mode1')
PARTS = ('retained', 'omitted')
SCOPE = ('Initial Gamma force noise in two fixed paid output-directed modes; '
         'not late Ni output, bath/error bound, 64 external q, damping, heat or loss.')


def _unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError('duplicate JSON key: ' + key)
        value[key] = item
    return value


def load_json(path):
    with Path(path).open(encoding='utf-8') as stream:
        return json.load(stream, object_pairs_hook=_unique_object,
                         parse_constant=lambda token: _nonfinite(token))


def _nonfinite(token):
    raise ValueError('nonfinite JSON number: ' + token)


def binding(path):
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    count = 0
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b''):
            digest.update(chunk)
            count += len(chunk)
    if path.stat().st_size != count:
        raise ValueError('file size changed while reading: ' + str(path))
    return dict(path=str(path), bytes=count, sha256=digest.hexdigest())


def same_binding(left, right):
    return (Path(left['path']).resolve() == Path(right['path']).resolve()
            and int(left['bytes']) == int(right['bytes'])
            and left['sha256'] == right['sha256'])


def authenticate(record, expected_path=None):
    if expected_path is not None and Path(record['path']).resolve() != Path(expected_path).resolve():
        raise ValueError('receipt binds a different result path')
    actual = binding(record['path'])
    if not same_binding(actual, record):
        raise ValueError('file byte/SHA binding differs: ' + str(record['path']))
    return actual


def save_json(path, value):
    """Exclusive publication with data and directory fsync; never overwrite."""
    path = Path(path)
    data = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    with path.open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return binding(path)


def number(value, label, nonnegative=False):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError('finite numeric ' + label + ' required')
    if nonnegative and value < 0:
        raise ValueError('nonnegative ' + label + ' required')
    return value


def complex_value(value):
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        raise ValueError('exact real/imag complex record required')
    return complex(number(value['real'], 'real'), number(value['imag'], 'imag'))


def complex_record(value):
    return dict(real=value.real, imag=value.imag)


def complex_sum(values):
    values = list(values)
    return complex(math.fsum(v.real for v in values), math.fsum(v.imag for v in values))


def parse_covariance(rows):
    if not isinstance(rows, list) or len(rows) != 4:
        raise ValueError('all four ordered 2x2 covariance rows required')
    full, blocks, defects = {}, {}, {}
    for row in rows:
        pair = (row['left'], row['right'])
        if pair[0] not in MODES or pair[1] not in MODES or pair in full:
            raise ValueError('duplicate or unknown ordered mode row')
        full[pair] = complex_value(row['full'])
        defects[pair] = complex_value(row['reconstruction_defect'])
        if not isinstance(row['blocks'], list) or len(row['blocks']) != 4:
            raise ValueError('all four RR/RE/ER/EE blocks required per row')
        for block in row['blocks']:
            key = pair + (block['part_left'], block['part_right'])
            if key[2] not in PARTS or key[3] not in PARTS or key in blocks:
                raise ValueError('duplicate or unknown retained/omitted block')
            blocks[key] = complex_value(block['value'])
    return full, blocks, defects


def load_context(initial_result, coefficient_manifest, channel_archive, contract_path):
    contract_binding = binding(contract_path)
    compatibility = load_compatibility()
    if not same_binding(contract_binding, compatibility['implementation_contract']):
        raise ValueError('scientific profile binds a different original A34 plan')
    contract = load_json(contract_path)
    if contract['schema'] != 'P254-A34-full-source-two-mode-initial-force-noise-plan-v1':
        raise ValueError('A34 implementation contract001 required')
    inputs = dict(initial_result=binding(initial_result),
                  coefficient_manifest=binding(coefficient_manifest),
                  channel_archive=binding(channel_archive))
    for key, authority in (('initial_result', 'original_initial'),
                           ('coefficient_manifest', 'original_coefficients'),
                           ('channel_archive', 'mode_archive')):
        if not same_binding(inputs[key], contract[authority]):
            raise ValueError('input differs from original A34 authority: ' + key)
    initial = load_json(initial_result)
    if initial['status'] != 'complete' or initial['mode'] != 'initialize':
        raise ValueError('complete original initialization required')
    saved = initial['saved']
    original_binding = authenticate(saved['original_boundary'])
    preparation = authenticate(saved['preparation'])
    boundary = load_json(original_binding['path'])
    slots = boundary['slots']
    if (boundary['schema'] != 'streamed-material-mean-v1' or len(slots) != 128
            or boundary['time_ry_inverse'] != 0
            or any(boundary['u_bohr']) or any(boundary['p_canonical'])
            or boundary['model_identity'] != saved['original_model_identity']
            or len({slot['identity'] for slot in slots}) != 128):
        raise ValueError('literal original complete128 t=0 boundary required')
    manifest = load_json(coefficient_manifest)
    if (manifest['schema'] != 'p254-physical-activeband-cache-v1'
            or manifest['status'] != 'complete' or manifest['specification']['bands'] != 12
            or not same_binding(manifest['specification']['original_boundary'], original_binding)):
        raise ValueError('matching original saved12 coefficient manifest required')
    selected = {}
    for key, item in manifest['items'].items():
        if not key.startswith('slot/'):
            continue
        identity = item['metadata']['identity']
        if identity in selected:
            raise ValueError('duplicate coefficient sector identity')
        selected[identity] = (key, item)
    if set(selected) != {slot['identity'] for slot in slots}:
        raise ValueError('coefficient manifest does not cover original128 exactly')
    for slot in slots:
        if number(slot['weight'], 'original quadrature') != 1 / 64:
            raise ValueError('expected actual original 1/64 sector quadrature')
        _, item = selected[slot['identity']]
        metadata = item['metadata']
        if (metadata['basis_identity'] != slot['basis_identity']
                or metadata['frame_identity'] != slot['gauge_identity']
                or not same_binding(metadata['frame'], slot['frame'])
                or not same_binding(metadata['raw_d'], slot['coefficient'])):
            raise ValueError('saved coefficients differ from original frame/D')
    return dict(contract=contract, contract_binding=contract_binding, inputs=inputs,
                initial=initial, boundary=boundary, slots=slots, selected=selected,
                original_boundary=original_binding, preparation=preparation,
                scientific_compatibility=compatibility,
                scientific_compatibility_identity=compatibility['scientific_compatibility_identity'])


def load_case(done_path, context, expected_record=None):
    """The completion receipt, not a partial result filename, is authoritative."""
    done_path = Path(done_path).resolve(strict=True)
    if done_path.name != 'done.json':
        raise ValueError('case completion must be named done.json')
    done_binding = binding(done_path)
    done = load_json(done_path)
    result_binding = authenticate(done['result'], done_path.parent / 'result.json')
    result = load_json(result_binding['path'])
    if result['schema'] != RESULT_SCHEMA or result['status'] != 'complete':
        raise ValueError('complete A34 original-sector result required')
    record = result['original_record']
    if type(record) is not int or record not in range(128):
        raise ValueError('exact original record in 0..127 required')
    if expected_record is not None and record != expected_record:
        raise ValueError('completion belongs to another original record')
    scientific_admission = admit_result(result, result_binding,
                                        context['scientific_compatibility'])
    if (result['quadrature_applied'] is not False or result['mode_indices'] != [0, 1]
            or any(type(index) is not int for index in result['mode_indices'])):
        raise ValueError('unweighted exact modes [0,1] required')
    if not same_binding(result['mode_archive'], context['inputs']['channel_archive']):
        raise ValueError('different paid output-directed modes')
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
    frame = result.get('frame', source.get('frame'))
    if frame is None or not same_binding(frame, slot['frame']):
        raise ValueError('original full frame binding required')
    expected_metadata = {key: slot[key] for key in (
        'identity', 'component', 'weight', 'shape', 'basis_identity',
        'gauge_identity', 'source_identity')}
    if source['original_metadata'] != expected_metadata:
        raise ValueError('source metadata differs from literal original slot')
    key, item = context['selected'][slot['identity']]
    if (source['coefficient_item'] != key
            or not same_binding(source['coefficient_archive'], item['archive'])
            or source['selected'] != item['metadata']['selected']):
        raise ValueError('source active selection differs from original coefficients')
    if (source['case'] != context['initial']['case']
            or source['model_identity'] != context['initial']['saved']['model_identity']
            or source['original_model_identity'] != context['boundary']['model_identity']):
        raise ValueError('source model/preparation mismatch')
    if type(source['spin']) is not int or source['spin'] not in (1, 2):
        raise ValueError('literal explicit spin 1/2 required')
    mesh = source['mesh_index']
    if not isinstance(mesh, list) or len(mesh) != 3 or any(type(i) is not int for i in mesh):
        raise ValueError('exact three-index internal mesh required')
    explicit_identity = str(source['spin']) + ':' + ','.join(str(i) for i in mesh)
    if (explicit_identity != slot['identity']
            or source['spin'] - 1 != slot['component']):
        raise ValueError('spin/internal mesh differs from original record identity')
    work = result['work']
    for key in ('elapsed_seconds', 'peak_rss_bytes', 'physical_force_forward_actions',
                'physical_force_adjoint_actions', 'hamiltonian_actions', 'trajectory_actions'):
        number(work[key], key, nonnegative=True)
    expected_actions = 0 if record == 0 else 2
    if (work['physical_force_forward_actions'] != expected_actions
            or work['physical_force_adjoint_actions'] != expected_actions
            or work['hamiltonian_actions'] != 0 or work['trajectory_actions'] != 0):
        raise ValueError('case violates paid-source0/two-pair/no-H/no-trajectory work contract')
    parsed = parse_covariance(result['covariance'])
    return dict(original_record=record, done=done_binding, result=result_binding,
                payload=result, parsed=parsed, scientific_admission=scientific_admission)


def aggregate(cases, context):
    by_record, spin_mesh = {}, set()
    for case in cases:
        record = case['original_record']
        if record in by_record:
            raise ValueError('duplicate original record: ' + str(record))
        source = case['payload']['source']
        identity = (source['spin'], tuple(source['mesh_index']))
        if identity in spin_mesh:
            raise ValueError('duplicate explicit spin/internal-k sector')
        spin_mesh.add(identity)
        by_record[record] = case
    present = sorted(by_record)
    missing = sorted(set(range(128)) - set(present))
    if not missing:
        meshes = [{mesh for spin, mesh in spin_mesh if spin == s} for s in (1, 2)]
        if len(meshes[0]) != 64 or meshes[0] != meshes[1]:
            raise ValueError('full coverage must contain both explicit spins on the same64 mesh')
    rows = []
    for a in MODES:
        for b in MODES:
            pair = (a, b)
            full = complex_sum(case['payload']['quadrature_weight'] * case['parsed'][0][pair]
                               for case in by_record.values())
            blocks = []
            for pa in PARTS:
                for pb in PARTS:
                    value = complex_sum(case['payload']['quadrature_weight']
                                        * case['parsed'][1][pair + (pa, pb)]
                                        for case in by_record.values())
                    blocks.append(dict(part_left=pa, part_right=pb, value=complex_record(value)))
            reconstructed = complex_sum(complex_value(block['value']) for block in blocks)
            raw_defect = complex_sum(case['payload']['quadrature_weight'] * case['parsed'][2][pair]
                                     for case in by_record.values())
            rows.append(dict(left=a, right=b, full=complex_record(full), blocks=blocks,
                             block_reconstruction=complex_record(reconstructed),
                             reconstruction_defect=complex_record(full - reconstructed),
                             weighted_raw_reconstruction_defect=complex_record(raw_defect)))
    return dict(schema='P254-A34-initial-Gamma-two-mode-force-noise-aggregate-v1',
                status='complete' if not missing else 'partial',
                original_records=present, missing_original_records=missing,
                required_original_records=list(range(128)), sector_count=len(present),
                contract=context['contract_binding'], inputs=context['inputs'],
                scientific_compatibility_identity=context['scientific_compatibility_identity'],
                original_boundary=context['original_boundary'], preparation=context['preparation'],
                mode_indices=[0, 1], mode_archive=context['inputs']['channel_archive'],
                normalization=dict(rule='sum_k actual_quadrature_weight_k * unweighted_C_k',
                                   quadrature_applied_once=True, quadrature_weight=1 / 64,
                                   explicit_spin_sectors=128, extra_spin_factor=1,
                                   final_spatial_factor=1, partial_renormalization=False),
                covariance=rows,
                sectors=[dict(original_record=r, done=by_record[r]['done'],
                              result=by_record[r]['result'],
                              scientific_admission=by_record[r]['scientific_admission'],
                              quadrature_weight=by_record[r]['payload']['quadrature_weight'],
                              source=by_record[r]['payload']['source'],
                              work=by_record[r]['payload']['work'],
                              raw_unweighted_covariance=by_record[r]['payload']['covariance'])
                         for r in present],
                numerical_repairs_applied=False, scope=SCOPE)


def add_input_arguments(parser):
    parser.add_argument('--initial-result', type=Path, required=True)
    parser.add_argument('--coefficient-manifest', type=Path, required=True)
    parser.add_argument('--channel-archive', type=Path, required=True)
    parser.add_argument('--contract', type=Path,
                        default=Path(__file__).with_name('implementation-contract-001.json'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_input_arguments(parser)
    parser.add_argument('--done', type=Path, action='append', default=[],
                        help='explicit committed case done.json; repeat for each sector')
    parser.add_argument('--batch-summary', type=Path, action='append', default=[],
                        help='batch run summary; uses only its explicit completed receipt list')
    parser.add_argument('--output', type=Path, required=True,
                        help='new exclusive aggregate JSON path (parent directory must exist)')
    args = parser.parse_args()
    try:
        context = load_context(args.initial_result, args.coefficient_manifest,
                               args.channel_archive, args.contract)
        paths = list(args.done)
        for path in args.batch_summary:
            summary = load_json(path)
            if summary['schema'] != 'P254-A34-mode-noise-batch-summary-v1':
                raise ValueError('wrong batch summary schema')
            for item in summary['completed']:
                authenticate(item['done'])
                paths.append(Path(item['done']['path']))
        if not paths:
            raise ValueError('supply --done or --batch-summary committed cases')
        result = aggregate([load_case(path, context) for path in paths], context)
        code = 0 if result['status'] == 'complete' else 2
    except Exception as error:
        result = dict(schema='P254-A34-initial-Gamma-two-mode-force-noise-aggregate-v1',
                      status='invalid', error=repr(error), scope=SCOPE)
        code = 1
    receipt = save_json(args.output, result)
    print(json.dumps(dict(event='aggregate-exit', status=result['status'],
                          output=receipt, exit_code=code), allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
