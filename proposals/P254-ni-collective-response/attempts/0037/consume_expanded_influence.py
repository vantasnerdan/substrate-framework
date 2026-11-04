"""Consume the complete A37 basis bank: no physical context or action replay.

The unmodified A31 integrator and CAR consumer compute one source0/Gamma/Ni1-x
nonstationary, frozen-H contribution. The old8/new692 contrast is not an error
bound or a correction to the moving A32 total.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / 'src')] + [str(HERE.parent / name)
    for name in ('0031', '0030', '0029', '0028')]

import numpy as np
from export_frozen_bath_comparator import car_moments, projected_factors
from filtered_quadratic_bath import integrate_filtered_operators
from substrate_framework.streamed_material_state import authenticate, file_binding
from expanded_influence_readout import (SCOPE, VERTICES, model_contrast,
    raw_bank_defects, reconstruction, save_json, save_npz)

CONTRACT = HERE / 'implementation-contract-001.json'
CONTRACT_SHA = '037bfd5014b6f9ce994a0658af5f0e8d4b36031fd5c3fb92238ced5a197a61f1'
NORMALIZATION = dict(mesh=[4, 4, 4], ncells=64,
    force_insertion='unweighted F_k/sqrt(N)', further_source_weight_applied=False,
    spin_factor_applied=False, final_real_space_inverse_fourier_applied=False)


def same_binding(left, right):
    return (left['sha256'] == right['sha256']
            and int(left['bytes']) == int(right['bytes']))


def read_bound(binding):
    return json.loads(authenticate(binding).read_bytes())


def authenticate_handoff(basis_complete, reference_result):
    """Bind the complete bank and its original observable, without restoring it."""
    contract_binding = file_binding(CONTRACT)
    if contract_binding['sha256'] != CONTRACT_SHA or contract_binding['bytes'] != 9415:
        raise ValueError('the settled A37 implementation contract is required')
    contract = json.loads(CONTRACT.read_bytes())
    settled = {Path(row['path']).name: row for row in contract['settled_inputs']}
    # These are the unchanged numerical functions actually imported above.
    for name in ('cached_bath_inputs.py', 'noise_subspace.py',
                 'filtered_quadratic_bath.py', 'export_frozen_bath_comparator.py'):
        authenticate(settled[name])
    manifest_binding = file_binding(basis_complete)
    manifest = read_bound(manifest_binding)
    if manifest['status'] != contract['handoff']['status']:
        raise ValueError('basis-complete must have complete status')
    for field in contract['handoff']['metadata_fields']:
        if field not in manifest:
            raise ValueError('missing handoff field: ' + field)
    metadata = manifest['basis_metadata']
    if (metadata['original_rank'] != 346 or metadata['extra_rank_requested'] != 692
            or metadata['extra_rank_selected'] != 692 or metadata['dimension'] != 1038):
        raise ValueError('require all346 original raw columns and exactly692 extra columns')
    reference_binding = file_binding(reference_result)
    if not same_binding(reference_binding, settled['finite-bath-material-result-001.json']):
        raise ValueError('require the actual complete oldA31 material result')
    reference = read_bound(reference_binding)
    _match_reference(manifest, reference)
    expected_paths = {
        'action_packet': '/var/tmp/p254-a30-output-directed-source0-actions-001/packet.json',
        'mode1_packet': '/var/tmp/p254-a31-source0-two-mode-bath8-001/mode1-complete.json',
        'filter_packet': '/var/tmp/p254-a30-loaded-gamma-ni1-causal-filter-001/packet.json',
        'filter_archive': '/var/tmp/p254-a30-loaded-gamma-ni1-causal-filter-001/canonical-force-filter.npz'}
    by_path = {row['path']: row for row in contract['settled_inputs']}
    for key, path in expected_paths.items():
        if not same_binding(manifest[key], by_path[path]):
            raise ValueError('settled original binding differs: ' + key)
        authenticate(manifest[key])
    filter_packet = read_bound(manifest['filter_packet'])
    if not same_binding(filter_packet['archive'], manifest['filter_archive']):
        raise ValueError('canonical filter archive differs from its packet')
    source_binding = filter_packet['source']['original_result']
    if not same_binding(source_binding, settled['result.json']):
        raise ValueError('original normalization source differs')
    original = read_bound(source_binding)
    mesh = original['configuration']['coefficient_specification']['mesh']
    if (mesh != manifest['mesh'] or mesh != NORMALIZATION['mesh']
            or manifest['ncells'] != 64 or int(np.prod(mesh)) != 64):
        raise ValueError('original actual filter mesh/N64 differs')
    observable = filter_packet['filter']
    if (observable['coordinates'] != [0] or observable['q_index'] != 0
            or observable['final_time_inverse_ry'] != 30.0 or observable['nprobes'] != 1):
        raise ValueError('require the original first Ni1/x Gamma output filter')
    for binding in manifest['implementation']:
        authenticate(binding)
    bank_path = authenticate(manifest['basis'])
    with np.load(bank_path, allow_pickle=False) as archive:
        if set(archive.files) != set(contract['handoff']['npz_keys']):
            raise ValueError('saved bank keys differ from the complete handoff')
        bank = {key: archive[key] for key in archive.files}
    _check_bank(bank, metadata, reference)
    return manifest, bank, reference, dict(basis_complete=manifest_binding,
        reference_result=reference_binding, implementation_contract=contract_binding,
        normalization_source=source_binding, original_filter_observable=observable)


def _match_reference(manifest, reference):
    if reference['status'] != 'complete' or reference['basis_metadata']['extra_rank_selected'] != 8:
        raise ValueError('reference must be the complete old8 model')
    if reference['normalization'] != NORMALIZATION:
        raise ValueError('oldA31 observable normalization differs')
    for field in ('original_state_authority', 'original_factor_identity', 'original_pw_basis_identity'):
        if manifest[field] != reference[field]:
            raise ValueError('original state/factor/basis differs: ' + field)
    if manifest['original_source'] != reference['source']:
        raise ValueError('source0 original provenance differs from oldA31')
    if manifest['basis_metadata']['operators'] != reference['basis_metadata']['operators']:
        raise ValueError('literal H/two-mode operator identities differ from oldA31')
    for field, old in (('action_packet', 'original_action_packet'),
                       ('filter_packet', 'filter_packet'), ('filter_archive', 'filter_archive')):
        if not same_binding(manifest[field], reference[old]):
            raise ValueError('same-observable reference binding differs: ' + field)


def _check_bank(bank, metadata, reference):
    shapes = dict(extra_basis=(metadata['full_pw_dimension'], 692),
        inverse_sqrt_gram=(346, 346), projected_frame=(1038, 346), gram=(346, 346),
        density_coefficients=(346, 346), active_coordinates=(346, metadata['active_count']),
        projected_active_coordinates=(1038, metadata['active_count']))
    for key in ('hamiltonian', 'hamiltonian_adjoint', 'projected_density',
                'active_projector', 'basis_metric', 'force_mode0', 'force_mode1',
                'adjoint_force_mode0', 'adjoint_force_mode1'):
        shapes[key] = (1038, 1038)
    for key, shape in shapes.items():
        if bank[key].shape != shape:
            raise ValueError('complete raw bank shape differs: ' + key)
    coefficients = bank['modal_coefficients']
    if coefficients.ndim != 4 or coefficients.shape[1:] != (8, 2, 1):
        raise ValueError('require original coefficients[segment,degree,2,first Ni1/x probe]')
    intervals = np.asarray(reference['integration_metadata']['intervals'])
    if (bank['segment_bounds'].shape != intervals.shape
            or not np.array_equal(bank['segment_bounds'], intervals)
            or coefficients.shape[0] != len(intervals)):
        raise ValueError('original accepted filter intervals differ from oldA31')


def integrate_bank(bank, ncells):
    """The exact old comparator: four force channels and two output probes."""
    projector = bank['active_projector']
    force = np.stack([bank['force_mode0'], bank['force_mode1']])
    adjoint = np.stack([bank['adjoint_force_mode0'], bank['adjoint_force_mode1']])
    retained = np.stack([projector @ matrix @ projector for matrix in force])
    retained_adjoint = np.stack([projector @ matrix @ projector for matrix in adjoint])
    coefficients = bank['modal_coefficients']
    batched = np.zeros(coefficients.shape[:2] + (4, 2), dtype=np.complex128)
    batched[:, :, :2, 0] = coefficients[:, :, :, 0]
    batched[:, :, 2:, 1] = coefficients[:, :, :, 0]
    # Insert F/sqrt(N64) ONCE; no source, spin or final inverse-Fourier weight.
    return integrate_filtered_operators(bank['hamiltonian'],
        np.concatenate((force, retained)) / np.sqrt(ncells),
        np.concatenate((adjoint, retained_adjoint)) / np.sqrt(ncells),
        bank['segment_bounds'], batched)


def save_filtered(directory, filtered, input_bindings, basis_binding):
    matrices = dict(zip(VERTICES, (filtered['operators'][0], filtered['operators'][1],
        filtered['operators'][0] - filtered['operators'][1])))
    adjoints = dict(zip(VERTICES, (filtered['adjoint_operators'][0],
        filtered['adjoint_operators'][1],
        filtered['adjoint_operators'][0] - filtered['adjoint_operators'][1])))
    arrays = {'operator_' + key: value for key, value in matrices.items()}
    arrays.update({'adjoint_operator_' + key: value for key, value in adjoints.items()})
    # Actual integrator coordinates are saved, never regenerated from operators.
    arrays.update({'spectral_' + key: value for key, value in filtered['spectral'].items()})
    binding = save_npz(directory / 'filtered-operators.npz', arrays)
    complete = save_json(directory / 'filtered-complete.json', dict(status='complete',
        basis=basis_binding, inputs=input_bindings, filtered_operators=binding,
        integration_metadata=filtered['metadata'], integration_diagnostics=filtered['diagnostics'],
        normalization=NORMALIZATION, scope=SCOPE))
    return matrices, adjoints, binding, complete


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basis-complete', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference-result', required=True)
    args = parser.parse_args()
    started = perf_counter()
    manifest, bank, reference, bindings = authenticate_handoff(
        args.basis_complete, args.reference_result)
    loading_seconds = perf_counter() - started
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    print('[physics] authenticated expanded1038 basis; integrating all77 original '
          'degree-seven segments, four channels/two probes', flush=True)
    integration_started = perf_counter()
    filtered = integrate_bank(bank, manifest['ncells'])
    integration_seconds = perf_counter() - integration_started
    matrices, adjoints, operator_binding, filtered_complete = save_filtered(
        directory, filtered, bindings, manifest['basis'])
    print('[physics] actual spectral/eigencoordinate bank durably saved; '
          'contracting all9 ordered covariances', flush=True)
    basis_identity = 'A37 finite-bath coordinate basis ' + manifest['basis']['sha256']
    contraction_started = perf_counter()
    moments = car_moments(bank['projected_frame'], bank['gram'], bank['density_coefficients'],
        projected_factors(matrices, adjoints, bank['projected_frame'], basis_identity,
            manifest['original_state_authority'], 'bohr'), manifest['original_factor_identity'],
        basis_identity, manifest['original_state_authority'])
    contraction_seconds = perf_counter() - contraction_started
    result = dict(schema='P254-A37-expanded-frozen-H-temporal-influence-v1', status='complete',
        **bindings, basis=manifest['basis'], filtered_operators=operator_binding,
        filtered_complete=filtered_complete, source=manifest['original_source'],
        original_state_authority=manifest['original_state_authority'],
        original_factor_identity=manifest['original_factor_identity'],
        original_pw_basis_identity=manifest['original_pw_basis_identity'],
        original_action_packet=manifest['action_packet'], mode1_packet=manifest['mode1_packet'],
        filter_packet=manifest['filter_packet'], filter_archive=manifest['filter_archive'],
        basis_metadata=manifest['basis_metadata'], basis_diagnostics=manifest['basis_diagnostics'],
        initial_unweighted_force_covariance=manifest['initial_force_covariance'],
        observable=dict(output='Ni1/x', q='Gamma', original_probe_index=0,
            final_time_inverse_ry=30.0, mean_units='bohr', ordered_covariance_units='bohr^2',
            covariance_convention='connected <delta A delta B>; no conjugation of B'),
        filtered_covariance=moments, **reconstruction(moments, matrices, adjoints),
        **raw_bank_defects(bank), normalization=NORMALIZATION,
        old8_new692_model_enrichment=model_contrast(reference, moments, bindings['reference_result']),
        integration_metadata=filtered['metadata'], integration_diagnostics=filtered['diagnostics'],
        work=dict(loading_seconds=loading_seconds, integration_seconds=integration_seconds,
            contraction_seconds=contraction_seconds, total_seconds=perf_counter() - started,
            peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            saved_bank_array_bytes=sum(value.nbytes for value in bank.values()),
            consumer_new_physical_actions=0, consumer_original_material_preparations=0,
            consumer_original_action_replays=0, producer_work=manifest['work'],
            producer_new_physical_actions=manifest['new_physical_actions']),
        implementation=[file_binding(path) for path in (
            Path(__file__), HERE / 'expanded_influence_readout.py',
            HERE.parent / '0031' / 'export_frozen_bath_comparator.py',
            HERE.parent / '0031' / 'filtered_quadratic_bath.py',
            ROOT / 'src/substrate_framework/physical_factor_preparation.py',
            ROOT / 'src/substrate_framework/streamed_material_state.py',
            HERE.parent / '0028' / 'clamped_scientific_readout.py')],
        producer_implementation=manifest['implementation'], scope=SCOPE)
    binding = save_json(directory / 'result.json', result)
    print(json.dumps(dict(result=binding, filtered_complete=filtered_complete,
        total_seconds=result['work']['total_seconds']), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
