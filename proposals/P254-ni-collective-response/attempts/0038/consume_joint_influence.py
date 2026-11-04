"""Consume A38 raw current matrices and paid A37 ionic/spectral banks only.

This program restores no material context and performs no physical action,
original material preparation, H/basis/spectral eigensolve, inverse or original
ionic integration. The unchanged CAR state-admissibility diagnostics are charged.
New field/joint operators are saved before the original-state ordered CAR readout.
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
    for name in ('0037', '0036', '0031', '0030', '0029', '0028')]

import numpy as np
from canonical_current_noise import current_force_matrix
from expanded_influence_readout import save_json, save_npz
from integrate_saved_spectral import load_saved_spectral
from joint_influence_readout import (PARTS, SCOPE, field_source_amplitude,
    integrate_field, joint_moments, joint_operators, joint_reconstruction)
from substrate_framework.streamed_material_state import authenticate, file_binding

CONTRACT = dict(path=str(HERE / 'implementation-contract-001.json'), bytes=15306,
    sha256='623ef2c7ebc3fec3b0166783b165458e349abaf6d4d06c4da8bb5c77433f2668')
CLARIFICATION = dict(path=str(HERE / 'unchanged-loader-and-CAR-cost-clarification-001.json'),
    bytes=1010, sha256='a0ab6089219256eb86d3d89ecd8b2dc661d25b1ce2920d96c4ae6e2ce76a49ba')
CURRENT_NAMES = ('current_x', 'current_y', 'current_z')


def read_bound(binding):
    return json.loads(authenticate(binding).read_bytes())


def same_binding(left, right):
    return (Path(left['path']).resolve() == Path(right['path']).resolve()
        and int(left['bytes']) == int(right['bytes']) and left['sha256'] == right['sha256'])


def require_binding(actual, expected, name):
    if not same_binding(actual, expected):
        raise ValueError('settled binding differs: ' + name)


def authenticate_handoff(current_projection_complete):
    """Authenticate original complete receipts and new producer, without restoration."""
    contract = read_bound(CONTRACT)
    clarification = read_bound(CLARIFICATION)
    require_binding(clarification['original_contract'], CONTRACT, 'CAR cost clarification')
    settled = {Path(row['path']).resolve(): row for row in contract['settled_source_inventory']}
    for binding in settled.values():
        authenticate(binding)
    by_name = {path.name: row for path, row in settled.items()}
    producer_binding = file_binding(current_projection_complete)
    producer = read_bound(producer_binding)
    if producer['schema'] != 'P254-A38-current-projection-v1' or producer['status'] != 'complete':
        raise ValueError('complete canonical A38 current projection required')
    for name in contract['handoff']['required_manifest_fields']:
        if name not in producer:
            raise ValueError('incomplete current handoff: ' + name)
    expected = dict(basis_complete='basis-complete.json', basis='finite-bath-basis.npz',
        a37_filtered_complete='filtered-complete.json', a37_filtered_operators='filtered-operators.npz')
    # result.json occurs for A35 and A37: bind the full settled path, not a basename.
    expected['a37_temporal_result'] = '/var/tmp/p254-a37-source0-two-mode-influence692-001/result.json'
    for field, name in expected.items():
        target = settled[Path(name)] if name.startswith('/') else by_name[name]
        require_binding(producer[field], target, field)
    require_binding(producer['current_reference'],
        settled[Path('/var/tmp/p254-a35-uniform-current-source0-002/result.json')], 'current_reference')
    for binding in producer['implementation']:
        authenticate(binding)
    basis = read_bound(producer['basis_complete'])
    complete = read_bound(producer['a37_filtered_complete'])
    temporal = read_bound(producer['a37_temporal_result'])
    if any(item['status'] != 'complete' for item in (basis, complete, temporal)):
        raise ValueError('all original A37 receipts must be complete')
    for item in (basis, complete, temporal):
        require_binding(item['basis'], producer['basis'], 'original basis')
    require_binding(complete['filtered_operators'], producer['a37_filtered_operators'], 'filtered receipt')
    require_binding(temporal['filtered_operators'], producer['a37_filtered_operators'], 'temporal receipt')
    for field in ('original_state_authority', 'original_factor_identity', 'original_pw_basis_identity'):
        if producer[field] != basis[field] or temporal[field] != basis[field]:
            raise ValueError('unchanged original authority differs: ' + field)
    if producer['original_source'] != basis['original_source'] or temporal['source'] != basis['original_source']:
        raise ValueError('original computation provenance differs')
    if (basis['basis_metadata']['dimension'] != 1038
            or basis['basis_metadata']['original_rank'] != 346
            or basis['basis_metadata']['extra_rank_selected'] != 692):
        raise ValueError('all346 original columns and unchanged692 extras required')
    if (producer['source_record'] != 0 or producer['ncells'] != 64
            or producer['quadrature_weight'] != 1 / 64 or basis['ncells'] != 64):
        raise ValueError('source0/N64/w1over64 canonical normalization required')
    field_result = read_bound(by_name['canonical-field-force-source0-partial-material-result-002.json'])
    polarizations = field_result['field_normalization']['polarizations']
    matrix = current_force_matrix(polarizations)
    amplitude = field_source_amplitude(64, 1 / 64, matrix)
    if (producer['polarizations'] != polarizations
            or not np.array_equal(matrix, field_result['normalization']['current_force_matrix'])
            or not np.array_equal(matrix, producer['current_force_matrix'])
            or not np.array_equal(amplitude, producer['source_field_amplitude'])):
        raise ValueError('original A36 polarizations/map/source amplitude differs')
    filter_binding = by_name['packet.json']
    require_binding(basis['filter_packet'], filter_binding, 'original A30 filter packet')
    packet = read_bound(filter_binding)
    require_binding(packet['archive'], by_name['canonical-force-filter.npz'], 'original filter archive')
    require_binding(basis['filter_archive'], packet['archive'], 'basis filter archive')
    observable = packet['filter']
    if (observable['coordinates'] != [0] or observable['q_index'] != 0
            or observable['nprobes'] != 1 or observable['accepted_intervals'] != 77
            or observable['polynomial_degree'] != 7 or observable['final_time_inverse_ry'] != 30.0):
        raise ValueError('original Ni1/x Gamma degree7/all77 filter required')
    aliases = producer['current_adjoint_aliases']
    if aliases != {name: name for name in CURRENT_NAMES}:
        raise ValueError('qualified true-current adjoint aliases required')
    return producer, basis, temporal, amplitude, dict(implementation_contract=dict(CONTRACT),
        unchanged_CAR_cost_clarification=dict(CLARIFICATION),
        current_projection_complete=producer_binding, a37_basis_complete=producer['basis_complete'],
        a37_temporal_result=producer['a37_temporal_result'], a37_filtered_complete=producer['a37_filtered_complete'],
        a37_filtered_operators=producer['a37_filtered_operators'], original_filter_packet=filter_binding,
        original_filter_archive=packet['archive'], a36_canonical_force_result=
            by_name['canonical-field-force-source0-partial-material-result-002.json'])


def load_numeric(producer, temporal, bindings):
    """Read only required original state/projector/filter and paid ionic operators."""
    with np.load(authenticate(producer['basis']), allow_pickle=False) as archive:
        bank = {key: archive[key] for key in ('projected_frame', 'gram',
            'density_coefficients', 'active_projector', 'segment_bounds')}
    shapes = dict(projected_frame=(1038, 346), gram=(346, 346),
        density_coefficients=(346, 346), active_projector=(1038, 1038), segment_bounds=(77, 2))
    for key, shape in shapes.items():
        if bank[key].shape != shape or not np.all(np.isfinite(bank[key])):
            raise ValueError('original state/projector/filter shape or finiteness differs: ' + key)
    if not np.array_equal(bank['segment_bounds'], temporal['integration_metadata']['intervals']):
        raise ValueError('actual A37 intervals differ from original temporal receipt')
    with np.load(authenticate(bindings['original_filter_archive']), allow_pickle=False) as archive:
        field = archive['field_chebyshev']
        if not np.array_equal(bank['segment_bounds'], archive['segment_bounds']):
            raise ValueError('actual A30/A37 accepted intervals differ')
    if field.shape != (77, 8, 3, 1) or not np.all(np.isfinite(field)):
        raise ValueError('original complete canonical field filter required')
    with np.load(authenticate(producer['current_bank']), allow_pickle=False) as archive:
        if set(archive.files) != set(CURRENT_NAMES):
            raise ValueError('current bank requires exactly current_x/y/z')
        currents = [archive[name] for name in CURRENT_NAMES]
    if any(value.shape != (1038, 1038) or value.dtype != np.complex128
           or not np.all(np.isfinite(value)) for value in currents):
        raise ValueError('raw current matrices must be finite complex128[1038,1038]')
    with np.load(authenticate(producer['a37_filtered_operators']), allow_pickle=False) as archive:
        ionic = {name: archive['operator_' + name] for name in PARTS}
        adjoints = {name: archive['adjoint_operator_' + name] for name in PARTS}
    if any(value.shape != (1038, 1038) or value.dtype != np.complex128
           or not np.all(np.isfinite(value)) for value in (*ionic.values(), *adjoints.values())):
        raise ValueError('paid ionic operators require finite complex128[1038,1038]')
    saved = load_saved_spectral(producer['a37_filtered_operators'])
    if saved['eigenvalues'].size != 1038:
        raise ValueError('actual A37 spectral dimension differs')
    return bank, field, currents, ionic, adjoints, saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--current-projection-complete', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = perf_counter()
    producer, basis, temporal, amplitude, bindings = authenticate_handoff(args.current_projection_complete)
    bank, field, currents, ionic, ionic_adjoints, saved = load_numeric(producer, temporal, bindings)
    loading_seconds = perf_counter() - started
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    integration_started = perf_counter()
    # True adjoints alias the qualified same raw J actions, not matrix.conj().T.
    filtered = integrate_field(saved, currents, currents, bank['active_projector'],
        bank['segment_bounds'], field, amplitude)
    integration_seconds = perf_counter() - integration_started
    matrices, adjoints, arrays = joint_operators(ionic, ionic_adjoints, filtered)
    operator_binding = save_npz(directory / 'joint-operators.npz', arrays)
    complete_binding = save_json(directory / 'joint-operators-complete.json', dict(
        schema='P254-A38-joint-operators-v1', status='complete', **bindings,
        basis=producer['basis'], current_bank=producer['current_bank'], operators=operator_binding,
        spectral_reference=producer['a37_filtered_operators'], spectral_members='spectral_*; not duplicated',
        integration_metadata=filtered['metadata'], integration_diagnostics=filtered['diagnostics'], scope=SCOPE))
    contraction_started = perf_counter()
    moments = joint_moments(bank, matrices, adjoints, basis['basis_identity'],
        producer['original_factor_identity'], producer['original_state_authority'])
    reconstruction = joint_reconstruction(moments, matrices, adjoints)
    contraction_seconds = perf_counter() - contraction_started
    implementation = [file_binding(path) for path in (Path(__file__),
        HERE / 'joint_influence_readout.py', HERE / 'integrate_saved_spectral.py',
        HERE.parent / '0037/expanded_influence_readout.py',
        HERE.parent / '0036/canonical_current_noise.py',
        HERE.parent / '0031/filtered_quadratic_bath.py',
        HERE.parent / '0031/export_frozen_bath_comparator.py',
        ROOT / 'src/substrate_framework/physical_factor_preparation.py',
        ROOT / 'src/substrate_framework/streamed_material_state.py')]
    result = dict(schema='P254-A38-joint-frozen-H-temporal-influence-v1', status='complete',
        **bindings, basis=producer['basis'], current_bank=producer['current_bank'],
        joint_operators=operator_binding, joint_operators_complete=complete_binding,
        spectral_reference=producer['a37_filtered_operators'], source=producer['original_source'],
        original_state_authority=producer['original_state_authority'],
        original_factor_identity=producer['original_factor_identity'],
        original_pw_basis_identity=producer['original_pw_basis_identity'], basis_metadata=basis['basis_metadata'],
        observable=dict(output='Ni1/x', q='Gamma', original_probe_index=0, final_time_inverse_ry=30.0,
            mean_units='bohr', ordered_covariance_units='bohr^2',
            covariance_convention='connected <delta A delta B>; no conjugation of B'),
        filtered_covariance=moments, **reconstruction,
        normalization=dict(ncells=64, source_record=0, quadrature_weight=1 / 64,
            polarizations=producer['polarizations'], current_force_matrix=producer['current_force_matrix'],
            source_field_amplitude=amplitude, ionic_policy='reuse paid A37 F/sqrt(N) operators unchanged',
            field_policy='field_chebyshev times sqrt(N*w0)*M times raw J; M=2e; N*w0=1',
            extra_spin_factor=1, extra_spatial_factor=1, mass_division=False,
            inverse_Fourier_factor=False, independent_noise_addition=False),
        initial_covariance=producer['initial_covariance'],
        integration_metadata=filtered['metadata'], integration_diagnostics=filtered['diagnostics'],
        carried_owning_proofs=[row for row in json.loads(authenticate(CONTRACT).read_bytes())['settled_source_inventory']
            if Path(row['path']).name in ('owning-block-oracles-001.json', 'owning-canonical-field-force-oracle-001.json')],
        work=dict(loading_seconds=loading_seconds, field_integration_seconds=integration_seconds,
            contraction_seconds=contraction_seconds, total_seconds=perf_counter() - started,
            peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            consumer_physical_actions=0, consumer_material_preparations=0, consumer_H_F_replays=0,
            consumer_H_basis_spectral_eigensolves=0, producer_work=producer['work'],
            unchanged_CAR_state_admissibility='existing small Gram/state eigendiagnostics charged in contraction_seconds',
            new_operator_payload_bytes=sum(value.nbytes for value in arrays.values()),
            spectral_bank_copied=False, original_A37_wholecost_seconds_carried_closed=1840.656),
        implementation=implementation, producer_implementation=producer['implementation'], scope=SCOPE)
    binding = save_json(directory / 'result.json', result)
    print(json.dumps(dict(result=binding, joint_operators_complete=complete_binding,
        total_seconds=result['work']['total_seconds']), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
