"""Consume paid moving readouts and finite-bath data; never evolve a trajectory."""
import argparse
import json
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src')]+[str(HERE.parent/name)
    for name in ('0032', '0031', '0030', '0029', '0028')]

import numpy as np
from export_hybrid_gamma import _arrays, _json, publish_json
from substrate_framework.streamed_material_state import authenticate, file_binding
from reduced_coherence_readout import normalized_ordered
from hybrid_noise_components import decompose_initial_noise
from finite_bath_influence import frozen_force_influence


def common_normalized(matrix, full_variances):
    values = np.empty_like(matrix, dtype=complex)
    for left in range(len(full_variances)):
        for right in range(len(full_variances)):
            row = normalized_ordered(matrix[left, right], full_variances[left],
                                     full_variances[right])
            if not row['defined']:
                raise ValueError('the actual full ordered normalization is undefined')
            values[left, right] = row['complex_degree']
    return values


def run(moving_result_path, output):
    started = perf_counter()
    directory = Path(output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    moving_binding = file_binding(moving_result_path)
    moving = _json(moving_binding)
    if moving['status'] != 'complete' or moving['schema'] != 'P254-A32-moving-hybrid-gamma-v1':
        raise ValueError('actual complete moving-hybrid material result required')
    input_binding = moving['input_bindings']
    inputs = _json(input_binding)
    source = inputs['source']
    metadata = inputs['model_metadata']
    if metadata['source_slot'] != 0 or source['q_index'] != 0:
        raise ValueError('the published original source0/Gamma boundary is required')
    readout_binding = moving['artifacts']['readout']
    readout = _arrays(readout_binding)
    coordinates = _arrays(inputs['interaction_coordinates'])
    bank = _arrays(source['basis'])
    initial = _arrays(source['baseline0'])
    # This is the exact original PRODUCT covariance assembly, as in
    # material_clamped_array.py:79-103. No force lattice, FFT or generator needed.
    mechanical_binding = source['coefficient_source_archives']['mechanical']
    with np.load(authenticate(mechanical_binding), allow_pickle=False) as archive:
        ionic = archive['ionic_ordered']
        field = archive['field_ordered']
    nx = metadata['mechanical_size']
    if ionic.shape[0]+field.shape[0] != nx or ionic.shape[0] % 2:
        raise ValueError('saved canonical ionic/field product layout differs')
    w0 = np.zeros((nx, nx), complex)
    ni = len(ionic)
    w0[:ni, :ni] = ionic
    w0[ni:, ni:] = field
    if (not np.array_equal(readout['times_inverse_ry'], moving['times_inverse_ry'])
            or initial['rho'].shape[0] != metadata['source_count']):
        raise ValueError('actual readout time/source axes differ')
    basis_identity = 'finite-bath:'+source['basis']['sha256']
    factor_identity = 'projected-original-frame:'+source['state_authority']
    authority = source['state_authority']
    identities = dict(factor_identity=factor_identity, basis_identity=basis_identity,
                      state_authority=authority)
    print('[physics] consuming original saved duals, D and actual coordinate bank; no propagation', flush=True)
    budget_started = perf_counter()
    budget = decompose_initial_noise(readout['preparation_duals'], w0, ni//2,
        initial['rho'], bank, coordinates, source_slot=0, **identities)
    budget_seconds = perf_counter()-budget_started
    full_variances = np.diag(readout['ordered_covariance'])
    normalized = dict(mechanical={name:common_normalized(value, full_variances)
        for name,value in budget['mechanical'].items()},
        electronic_by_slot=np.stack([common_normalized(value, full_variances)
                                    for value in budget['electronic_by_slot']]),
        source0={name:common_normalized(value, full_variances)
                 for name,value in budget['source0'].items()},
        total=common_normalized(budget['total'], full_variances))
    original_manifest = _json(source['original_coefficient_manifest_snapshot'])
    mesh = original_manifest['specification']['mesh']
    ncells = int(np.prod(mesh))
    influence_started = perf_counter()
    influence = frozen_force_influence(bank, coordinates, readout['times_inverse_ry'],
        source_factor=1/np.sqrt(ncells), **identities)
    influence_seconds = perf_counter()-influence_started
    result = dict(schema='P254-A33-saved-noise-budget-and-influence-v1', status='complete',
        original_moving_result=moving_binding, original_input_bindings=input_binding,
        actual_saved_readout=readout_binding,
        actual_interaction_coordinates=inputs['interaction_coordinates'],
        actual_finite_bath=source['basis'], original_preparation=source['baseline0'],
        original_mechanical_product_covariance=mechanical_binding,
        original_manifest=source['original_coefficient_manifest_snapshot'],
        times_inverse_ry=readout['times_inverse_ry'], atom=moving['atom'], axis=moving['axis'],
        q_index=0, source_slot=0, ncells=ncells, source_count=metadata['source_count'],
        noise_budget=budget, common_normalized_contributions=normalized,
        normalization='Every complex contribution uses the SAME full-output principal-root '
                      'denominator; component magnitudes/fractions are not independently normalized or additive.',
        saved_covariance=readout['ordered_covariance'],
        raw_total_minus_saved_covariance=budget['total']-readout['ordered_covariance'],
        raw_normalized_total_minus_saved=normalized['total']-readout['normalized_gamma_coherence'],
        frozen_two_mode_influence=influence,
        interpretation='Moving initial-source budget and frozen-H force kernel are DISTINCT '
            'physical components. The kernel is not added independently atop the coupled response. '
            'No continuum/irreversible linewidth, heat, physical uncertainty bound or full-source/q calibration.',
        implementation={name:file_binding(HERE/name) for name in
            ('hybrid_noise_components.py','finite_bath_influence.py',Path(__file__).name)},
        inherited_car_consumer=file_binding(HERE.parent/'0031/export_frozen_bath_comparator.py'),
        implementation_contract=file_binding(HERE/'implementation-contract-001.json'),
        work=dict(total_seconds=perf_counter()-started, source_budget_seconds=budget_seconds,
            frozen_influence_seconds=influence_seconds,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024,
            trajectory_integrations=0, hamiltonian_actions=0, force_actions=0,
            eigendecompositions_of_H=0, state_refills=0, dft_actions=0,
            newly_built_gamma_caches=0, original_namespaces_written=False))
    binding = publish_json(directory/'result.json', result)
    print(json.dumps(dict(result=binding, work=result['work']), sort_keys=True), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--moving-result', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run(args.moving_result, args.output)
