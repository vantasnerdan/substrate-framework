"""Actual saved-data original-output/omitted-forcing cross contraction.

This defines a conditional forced-output comparator. It does not turn the
frozen force-vertex split into dynamically matched passive feedback.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter

from cached_bath_inputs import HERE
from export_frozen_bath_comparator import car_moments, projected_factors, save_json
from export_output_directed_actions import restore_filter
from homogeneous_electronic_output import homogeneous_electronic_operators
import numpy as np
from substrate_framework.streamed_material_state import authenticate, file_binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bath-result',required=True)
    parser.add_argument('--filter-packet',required=True)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    started = perf_counter()
    result_binding = file_binding(args.bath_result)
    if result_binding['sha256'] != '3d2fdab5617af7c9b7fe8e6afdd3cde6ea30e074cad4e2ad481ca86875e57d93':
        raise ValueError('the settled actual A31 bath result is required')
    bath = json.loads(Path(args.bath_result).read_bytes())
    filter_,filter_binding,filter_archive = restore_filter(args.filter_packet)
    packet = json.loads(Path(args.filter_packet).read_bytes())
    source = bath['source']
    if (packet['filter']['q_index'] != 0
            or filter_binding['sha256'] != bath['filter_packet']['sha256']
            or source['coefficient_manifest']['sha256'] != packet['source']['original_coefficient_manifest_snapshot']['sha256']
            or source['original_boundary']['sha256'] != packet['source']['original_boundary']['sha256']):
        raise ValueError('same Gamma filter,original state and coefficient-manifest authority required')
    kind,slot_text = source['coefficient_item'].split('/')
    if kind != 'slot':
        raise ValueError('saved coefficient item must identify the source slot')
    slot = int(slot_text)
    keys = ('projected_frame','gram','density_coefficients','projected_active_coordinates','projected_density')
    with np.load(authenticate(bath['basis']),allow_pickle=False) as archive:
        bank = {key:archive[key] for key in keys}
    with np.load(authenticate(bath['filtered_operators']),allow_pickle=False) as archive:
        omitted = archive['operator_omitted_vertex']
        omitted_adjoint = archive['adjoint_operator_omitted_vertex']
    active = bank['projected_active_coordinates']
    nx = 2*packet['filter']['ndof']+2*packet['filter']['nfield']
    ns,remainder = divmod(filter_.initial_dual.shape[0]-nx,active.shape[1]**2)
    if remainder:
        raise ValueError('saved initial dual does not have the admitted coupled layout')
    homogeneous = homogeneous_electronic_operators(filter_.initial_dual,active,
        mechanical_size=nx,source_count=ns,source_slot=slot)[0]
    matrices = dict(homogeneous_electronic=homogeneous,omitted_forced=omitted,
                    conditional_sum=homogeneous+omitted)
    adjoints = dict(homogeneous_electronic=homogeneous.conj().T,
                   omitted_forced=omitted_adjoint,
                   conditional_sum=homogeneous.conj().T+omitted_adjoint)
    basis_identity = 'A31 finite-bath coordinate basis '+bath['basis']['sha256']
    moments = car_moments(bank['projected_frame'],bank['gram'],bank['density_coefficients'],
        projected_factors(matrices,adjoints,bank['projected_frame'],basis_identity,
                          bath['original_state_authority'],'bohr'),
        bath['original_factor_identity'],basis_identity,bath['original_state_authority'])
    covariance = {(row['left'],row['right']):row['value'] for row in moments['covariance']}
    h,e = 'homogeneous_electronic','omitted_forced'
    difference = covariance['conditional_sum','conditional_sum']-covariance[h,h]
    reconstruction = difference-covariance[e,e]-covariance[h,e]-covariance[e,h]
    with np.load(authenticate(source['coefficient_archive']),allow_pickle=False) as archive:
        saved_rho = archive['rho0']
    compression = active.conj().T@bank['projected_density']@active
    output = dict(schema='P254-original-homogeneous-output-cross-v1',status='complete',
        bath_result=result_binding,filter_packet=filter_binding,filter_archive=filter_archive,
        source_identity=source['original_metadata']['identity'],saved_coefficient_item=source['coefficient_item'],
        original_source_record=source['original_record'],source_slot=slot,source_count=ns,
        original_state_authority=bath['original_state_authority'],basis=bath['basis'],
        filtered_operators=bath['filtered_operators'],coefficient_archive=source['coefficient_archive'],
        homogeneous_operator='Qa @ Zslot.T @ Qa†;ordinarytranspose;already pulled back;not repropagated',
        original_preparation_compression_defect_frobenius=float(np.linalg.norm(compression-saved_rho)),
        covariance=moments,conditional_added_variance=difference,
        raw_added_variance_minus_omitted_and_both_cross_terms=reconstruction,
        normalization='No extra homogeneous Fourier/source/spin factor;omitted forcing already hasF/sqrtN;no final realspace1/N.',
        elapsed_seconds=perf_counter()-started,
        implementation=[file_binding(HERE/name) for name in ('homogeneous_electronic_output.py',
            'consume_homogeneous_output_cross.py','export_frozen_bath_comparator.py')],
        scope='One original source0/Gamma electronic component of reconstructed movingactive output plus frozenH/twomode/8extra-direction omittedforce forcing,in unchanged originalstate. Includes both ordered original-output/source cross terms. NOT dynamically matched complementfeedback,total ionic/field/allsource covariance,corrected A29 spatial coherence,physical loss,lifetime or completecampaign.')
    binding = save_json(Path(args.output),output)
    print(json.dumps(dict(result=binding,elapsed_seconds=output['elapsed_seconds']),allow_nan=False))


if __name__ == '__main__':
    main()
