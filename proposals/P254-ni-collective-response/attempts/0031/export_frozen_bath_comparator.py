"""Actual source0/two-mode frozen-initial-H NONSTATIONARY comparator.

All346 original state columns survive. Eight new empty basis directions are
an explicit approximation, not an accuracy criterion. No full moving/RPA
closure, all-source noise, damping rate or corrected lattice coherence follows.
"""
import argparse
import json
from pathlib import Path
import resource
from time import perf_counter

from cached_bath_inputs import CachedBathInputs, HERE
from clamped_scientific_readout import encode
from filtered_quadratic_bath import integrate_filtered_operators
from noise_subspace import build_noise_subspace
import numpy as np
from substrate_framework.physical_factor_preparation import ActionFactors, FactorFunctional, FactorPreparation
from substrate_framework.streamed_material_state import file_binding


def save_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(encode(value), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    return file_binding(path)


def save_npz(path, arrays):
    with Path(path).open('xb') as stream:
        np.savez(stream, **arrays)
    return file_binding(path)


def car_moments(frame, gram, density, factors, factor_identity, basis_identity, authority):
    """Canonical CAR consumer, retaining supplied S and both true action legs."""
    for probe in factors.values():
        probe.validate(frame.shape, basis_identity)
    overlaps = {name: frame.conj().T@probe.action for name, probe in factors.items()}
    products = {(a,b): left.adjoint_action.conj().T@right.action
                for a,left in factors.items() for b,right in factors.items()}
    metadata = {name: {key: getattr(probe,key) for key in ('operator_identity',
        'adjoint_identity','units','pw_basis_identity','source_authority')}
        for name,probe in factors.items()}
    functional = FactorFunctional(gram, overlaps, products, metadata,
                                  factor_identity, basis_identity, frame.shape[0])
    preparation = FactorPreparation(functional, density, authority)
    return dict(means={name: preparation.mean(name) for name in factors},
        covariance=[dict(left=a,right=b,value=preparation.ordered(a,b))
                    for a in factors for b in factors],
        admissibility=preparation.admissibility)


def projected_factors(operators, adjoints, frame, basis_identity, authority, units):
    return {name: ActionFactors(matrix@frame, adjoints[name]@frame,
        name, name+' independently supplied adjoint', units, basis_identity, authority)
        for name,matrix in operators.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--action-packet', required=True)
    parser.add_argument('--filter-packet', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--extra-rank', type=int, required=True)
    args = parser.parse_args()
    from export_full_noise_probe import runtime
    runtime()
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    inputs = CachedBathInputs(args.action_packet,args.filter_packet,directory,args.extra_rank)
    context = inputs.context
    print('[physics] second physical mode saved;cached literal H/mode0 not replayed',flush=True)
    basis_started = perf_counter()
    bank = build_noise_subspace(context.psi,context.gram,context.density_coefficients,
        context.active_coordinates,inputs.hamiltonian,inputs.forces,inputs.apply_extra,
        extra_rank=args.extra_rank)
    basis_seconds = perf_counter()-basis_started
    bank_arrays = {key:bank[key] for key in ('extra_basis','inverse_sqrt_gram',
        'projected_frame','gram','density_coefficients','active_coordinates',
        'projected_density','projected_active_coordinates','active_projector',
        'hamiltonian','hamiltonian_adjoint','singular_values','state_eigenvalues')}
    bank_arrays['basis_metric'] = bank['diagnostics']['basis_metric']
    for name in inputs.forces:
        bank_arrays['force_'+name] = bank['forces'][name]
        bank_arrays['adjoint_force_'+name] = bank['force_adjoints'][name]
    basis_binding = save_npz(directory/'finite-bath-basis.npz',bank_arrays)
    basis_identity = 'A31 finite-bath coordinate basis '+basis_binding['sha256']
    diagnostics = dict(bank['diagnostics'])
    del diagnostics['basis_metric']
    save_json(directory/'basis-complete.json',dict(status='complete',basis=basis_binding,
        metadata=bank['metadata'],diagnostics=diagnostics,seconds=basis_seconds,
        new_physical_actions=inputs.extra_calls))
    print('[physics] finite bath basis saved;all original state columns retained',flush=True)
    initial_full = car_moments(context.psi,context.gram,context.density_coefficients,
        inputs.forces,context.factor_identity,context.pw_basis_identity,context.state_authority)
    initial_projected = car_moments(bank['projected_frame'],bank['gram'],bank['density_coefficients'],
        projected_factors(bank['forces'],bank['force_adjoints'],bank['projected_frame'],
                          basis_identity,context.state_authority,'Ry/bohr'),
        context.factor_identity,basis_identity,context.state_authority)
    save_json(directory/'initial-force-projection-complete.json',dict(
        original_full_pw=initial_full,finite_bath=initial_projected,
        scope='Initial two-mode unweighted force-noise projection comparison only;not a late error bound'))
    projector = bank['active_projector']
    names = list(inputs.forces)
    force = np.stack([bank['forces'][name] for name in names])
    adjoint = np.stack([bank['force_adjoints'][name] for name in names])
    retained = np.stack([projector@matrix@projector for matrix in force])
    retained_adjoint = np.stack([projector@matrix@projector for matrix in adjoint])
    # Tangent p'_q source is F_k/sqrt(N);CAR adds no further sector/spin/cell weight.
    # Final real-space inverse-Fourier1/N and the other127 sources are NOT applied.
    coefficients = inputs.modal_coefficients
    batched = np.zeros(coefficients.shape[:2]+(4,2),dtype=np.complex128)
    batched[:,:,:2,0] = coefficients[:,:,:,0]
    batched[:,:,2:,1] = coefficients[:,:,:,0]
    integration_started = perf_counter()
    filtered = integrate_filtered_operators(bank['hamiltonian'],
        np.concatenate((force,retained))/np.sqrt(inputs.ncells),
        np.concatenate((adjoint,retained_adjoint))/np.sqrt(inputs.ncells),
        inputs.segment_bounds,batched)
    integration_seconds = perf_counter()-integration_started
    matrices = dict(full_vertex=filtered['operators'][0],
                    retained_vertex=filtered['operators'][1],
                    omitted_vertex=filtered['operators'][0]-filtered['operators'][1])
    matrix_adjoints = dict(full_vertex=filtered['adjoint_operators'][0],
                          retained_vertex=filtered['adjoint_operators'][1],
                          omitted_vertex=filtered['adjoint_operators'][0]-filtered['adjoint_operators'][1])
    arrays = {'operator_'+name:value for name,value in matrices.items()}
    arrays.update({'adjoint_operator_'+name:value for name,value in matrix_adjoints.items()})
    arrays.update({'spectral_'+key:value for key,value in filtered['spectral'].items()})
    operator_binding = save_npz(directory/'filtered-operators.npz',arrays)
    moments = car_moments(bank['projected_frame'],bank['gram'],bank['density_coefficients'],
        projected_factors(matrices,matrix_adjoints,bank['projected_frame'],basis_identity,
                          context.state_authority,'bohr'),
        context.factor_identity,basis_identity,context.state_authority)
    cov = {(row['left'],row['right']):row['value'] for row in moments['covariance']}
    reconstruction = cov['full_vertex','full_vertex']-sum(cov[a,b]
        for a in ('retained_vertex','omitted_vertex') for b in ('retained_vertex','omitted_vertex'))
    h,rho = bank['hamiltonian'],bank['projected_density']
    result = dict(schema='P254-A31-actual-frozen-bath-comparator-v1',status='complete',
        original_action_packet=inputs.action_binding,filter_packet=inputs.filter_binding,
        filter_archive=inputs.filter_archive,new_mode_packet=file_binding(directory/'mode1-complete.json'),
        basis=basis_binding,filtered_operators=operator_binding,
        source=context.source,original_state_authority=context.state_authority,
        original_factor_identity=context.factor_identity,
        original_pw_basis_identity=context.pw_basis_identity,
        basis_metadata=bank['metadata'],basis_diagnostics=diagnostics,
        initial_unweighted_force_covariance=dict(full_pw=initial_full,finite_bath=initial_projected),
        filtered_covariance=moments,raw_full_minus_retained_omitted_sum=reconstruction,
        raw_state_commutator_frobenius=float(np.linalg.norm(h@rho-rho@h)),
        integration_metadata=filtered['metadata'],integration_diagnostics=filtered['diagnostics'],
        normalization=dict(mesh=inputs.mesh,ncells=inputs.ncells,force_insertion='unweighted F_k/sqrt(N)',
            further_source_weight_applied=False,spin_factor_applied=False,
            final_real_space_inverse_fourier_applied=False),
        work=dict(loading_seconds=inputs.loading_seconds,basis_seconds=basis_seconds,
            integration_seconds=integration_seconds,total_seconds=perf_counter()-started,
            peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            new_physical_actions=inputs.extra_calls,cached_original_H_mode0_replays=0),
        implementation=[file_binding(HERE/name) for name in ('cached_bath_inputs.py',
            'export_frozen_bath_comparator.py','noise_subspace.py','filtered_quadratic_bath.py')],
        scope='Original source0/Gamma;two selected ionic force modes;8 extra empty directions;'
            'frozen literal initial Hamiltonian;NONSTATIONARY unchanged original state. '
            'Initial retained/omitted FORCE-VERTEX split propagated in the same finite bath;'
            'both cross orders retained. Conditional forced-output electronic contribution, '
            'NOT corrected A29 total coherence,full moving/RPA response,all128 sources/64q, '
            'field/collision loss,accuracy certificate,physical damping,lifetime or campaign completion.')
    binding = save_json(directory/'result.json',result)
    print(json.dumps(dict(result=binding,total_seconds=result['work']['total_seconds']),allow_nan=False))


if __name__ == '__main__':
    main()
