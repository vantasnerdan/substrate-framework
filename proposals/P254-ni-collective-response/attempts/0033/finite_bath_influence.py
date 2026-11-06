"""Two-time electronic influence of the paid frozen-initial-H finite bath.

Keep the original nonstationary frame/D and independent physical force legs.
No stationary spectrum, FDT, Markov collision rate or moving-history recovery
is inferred. Raw H and H-adjoint similarity evolutions follow the existing A31
comparator; their defects are exposed, not repaired into physical unitarity.
"""
import numpy as np

from export_frozen_bath_comparator import car_moments, projected_factors
from substrate_framework.nonstationary_force_kernel import OrderedForcePair


def frozen_force_influence(bank, coordinate_bank, times, *, source_factor,
        factor_identity, basis_identity, state_authority):
    """Compute the actual two-mode ordered noise and -uF causal commutator.

    source_factor multiplies each force ONCE before the CAR contraction. The
    material caller supplies the existing canonical 1/sqrtN, without an extra
    spin or source measure. Instantaneous scalar/contact derivatives are NOT
    reconstructed from this noninstantaneous response.
    """
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or not len(times) or not np.isfinite(times).all()
            or not np.isfinite(source_factor) or source_factor <= 0):
        raise ValueError('finite saved times and explicit positive canonical source factor required')
    energies = coordinate_bank['eigenvalues']
    v = coordinate_bank['right_eigenvectors']
    inverse = coordinate_bank['inverse_right_eigenvectors']
    h = bank['hamiltonian']
    if (v.shape != h.shape or inverse.shape != h.shape
            or energies.shape != (len(h),)):
        raise ValueError('literal-H persisted eigencoordinate bank required')
    forces = np.stack([bank['force_mode0'], bank['force_mode1']])*source_factor
    adjoints = np.stack([bank['adjoint_force_mode0'], bank['adjoint_force_mode1']])*source_factor
    eigen_forces = inverse @ forces @ v
    eigen_adjoints = v.conj().T @ adjoints @ inverse.conj().T
    operators, operator_adjoints = {}, {}
    for ti, time in enumerate(times):
        phase = np.exp(1j*energies*time)
        undo = np.exp(-1j*energies*time)
        phase_adj = np.exp(1j*energies.conj()*time)
        undo_adj = np.exp(-1j*energies.conj()*time)
        for mode in range(2):
            name = f'{ti}:{mode}'
            operators[name] = v @ (phase[:, None]*eigen_forces[mode]*undo[None, :]) @ inverse
            operator_adjoints[name] = inverse.conj().T @ (
                phase_adj[:, None]*eigen_adjoints[mode]*undo_adj[None, :]) @ v.conj().T
    factors = projected_factors(operators, operator_adjoints, bank['projected_frame'],
        basis_identity, state_authority, 'Ry/bohr, canonical Gamma force')
    moments = car_moments(bank['projected_frame'], bank['gram'],
        bank['density_coefficients'], factors, factor_identity, basis_identity,
        state_authority)
    values = {(entry['left'], entry['right']): entry['value']
              for entry in moments['covariance']}
    shape = (len(times), 2, len(times), 2)
    ordered = np.empty(shape, complex)
    symmetrized = np.empty(shape, complex)
    commutator = np.empty(shape, complex)
    retarded = np.empty(shape, complex)
    for ti, time in enumerate(times):
        for mode in range(2):
            for tj, second_time in enumerate(times):
                for second_mode in range(2):
                    a, b = f'{ti}:{mode}', f'{tj}:{second_mode}'
                    pair = OrderedForcePair(float(time), float(second_time),
                        moments['means'][a], moments['means'][b], values[a, b], values[b, a])
                    index = ti, mode, tj, second_mode
                    ordered[index] = pair.ordered_ab
                    symmetrized[index] = pair.symmetrized
                    commutator[index] = pair.commutator
                    retarded[index] = pair.retarded(perturbation_sign=-1, equal_time_theta=.5)
    phase = np.exp(-1j*energies*times[-1])
    propagator = (v*phase) @ inverse
    return dict(times_inverse_ry=times, ordered=ordered,
        reverse_ordered=ordered.transpose(2, 3, 0, 1),
        symmetrized=symmetrized, commutator=commutator, retarded=retarded,
        means=moments['means'], original_state_admissibility=moments['admissibility'],
        diagnostics=dict(
            original_state_commutator_frobenius=float(np.linalg.norm(
                h@bank['projected_density']-bank['projected_density']@h)),
            raw_hamiltonian_hermiticity_frobenius=float(np.linalg.norm(h-h.conj().T)),
            supplied_force_adjoint_defect_by_mode=np.linalg.norm(
                adjoints-forces.conj().transpose(0, 2, 1), axis=(1, 2)),
            final_raw_propagator_unitarity_defect_frobenius=float(np.linalg.norm(
                propagator.conj().T@propagator-np.eye(len(h)))),
            maximum_evolved_adjoint_defect=float(max(np.max(np.abs(
                operator_adjoints[name]-operators[name].conj().T)) for name in operators))),
        normalization='each force has the supplied canonical source_factor once; no extra spin/cell measure',
        source_factor=source_factor, equal_time_theta=.5, perturbation_sign=-1,
        instantaneous_scalar_contact_terms_included=False,
        scope='source0/two paid force modes in the frozen literal-H finite bath; '
              'NONSTATIONARY two-time kernel, not a stationary positive spectrum, '
              'moving/RPA replacement, all-source/q noise, irreversible linewidth or heat')
