"""A31 noise-selected finite subspace; no dynamics or physical-loss closure.

All original raw columns survive. Only a temporary force-selection proxy is
physical-row sized; neither the symmetric original basis nor a PW-square
operator is constructed. The returned state is the original raw-D state, not
the absolute-state weighting used to select additional directions.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np
from scipy.linalg import svd

from substrate_framework.physical_factor_preparation import (
    ActionFactors, _gram_root, _hermitian, _tolerance,
)


def _matrix(value, shape, name):
    # In particular, do not coerce/copy the large admitted physical-row banks.
    value = np.asarray(value)
    if value.shape != shape or not np.isfinite(value).all():
        raise ValueError(name + ' must be a finite matrix of shape ' + str(shape))
    return value


def _max_abs(value):
    return float(np.max(np.abs(value), initial=0))


def _identity_metadata(factors):
    return {key: getattr(factors, key) for key in (
        'operator_identity', 'adjoint_identity', 'units',
        'pw_basis_identity', 'source_authority')}


def build_noise_subspace(psi, gram, density_coefficients, active_coordinates,
                         hamiltonian, forces, apply_extra, *, extra_rank):
    """Return the A31 frozen-initial-H comparator's numeric coordinate bank.

    ``projected_frame[n,r]`` represents the unchanged Psi in the basis
    B=[Psi S^-1/2, extra_basis]; B itself is NOT returned or stored.
    ``projected_density = projected_frame D projected_frame^dagger`` and
    ``projected_active_coordinates = projected_frame R`` explicitly transport
    state and active projector. Original S,D,R are returned without alteration.

    The temporary selection proxy concatenates (I-Q0 Q0^dagger) A Q0 |C0|^1/2
    for EVERY forward/independent-adjoint force leg. Its full thin singular
    spectrum is retained; extra_rank is a requested count, never a tolerance
    cutoff. Degenerate/null requested directions cannot certify completeness.
    Gram/state roundoff and Pauli defects are measured, never repaired.
    apply_extra(name,E,adjoint=bool) is called only for the selected new E;
    no callback occurs at rank zero. No physical operator acts on Psi again.
    """
    psi = np.asarray(psi)
    if psi.ndim != 2 or not 0 < psi.shape[1] <= psi.shape[0]:
        raise ValueError('Psi must have independent-factor dimensions [P,r]')
    rows, rank = psi.shape
    psi = _matrix(psi, (rows, rank), 'original Psi')
    gram = _matrix(gram, (rank, rank), 'supplied original Gram')
    root, gram_values, gram_vectors = _gram_root(gram)
    inverse_root = (gram_vectors / np.sqrt(gram_values)) @ gram_vectors.conj().T
    density = _hermitian(density_coefficients, 'original raw D')
    density = _matrix(density, (rank, rank), 'original raw D')
    coordinates = np.asarray(active_coordinates)
    if coordinates.ndim != 2 or not 0 <= coordinates.shape[1] <= rank:
        raise ValueError('R must have shape [r,active_count]')
    coordinates = _matrix(coordinates, (rank, coordinates.shape[1]), 'original R')
    active_metric = coordinates.conj().T @ gram @ coordinates
    active_defect = active_metric - np.eye(coordinates.shape[1])
    if _max_abs(active_defect) > _tolerance(gram):
        raise ValueError('R^dagger S R must be identity; no projector repair')
    if not isinstance(forces, Mapping) or not forces:
        raise ValueError('nonempty force ActionFactors mapping required')
    if any(not isinstance(key, str) or not key.strip() or key == 'hamiltonian'
           for key in forces):
        raise ValueError('force names must be nonempty and not reserved hamiltonian')
    if not isinstance(hamiltonian, ActionFactors):
        raise TypeError('hamiltonian must be canonical ActionFactors')
    operators = {'hamiltonian': hamiltonian, **forces}
    for factor in operators.values():
        if not isinstance(factor, ActionFactors):
            raise TypeError('canonical ActionFactors required')
        factor.validate(psi.shape, hamiltonian.pw_basis_identity)
    maximum_extra = min(rows-rank, 2*len(forces)*rank)
    if type(extra_rank) is not int or not 0 <= extra_rank <= maximum_extra:
        raise ValueError('explicit extra_rank must lie in [0,min(P-r,2*m*r)]')
    if not callable(apply_extra):
        raise TypeError('apply_extra callback required')

    represented = root @ density @ root
    state_values, state_vectors = np.linalg.eigh(represented)
    # Absolute values ONLY choose a basis. Actual covariance consumers use D.
    absolute_state_root = (state_vectors * np.sqrt(np.abs(state_values))) @ state_vectors.conj().T
    weighting = inverse_root @ absolute_state_root
    original_overlaps = {}
    for name, factor in operators.items():
        original_overlaps[name] = tuple(psi.conj().T @ action for action in
                                       (factor.action, factor.adjoint_action))
    width = 2*len(forces)*rank
    # One temporary bank, no copied full Q0. Overwrite is permitted in LAPACK;
    # original Psi/actions are never used as SVD workspace or modified.
    proxy = np.empty((rows, width), dtype=np.complex128, order='F')
    leg_order = []
    for mode, (name, factor) in enumerate(forces.items()):
        for leg, action in enumerate((factor.action, factor.adjoint_action)):
            columns = slice((2*mode+leg)*rank, (2*mode+leg+1)*rank)
            overlap = original_overlaps[name][leg]
            proxy[:, columns] = action @ weighting
            correction = inverse_root @ (inverse_root @ overlap) @ weighting
            proxy[:, columns] -= psi @ correction
            leg_order.append((name, 'adjoint' if leg else 'forward'))
    left, singular_values, right_adjoint = svd(
        proxy, full_matrices=False, overwrite_a=True, check_finite=False,
        lapack_driver='gesdd')
    del proxy, right_adjoint
    extras = left[:, :extra_rank].copy()
    del left
    # Only the explicitly requested vectors are reorthogonalized. No hidden
    # SVD/state-eigenvalue gate changes the requested count.
    for _ in range(2):
        extras -= psi @ (inverse_root @ (inverse_root @ (psi.conj().T @ extras)))
    if extra_rank:
        extras, triangular = np.linalg.qr(extras, mode='reduced')
        if np.any(np.abs(np.diag(triangular)) <= 128*np.finfo(float).eps*max(rows, rank)):
            raise ValueError('requested singular directions do not resolve independent extras; no rank fallback')
    else:
        extras = np.empty((rows, 0), dtype=np.complex128)

    dimension = rank + extra_rank
    frame = np.zeros((dimension, rank), dtype=np.complex128)
    frame[:rank] = root
    projected_density = frame @ density @ frame.conj().T
    projected_active = frame @ coordinates
    active_projector = projected_active @ projected_active.conj().T
    overlap_extra = inverse_root @ (psi.conj().T @ extras)
    basis_metric = np.zeros((dimension, dimension), dtype=np.complex128)
    basis_metric[:rank, :rank] = inverse_root @ gram @ inverse_root
    basis_metric[:rank, rank:] = overlap_extra
    basis_metric[rank:, :rank] = overlap_extra.conj().T
    basis_metric[rank:, rank:] = extras.conj().T @ extras
    projected, operator_diagnostics = {}, {}
    callback_calls = []
    for name, factor in operators.items():
        matrices, leg_diagnostics = [], []
        for leg, action in enumerate((factor.action, factor.adjoint_action)):
            matrix = np.empty((dimension, dimension), dtype=np.complex128)
            matrix[:rank, :rank] = inverse_root @ original_overlaps[name][leg] @ inverse_root
            matrix[rank:, :rank] = (extras.conj().T @ action) @ inverse_root
            if extra_rank:
                new_action = _matrix(apply_extra(name, extras, adjoint=bool(leg)),
                                     extras.shape, name+' extra action')
                callback_calls.append((name, bool(leg), extra_rank))
                matrix[:rank, rank:] = inverse_root @ (psi.conj().T @ new_action)
                matrix[rank:, rank:] = extras.conj().T @ new_action
                extra_residual = new_action - psi @ (inverse_root @ matrix[:rank, rank:])
                extra_residual -= extras @ matrix[rank:, rank:]
                extra_residual_norm = float(np.linalg.norm(extra_residual))
                del new_action, extra_residual
            else:
                extra_residual_norm = 0.0
            # Measured residual on original columns, without a Q0 bank.
            residual = action @ inverse_root
            residual -= psi @ (inverse_root @ matrix[:rank, :rank])
            residual -= extras @ matrix[rank:, :rank]
            leg_diagnostics.append(dict(
                original_action_projection_residual_frobenius=float(np.linalg.norm(residual)),
                extra_action_projection_residual_frobenius=extra_residual_norm))
            del residual
            matrices.append(matrix)
        projected[name] = matrices
        operator_diagnostics[name] = dict(
            forward=leg_diagnostics[0], adjoint=leg_diagnostics[1],
            supplied_adjoint_defect_frobenius=float(np.linalg.norm(matrices[1]-matrices[0].conj().T)),
            forward_hermiticity_frobenius=float(np.linalg.norm(matrices[0]-matrices[0].conj().T)))
    discarded = singular_values[extra_rank:]
    diagnostics = dict(
        gram_eigenvalues=gram_values,
        supplied_gram_hermiticity_frobenius=float(np.linalg.norm(gram-gram.conj().T)),
        basis_metric=basis_metric,
        basis_metric_identity_max_abs=_max_abs(basis_metric-np.eye(dimension)),
        original_extra_overlap_frobenius=float(np.linalg.norm(overlap_extra)),
        original_raw_D_eigenvalues=np.linalg.eigvalsh(density),
        state_eigenvalues=state_values,
        state_hermiticity_frobenius=float(np.linalg.norm(represented-represented.conj().T)),
        state_pauli_roundoff_allowance=float(_tolerance(represented)),
        state_within_pauli_roundoff=bool(state_values.min() >= -_tolerance(represented)
                                       and state_values.max() <= 1+_tolerance(represented)),
        state_transfer_frobenius=float(np.linalg.norm(projected_density[:rank, :rank]-represented)),
        active_metric_identity_max_abs=_max_abs(active_defect),
        active_projector_idempotence_frobenius=float(np.linalg.norm(active_projector@active_projector-active_projector)),
        state_active_commutator_frobenius=float(np.linalg.norm(projected_density@active_projector-active_projector@projected_density)),
        singular_values=singular_values,
        selected_rank=extra_rank,
        discarded_proxy_frobenius=float(np.linalg.norm(discarded)),
        discarded_proxy_spectral=float(discarded[0]) if len(discarded) else 0.0,
        operators=operator_diagnostics)
    metadata = dict(
        schema='P254-noise-subspace-v1', full_pw_dimension=rows,
        original_rank=rank, extra_rank_requested=extra_rank, extra_rank_selected=extra_rank,
        dimension=dimension, active_count=coordinates.shape[1],
        basis='B=[Psi S^-1/2,E]; no copied Q0 bank returned',
        projected_frame='B^dagger Psi=[S^1/2;0]; coordinates of ALL original raw columns',
        state='rho_B=projected_frame D projected_frame^dagger; original raw D unchanged',
        active='Qactive_B=projected_frame R; Pactive_B=Qactive_B Qactive_B^dagger',
        gram_authority='supplied admitted original S; original Gram is not recomputed',
        selection='SVD of complementary forward/true-adjoint force actions weighted by |S^1/2 D S^1/2|^1/2',
        proxy_leg_order=leg_order, selection_proxy_shape=(rows, width),
        rank_policy='explicit requested count; no tolerance-selected rank or state clipping',
        discarded_proxy_scope='initial basis-selection proxy only; not a late-time error or physical tail bound',
        original_physical_action_replays=0, extra_callback_calls=callback_calls,
        operators={name: _identity_metadata(factor) for name, factor in operators.items()},
        scope='selected-mode finite-basis nonstationary frozen-initial-H comparator, not full moving/RPA loss closure')
    return dict(extra_basis=extras, inverse_sqrt_gram=inverse_root,
                projected_frame=frame, gram=gram, density_coefficients=density,
                active_coordinates=coordinates, projected_density=projected_density,
                projected_active_coordinates=projected_active, active_projector=active_projector,
                hamiltonian=projected['hamiltonian'][0],
                hamiltonian_adjoint=projected['hamiltonian'][1],
                forces={name: projected[name][0] for name in forces},
                force_adjoints={name: projected[name][1] for name in forces},
                singular_values=singular_values, state_eigenvalues=state_values,
                diagnostics=diagnostics, metadata=metadata)
