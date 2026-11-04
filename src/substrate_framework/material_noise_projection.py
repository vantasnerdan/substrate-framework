"""Original-preparation retained/omitted electronic covariance at t=0.

For the caller-authenticated full raw frame Psi, supplied Gram S and raw density
coefficients D, rho=Psi D Psi^dagger has an exactly empty orthogonal complement.
The supplied coordinates R define Q=Psi R, R^dagger S R=I, P=Q Q^dagger. Each
operator is split as F_R=P F P and F_E=F-F_R, WITHOUT restricting the state.
All directed blocks use the same canonical CAR contraction
C_AB=Tr(D T_AB)-Tr(D M_A D M_B), including initially empty PW states.

Only raw-factor-sized moments are formed. In particular, neither P nor a bank
of projected full-row actions is constructed. Independent F^dagger Psi inputs
are used even for non-Hermitian F; no overlap is repaired to match its adjoint.

Sources: physical_factor_preparation.py (FactorFunctional, FactorPreparation),
material_force_fluctuations.py (row-bounded full-PW cross-Gram contraction), and
P254 attempts/0029/reduced_active_basis.py (Q=Psi R coordinate convention).
This is one unweighted quasifree sector at the original time, not an all-source
or full-q result, a two-time kernel, a stationary bath, damping or a lifetime.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .material_force_fluctuations import _cross_gram
from .physical_factor_preparation import ActionFactors, FactorFunctional, FactorPreparation


_ROW_BLOCK_SIZE = 2048
_PARTS = ('retained', 'omitted')
_IDENTITY_FIELDS = ('operator_identity', 'adjoint_identity', 'units',
                    'pw_basis_identity', 'source_authority')


def _frame(value, shape, name):
    """Validate full-row inputs without making a full-frame dtype/layout copy."""
    value = np.asarray(value)
    if value.ndim != 2 or value.shape != shape or value.dtype.kind not in 'fc':
        raise ValueError(name+' must retain every full-PW row and raw source column')
    for start in range(0, shape[0], _ROW_BLOCK_SIZE):
        if not np.isfinite(value[start:start+_ROW_BLOCK_SIZE]).all():
            raise ValueError(name+' must be finite')
    return value


def _coordinates(value, shape, name):
    value = np.asarray(value, dtype=np.complex128)
    if value.ndim != 2 or value.shape != shape or not np.isfinite(value).all():
        raise ValueError(name+' has invalid dimensions or nonfinite entries')
    return value


def _operator_metadata(probe, basis_identity):
    if not isinstance(probe, ActionFactors):
        raise TypeError('probes must contain ActionFactors with actual independent adjoints')
    metadata = {key: getattr(probe, key) for key in _IDENTITY_FIELDS}
    if any(not isinstance(value, str) or not value.strip() for value in metadata.values()):
        raise ValueError('explicit operator, adjoint, units, basis and source identities required')
    if metadata['pw_basis_identity'] != basis_identity:
        raise ValueError('operator and factor PW basis identities differ')
    return metadata


def _moments(psi, coordinates, left, right, active_metric, probes, basis_identity):
    """Full moments and low-rank retained legs; never construct Q or PFP Psi."""
    overlaps, legs, metadata, actions = {}, {}, {}, {}
    for name, probe in probes.items():
        metadata[name] = _operator_metadata(probe, basis_identity)
        action = _frame(probe.action, psi.shape, name+' action')
        adjoint = action if probe.adjoint_action is probe.action else _frame(
            probe.adjoint_action, psi.shape, name+' adjoint action')
        full = _cross_gram(psi, action, _ROW_BLOCK_SIZE)
        adjoint_overlap = full if adjoint is action else _cross_gram(
            psi, adjoint, _ROW_BLOCK_SIZE)
        active_action = coordinates.conj().T @ full
        active_adjoint = coordinates.conj().T @ adjoint_overlap
        retained_vertex = active_action @ coordinates
        retained_adjoint = active_adjoint @ coordinates
        overlaps[name, 'full'] = full
        overlaps[name, 'retained'] = right @ retained_vertex @ left
        overlaps[name, 'omitted'] = full-overlaps[name, 'retained']
        # The independent adjoint leg is not inferred from full or its dagger.
        legs[name] = (active_action, active_adjoint.conj().T,
                      retained_vertex @ left,
                      left.conj().T @ retained_adjoint.conj().T)
        actions[name] = (action, adjoint)
    products = {}
    for a, (_, adjoint) in actions.items():
        _, full_to_active, _, retained_left = legs[a]
        for b, (action, _) in actions.items():
            active_to_full, _, retained_right, _ = legs[b]
            ff = _cross_gram(adjoint, action, _ROW_BLOCK_SIZE)
            fr = full_to_active @ retained_right
            rf = retained_left @ active_to_full
            rr = retained_left @ active_metric @ retained_right
            products[(a, 'full'), (b, 'full')] = ff
            products[(a, 'retained'), (b, 'retained')] = rr
            products[(a, 'retained'), (b, 'omitted')] = rf-rr
            products[(a, 'omitted'), (b, 'retained')] = fr-rr
            products[(a, 'omitted'), (b, 'omitted')] = ff-fr-rf+rr
    return overlaps, products, metadata


def projected_noise_decomposition(psi, gram, density_coefficients,
        active_coordinates, probes, *, factor_identity, pw_basis_identity,
        state_authority):
    """Return full and four ordered retained/omitted covariance blocks.

    Parameters are the ORIGINAL Psi[P,r], admitted S[r,r], raw D[r,r],
    R[r,m] and Mapping[str, ActionFactors] on all P rows and r columns. The
    caller authenticates S as the Gram of this unchanged Psi; it is deliberately
    not recomputed. R is already selected in the original saved coordinate
    convention, not a set of principal raw-D indices. An empty selection m=0
    is allowed; 0 <= m <= r. R^dagger S R must be identity within the declared
    roundoff allowance, otherwise the input does not define the requested P.
    FactorPreparation performs the canonical Gram/D admissibility checks.

    Return value is a dict with:
      full[(a,b)]: <delta F_a delta F_b> (not F_a^dagger F_b);
      blocks[(a,part_a,b,part_b)]: same ordered covariance of split operators;
      reconstruction_defect[(a,b)]: full minus the raw sum of all four blocks;
      means[(a,part)]: mean, with part 'full', 'retained' or 'omitted';
      metadata: state/operator/basis identities, units, raw projector diagnostics.
    Both mixed orders are computed, never assumed zero or conjugated. Original
    weights, spin multiplicity and cell normalization belong to the caller.
    No clipping, Hermitization, refill, scalar noise fit or time evolution occurs.
    Supplied action/adjoint identities are caller authority, not authenticated
    from these finitely many columns. Returned diagnostics are not repairs.
    """
    for identity in (factor_identity, pw_basis_identity, state_authority):
        if not isinstance(identity, str) or not identity.strip():
            raise ValueError('explicit original factor, PW-basis and state identities required')
    psi = np.asarray(psi)
    if psi.ndim != 2 or not 0 < psi.shape[1] <= psi.shape[0]:
        raise ValueError('nonempty full-PW raw frame with independent-factor dimensions required')
    psi = _frame(psi, psi.shape, 'original raw frame')
    rank = psi.shape[1]
    gram = _coordinates(gram, (rank, rank), 'original admitted Gram')
    coordinates = np.asarray(active_coordinates, dtype=np.complex128)
    if coordinates.ndim != 2 or not 0 <= coordinates.shape[1] <= rank:
        raise ValueError('active coordinates must have shape [raw factor count, active count]')
    coordinates = _coordinates(coordinates, (rank, coordinates.shape[1]), 'active coordinates')
    if not isinstance(probes, Mapping) or not probes:
        raise ValueError('nonempty mapping of original physical probe actions required')
    if any(not isinstance(name, str) or not name.strip() for name in probes):
        raise ValueError('nonempty string probe names required')

    left = coordinates.conj().T @ gram  # Q^dagger Psi
    right = gram @ coordinates         # Psi^dagger Q
    metric = coordinates.conj().T @ right
    metric_defect = metric-np.eye(coordinates.shape[1])
    projector_allowance = (128*np.finfo(float).eps*max(1, rank)
                           * (1+np.max(np.abs(metric), initial=0)))
    if np.max(np.abs(metric_defect), initial=0) > projector_allowance:
        raise ValueError('R^dagger S R is not identity within roundoff; no projector repair')
    overlaps, products, operators = _moments(
        psi, coordinates, left, right, metric, probes, pw_basis_identity)
    part_metadata = {(name, part): dict(value, operator_part=part)
                     for name, value in operators.items() for part in ('full',)+_PARTS}
    functional = FactorFunctional(gram, overlaps, products, part_metadata,
                                  factor_identity, pw_basis_identity, psi.shape[0])
    preparation = FactorPreparation(functional, density_coefficients, state_authority)
    full, blocks, defects = {}, {}, {}
    for a in probes:
        for b in probes:
            full[a, b] = preparation.ordered((a, 'full'), (b, 'full'))
            for part_a in _PARTS:
                for part_b in _PARTS:
                    blocks[a, part_a, b, part_b] = preparation.ordered(
                        (a, part_a), (b, part_b))
            defects[a, b] = full[a, b]-sum(
                blocks[a, part_a, b, part_b]
                for part_a in _PARTS for part_b in _PARTS)
    means = {(name, part): preparation.mean((name, part))
             for name in probes for part in ('full',)+_PARTS}
    metadata = {
        'factor_identity': factor_identity, 'pw_basis_identity': pw_basis_identity,
        'state_authority': state_authority, 'operators': operators,
        'full_pw_dimension': psi.shape[0], 'raw_factor_count': rank,
        'active_count': coordinates.shape[1], 'admissibility': preparation.admissibility,
        'state': 'original rho=Psi D Psi^dagger; no restriction, refill or change of D',
        'gram_authority': 'caller-supplied admitted same-frame Gram; not recomputed',
        'projector': 'Q=Psi R; P=Q Q^dagger; retained=P F P; omitted=F-P F P',
        'projector_diagnostics': {
            'active_metric_identity_max_abs': float(np.max(np.abs(metric_defect), initial=0)),
            'active_metric_identity_frobenius': float(np.linalg.norm(metric_defect)),
            'active_metric_hermiticity_frobenius': float(np.linalg.norm(metric-metric.conj().T)),
            'coordinate_idempotence_defect_frobenius': float(np.linalg.norm(
                coordinates @ metric_defect @ coordinates.conj().T)),
            'roundoff_allowance_max_abs': float(projector_allowance)},
        'ordered_covariance': 'Tr(D T_AB)-Tr(D M_A D M_B); no conjugation of B',
        'adjoints': 'independent supplied adjoint actions, including for retained legs',
        'normalization': 'one unweighted single sector; no spin, source measure or cell factor',
        'covariance_units': 'product of the two declared probe units',
        'reconstruction': 'raw full-minus-sum-four-blocks; no clipping or symmetrization',
        'scope': 'initial electronic quasifree component only; not full128/64q, a late-time loss or lifetime',
        'source_references': [
            'src/substrate_framework/physical_factor_preparation.py:FactorFunctional,FactorPreparation',
            'src/substrate_framework/material_force_fluctuations.py:_cross_gram',
            'proposals/P254-ni-collective-response/attempts/0029/reduced_active_basis.py:prepare_active_basis']}
    return dict(full=full, blocks=blocks, reconstruction_defect=defects,
                means=means, metadata=metadata)
