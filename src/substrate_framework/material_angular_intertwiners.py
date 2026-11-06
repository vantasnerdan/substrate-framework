"""Declared ordinary-NC s/p/d intertwiners, not fitted projector covariance.

For column Cartesian momenta pt = R ps, row channel arrays obey
B(pt) = B(ps) Q(R). This is precisely Bt = C Q + residual in
material_operator_covariance.covariance_bound. Q is real orthogonal, including
improper operations, and never mixes original radial-projector identities.
Translations, atomic phases, rounded geometry and reciprocal row maps are NOT
part of Q. Their defects remain in the complete covariance bound.

The polynomial basis is qe_material_field_derivatives.solid_harmonics and
qe_nc_projectors.qe_real_harmonics: s; (z,-x,-y); then
((2z²-x²-y²)/2,-sqrt(3)xz,-sqrt(3)yz,sqrt(3)(x²-y²)/2,sqrt(3)xy),
with the common normalization for each l suppressed in the algebra below.
No material data, projector fit, D diagonalization or numerical orthogonalization
enters this construction. The tensor contraction is an analytic coefficient
identity, not least squares on sampled directions.
"""
from __future__ import annotations

import numpy as np


def signed_permutation(rotation):
    """Validate an exactly represented Cartesian signed permutation (det ±1)."""
    r = np.asarray(rotation)
    if (r.shape != (3, 3) or r.dtype.kind not in 'iuf'
            or not np.isfinite(r).all() or not np.isin(r, (-1, 0, 1)).all()
            or not np.array_equal(r.T @ r, np.eye(3))):
        raise ValueError('exact Cartesian signed-permutation rotation required')
    return r.astype(np.int64)


def harmonic_intertwiner(rotation, ell):
    """Q_l with Y_l(p @ R.T) = Y_l(p) @ Q_l in QE ylmr2 order.

    Q(R2 R1) = Q(R1) Q(R2); this right-action convention is intentional.
    The result is floating evaluation of exact rational/sqrt(3) coefficients.
    """
    r = signed_permutation(rotation)
    if not isinstance(ell, (int, np.integer)) or ell not in (0, 1, 2):
        raise ValueError('ordinary original s/p/d angular momentum required')
    if ell == 0:
        return np.ones((1, 1))
    if ell == 1:
        a = np.array([[0, 0, 1], [-1, 0, 0], [0, -1, 0]])
        return (a @ r.T @ a.T).astype(float)
    # Frobenius-orthogonal traceless tensors, each squared norm 3/2.
    root3 = np.sqrt(3.)
    t = np.zeros((5, 3, 3))
    t[0] = np.diag([-.5, -.5, 1.])
    t[1, 0, 2] = t[1, 2, 0] = -root3/2
    t[2, 1, 2] = t[2, 2, 1] = -root3/2
    t[3] = np.diag([root3/2, -root3/2, 0.])
    t[4, 0, 1] = t[4, 1, 0] = root3/2
    rotated = r.T @ t @ r
    return np.einsum('iab,jab->ij', t, rotated) * (2./3.)


def original_channel_intertwiner(upf, rotation):
    """Return Q(nh,nh) in the unchanged radial-then-angular UPF order.

    Repeated l channels get separate blocks indexed by the original radial ID.
    No ordering by l, radial norm, or eigenvalue is allowed. A malformed channel
    ledger is rejected rather than reordered. No active D is read or modified.
    """
    angular = np.asarray(upf.angular_l)
    if (angular.ndim != 1 or not angular.size or angular.dtype.kind not in 'iu'
            or not np.isin(angular, (0, 1, 2)).all()):
        raise ValueError('nonempty original integer s/p/d radial ledger required')
    radial = np.repeat(np.arange(angular.size), 2*angular+1)
    lm = np.concatenate([np.arange(l*l, (l+1)**2) for l in angular])
    if (not np.array_equal(upf.channel_radial, radial)
            or not np.array_equal(upf.channel_l, angular[radial])
            or not np.array_equal(upf.channel_lm, lm) or upf.nh != radial.size):
        raise ValueError('channel ledger is not the original QE radial-then-angular order')
    blocks = {int(l): harmonic_intertwiner(rotation, int(l)) for l in np.unique(angular)}
    result = np.zeros((radial.size, radial.size))
    start = 0
    for ell in angular:
        stop = start+2*int(ell)+1
        result[start:stop, start:stop] = blocks[int(ell)]
        start = stop
    return result


def intertwiner_diagnostics(intertwiner, source_d, target_d=None):
    """Measure Q Dt Q† - Ds, preserving arbitrary finite complex active D.

    A nonzero defect is reported, not repaired or treated as a failing angular
    construction. The existing consumer includes it in the nonlocal bound.
    """
    q = np.asarray(intertwiner)
    ds = np.asarray(source_d)
    dt = ds if target_d is None else np.asarray(target_d)
    if (q.ndim != 2 or q.shape[0] != q.shape[1] or ds.shape != q.shape
            or dt.shape != q.shape or not all(np.isfinite(v).all() for v in (q, ds, dt))):
        raise ValueError('finite square original-channel Q and active D matrices required')
    defect = q @ dt @ q.conj().T-ds
    return dict(
        orthogonality_frobenius=float(np.linalg.norm(q.conj().T @ q-np.eye(len(q)))),
        d_conjugacy_frobenius=float(np.linalg.norm(defect)),
        d_conjugacy_spectral=float(np.linalg.norm(defect, 2)),
        d_source_norm_frobenius=float(np.linalg.norm(ds)),
        d_target_norm_frobenius=float(np.linalg.norm(dt)),
        d_source_hermiticity_frobenius=float(np.linalg.norm(ds-ds.conj().T)),
        d_target_hermiticity_frobenius=float(np.linalg.norm(dt-dt.conj().T)),
        policy='original active D unchanged; nonzero conjugacy defect remains in complete bound')
