"""Integrate streamed physical density residuals with one shared moment pass.

The finite bath is the paid raw frozen-H bath, not a moving full-PW mean.
All full/retained matching occurs at the operator level before Wick contraction.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / name) for name in ('0038', '0031')]
import numpy as np
from filtered_quadratic_bath import _fourier_legendre_moments, _chebyshev_to_legendre
from integrate_saved_spectral import adjoint_moments_from_forward


def integrate_density_residual(saved, feature_provider, feature_indices, projector,
                               bounds, coefficients, *, target_saved=None, target_projector=None):
    """Cache only selected residual eigenvertices; never persist the mode bank.

    feature_provider(r)->(physical O_r, independently applied O_r†), already
    normalized by 1/sqrt(N). Exactly zero other-spin channels are not allocated.
    Each physical feature is projected ONCE, each interval/gap moment ONCE.
    Resident cache: 2*m*b*b*16 bytes; moments: 2*15*b*b*16; outputs:2*p*b*b*16.
    The caller's profile generator uses bounded PW-column batches.
    For nonzero q supply target_saved and target_projector from the actual
    k+q bank; vertices are [source_basis,target_basis], adjoints reversed.
    """
    coefficients = np.asarray(coefficients, complex)
    bounds = np.asarray(bounds, float)
    indices = tuple(int(i) for i in feature_indices)
    if coefficients.ndim != 4 or coefficients.shape[1] != 15:
        raise ValueError('density coefficients require [segment,15,feature,probe]')
    if bounds.shape != (len(coefficients), 2) or np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError('positive density intervals required')
    v, vi = saved['right_eigenvectors'], saved['inverse_right_eigenvectors']
    target_saved = saved if target_saved is None else target_saved
    target_projector = projector if target_projector is None else target_projector
    tv, tvi = target_saved['right_eigenvectors'], target_saved['inverse_right_eigenvectors']
    b, bt, probes = len(v), len(tv), coefficients.shape[-1]
    gaps = (saved['energy_gaps'] if target_saved is saved else
            saved['eigenvalues'][:, None] - target_saved['eigenvalues'][None, :])
    eigen = np.empty((len(indices), b, bt), complex)
    adj = np.empty((len(indices), bt, b), complex)
    defects = []
    for mode, feature in enumerate(indices):
        full, full_adjoint = feature_provider(feature)
        residual = full - projector @ full @ target_projector
        residual_adjoint = full_adjoint - target_projector @ full_adjoint @ projector
        eigen[mode] = vi @ residual @ tv
        adj[mode] = tv.conj().T @ residual_adjoint @ vi.conj().T
        defects.append(float(np.linalg.norm(full_adjoint - full.conj().T)))
    forward = np.zeros((probes, b, bt), complex)
    backward = np.zeros((probes, bt, b), complex)
    conversion = _chebyshev_to_legendre(15)
    if not np.array_equal(conversion[:8, :8], saved['chebyshev_to_legendre']):
        raise ValueError('higher-order conversion must preserve paid degree7 block')
    for segment, (left, right) in enumerate(bounds):
        half = (right-left)/2
        legendre = np.einsum('dk,kmp->dmp', conversion,
                             np.take(coefficients[segment], indices, axis=1))
        moments = _fourier_legendre_moments(gaps, half, left+half, 15)
        adjoint_moments = adjoint_moments_from_forward(moments)
        for mode in range(len(indices)):
            active = np.flatnonzero(np.any(legendre[:, mode] != 0, axis=0))
            for probe in active:
                forward[probe] += np.einsum('dij,d->ij', moments,
                    legendre[:, mode, probe]) * eigen[mode]
                backward[probe] += np.einsum('dij,d->ij', adjoint_moments,
                    legendre[:, mode, probe].conj()) * adj[mode]
    operators = v @ forward @ tvi
    adjoints = tvi.conj().T @ backward @ v.conj().T
    return dict(operators=operators, adjoint_operators=adjoints,
        diagnostics=dict(all_outputs_finite=bool(np.isfinite(operators).all()
            and np.isfinite(adjoints).all()), moment_evaluations=len(bounds),
            physical_features=len(indices), supplied_adjoint_defects=defects,
            resident_eigenvertex_bytes=int(eigen.nbytes+adj.nbytes), eigensolves=0),
        matching='integral density_dual_r(t)*(O_full_r(t)-P O_full_r(t) P); retained physical vertex propagated by SAME paid H')
