"""Integrate streamed physical density residuals with one shared moment pass.

The finite bath is the paid raw frozen-H bath, not a moving full-PW mean.
All full/retained matching occurs at the operator level before Wick contraction.
"""
from pathlib import Path
import sys
from time import perf_counter
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / name) for name in ('0038', '0031')]
import numpy as np
from filtered_quadratic_bath import _chebyshev_to_legendre
from density_moments import fourier_legendre_moments


def integrate_density_residual(saved, feature_provider, feature_indices, projector,
                               bounds, coefficients, *, target_saved=None, target_projector=None,
                               progress=None, checkpoint=None, row_block=128):
    """Integrate scalar profiles per literal spectral row block, then apply vertices.

    A completed immutable block contains the forward rows and independent
    reversed-adjoint columns. No all-mode profile bank is persisted. Every
    interval, degree, feature and probe remains included. Only a missing block
    incurs convolution work; complete recovery does not invoke the provider.
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
    windows = [(start, min(start+row_block, b)) for start in range(0, b, row_block)]
    keys = ['spectral-rows/%06d-%06d' % window for window in windows]
    reused = [checkpoint is not None and checkpoint.has(key) for key in keys]
    emit = progress if progress is not None else lambda event: None
    defects = []
    eigen = adj = None
    if not all(reused):
        eigen = np.empty((len(indices), b, bt), complex)
        adj = np.empty((len(indices), bt, b), complex)
        for mode, feature in enumerate(indices):
            full, full_adjoint = feature_provider(feature)
            residual = full - projector @ full @ target_projector
            residual_adjoint = full_adjoint - target_projector @ full_adjoint @ projector
            eigen[mode] = vi @ residual @ tv
            adj[mode] = tv.conj().T @ residual_adjoint @ vi.conj().T
            defects.append(float(np.linalg.norm(full_adjoint-full.conj().T)))
            emit(dict(stage='density_eigenvertex', feature=feature, status='complete'))
    conversion = _chebyshev_to_legendre(15)
    if not np.array_equal(conversion[:8, :8], saved['chebyshev_to_legendre']):
        raise ValueError('higher-order conversion must preserve paid degree7 block')
    legendre = np.einsum('dk,skmp->sdmp', conversion,
                         np.take(coefficients, indices, axis=2))
    weights = np.ascontiguousarray(legendre.transpose(0, 3, 2, 1))
    forward = np.zeros((probes, b, bt), complex)
    backward = np.zeros((probes, bt, b), complex)
    new_blocks = moment_evaluations = 0
    for block, ((start, stop), key) in enumerate(zip(windows, keys)):
        started = perf_counter()
        if not reused[block]:
            emit(dict(stage='density_spectral_rows', block=block, blocks=len(windows),
                      row_start=start, row_stop=stop, status='started'))
            profiles = np.zeros((probes, len(indices), stop-start, bt), complex)
            flat = profiles.reshape(probes*len(indices), -1)
            moment_seconds = 0.
            for segment, (left, right) in enumerate(bounds):
                moment_started = perf_counter()
                moments = fourier_legendre_moments(gaps[start:stop], (right-left)/2,
                                                    (right+left)/2, 15)
                moment_seconds += perf_counter()-moment_started
                flat += weights[segment].reshape(-1, 15)@moments.reshape(15, -1)
                moment_evaluations += 1
            value = dict(
                forward=np.einsum('pmij,mij->pij', profiles, eigen[:, start:stop]),
                backward=np.einsum('pmij,mji->pji', profiles.conj(), adj[:, :, start:stop]))
            if not all(np.isfinite(array).all() for array in value.values()):
                raise FloatingPointError('nonfinite spectral convolution block; no clipping')
            if checkpoint is not None:
                checkpoint.put(key, value, metadata=dict(row_start=start, row_stop=stop,
                    supplied_adjoint_defects=defects),
                    timings=dict(seconds=perf_counter()-started, moment_seconds=moment_seconds))
            new_blocks += 1
            del profiles, flat, moments
        if checkpoint is not None:
            value = checkpoint.get(key)
            if not defects:
                defects = checkpoint.manifest['items'][key]['metadata']['supplied_adjoint_defects']
        if value['forward'].shape != (probes, stop-start, bt) or value['backward'].shape != (probes, bt, stop-start):
            raise ValueError('saved spectral row block axes differ')
        forward[:, start:stop] = value['forward']
        backward[:, :, start:stop] = value['backward']
        emit(dict(stage='density_spectral_rows', block=block, blocks=len(windows),
                  row_start=start, row_stop=stop, status='reused' if reused[block] else 'complete',
                  elapsed_seconds=perf_counter()-started))
    if checkpoint is not None and checkpoint.manifest['status'] != 'complete':
        checkpoint.manifest['status'] = 'complete'
        checkpoint.publish()
    operators = v @ forward @ tvi
    adjoints = tvi.conj().T @ backward @ v.conj().T
    return dict(operators=operators, adjoint_operators=adjoints,
        diagnostics=dict(all_outputs_finite=bool(np.isfinite(operators).all()
            and np.isfinite(adjoints).all()), moment_evaluations=moment_evaluations,
            physical_features=len(indices), supplied_adjoint_defects=defects,
            resident_eigenvertex_bytes=0 if eigen is None else int(eigen.nbytes+adj.nbytes),
            new_spectral_row_blocks=new_blocks, reused_spectral_row_blocks=sum(reused),
            complete_spectral_row_blocks=len(windows), row_block_size=row_block, eigensolves=0),
        matching='integral density_dual_r(t)*(O_full_r(t)-P O_full_r(t) P); retained physical vertex propagated by SAME paid H')
