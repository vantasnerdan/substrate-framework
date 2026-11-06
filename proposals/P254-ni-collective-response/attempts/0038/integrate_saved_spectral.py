"""Integrate new force channels in authenticated, already-paid spectral coordinates.

No state, Hamiltonian action, eigensolve, inverse or conditioning computation is
performed. Complex gaps and independent supplied adjoint forces remain literal.
"""
from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3] / 'src'), str(HERE.parent / '0031')]

import numpy as np
from filtered_quadratic_bath import _fourier_legendre_moments
from substrate_framework.streamed_material_state import authenticate

SPECTRAL_KEYS = ('eigenvalues', 'right_eigenvectors', 'inverse_right_eigenvectors',
                 'energy_gaps', 'chebyshev_to_legendre')


def load_saved_spectral(binding):
    """Load actual canonical A37 NPZ members, not coordinates rebuilt from H.

    The same public loader serves the serialized small oracle fixture. Other
    members of the authenticated original archive are deliberately not loaded.
    """
    with np.load(authenticate(binding), allow_pickle=False) as archive:
        saved = {key: archive['spectral_' + key] for key in SPECTRAL_KEYS}
    energies = saved['eigenvalues']
    size = energies.size
    if energies.shape != (size,) or not size or energies.dtype != np.complex128:
        raise ValueError('saved eigenvalues must be nonempty complex128 vector')
    for key in SPECTRAL_KEYS[1:4]:
        if saved[key].shape != (size, size) or saved[key].dtype != np.complex128:
            raise ValueError('saved spectral matrix shape/dtype differs: ' + key)
    conversion = saved['chebyshev_to_legendre']
    if conversion.shape != (8, 8) or conversion.dtype != np.float64:
        raise ValueError('require actual degree-seven float64 conversion matrix')
    if not all(np.all(np.isfinite(value)) for value in saved.values()):
        raise ValueError('nonfinite saved spectral coordinates')
    gaps = saved['energy_gaps']
    if (not np.array_equal(gaps, energies[:, None] - energies[None, :])
            or not np.array_equal(-gaps, gaps.T)):
        raise ValueError('saved gaps must be the literal ordered energy differences')
    for value in saved.values():
        value.flags.writeable = False
    return saved


def adjoint_moments_from_forward(moments):
    """For real intervals/P_n: M(g*)=conj(M(-g)), and -g=g.T.

    This is NOT conjugation alone and does not require real energy gaps. The
    adjoint eigenbasis force and coefficient conjugation are separate operations.
    """
    return moments.swapaxes(-1, -2).conj()


def integrate_saved_spectral(saved, forces, force_adjoints, segment_bounds, coefficients):
    """A31 loop algebra with saved V/Vinv/gaps/conversion and one moment pass.

    F/Fadj[mode,b,b], coefficients[segment,degree,mode,probe]. Independent
    Fadj evolves with raw H†, not an assumed Hermitian or manufactured F† leg.
    """
    f = np.asarray(forces, dtype=np.complex128)
    fa = np.asarray(force_adjoints, dtype=np.complex128)
    coefficients = np.asarray(coefficients, dtype=np.complex128)
    bounds = np.asarray(segment_bounds)
    if np.iscomplexobj(bounds):
        raise ValueError('real interval bounds required')
    intervals = np.asarray(bounds, dtype=float)
    size = saved['eigenvalues'].size
    if f.ndim != 3 or f.shape[1:] != (size, size) or fa.shape != f.shape:
        raise ValueError('forces and supplied adjoints require [mode,basis,basis]')
    if coefficients.ndim != 4 or coefficients.shape[2] != len(f):
        raise ValueError('coefficients require [segment,degree,mode,probe]')
    segments, count, modes, probes = coefficients.shape
    if not segments or not modes or not probes or not 1 <= count <= 8:
        raise ValueError('nonempty axes and polynomial degree at most seven required')
    if intervals.shape == (segments + 1,):
        intervals = np.column_stack((intervals[:-1], intervals[1:]))
    if intervals.shape != (segments, 2) or np.any(intervals[:, 1] <= intervals[:, 0]):
        raise ValueError('positive-width intervals[segment,2] required')
    if not all(np.all(np.isfinite(a)) for a in (f, fa, coefficients, intervals)):
        raise ValueError('all force, coefficient and interval inputs must be finite')
    vectors = saved['right_eigenvectors']
    inverse = saved['inverse_right_eigenvectors']
    eigen_forces = inverse @ f @ vectors
    eigen_adjoints = vectors.conj().T @ fa @ inverse.conj().T
    conversion = saved['chebyshev_to_legendre'][:count, :count]
    forward = np.zeros((probes, size, size), dtype=np.complex128)
    adjoint = np.zeros_like(forward)
    for segment, (left, right) in enumerate(intervals):
        half = (right - left) / 2.0
        legendre = np.einsum('dk,kmp->dmp', conversion, coefficients[segment])
        moments = _fourier_legendre_moments(saved['energy_gaps'], half, left + half, count)
        adjoint_moments = adjoint_moments_from_forward(moments)
        for mode in range(modes):
            weights = np.einsum('dij,dp->pij', moments, legendre[:, mode, :])
            forward += weights * eigen_forces[mode]
            weights = np.einsum('dij,dp->pij', adjoint_moments, legendre[:, mode, :].conj())
            adjoint += weights * eigen_adjoints[mode]
    operators = vectors @ forward @ inverse
    adjoint_operators = inverse.conj().T @ adjoint @ vectors.conj().T
    return dict(operators=operators, adjoint_operators=adjoint_operators,
        metadata=dict(basis_count=size, segment_count=segments, polynomial_degree=count - 1,
            mode_count=modes, probe_count=probes, intervals=intervals.tolist(),
            spectral_policy='authenticated original actual V/Vinv/eigenvalues/gaps/conversion; no regeneration',
            hamiltonian_policy='literal raw H; no Hermitization or imaginary-gap clipping',
            adjoint_evolution='raw H† with independent supplied adjoints and conjugated coefficients',
            moment_identity='M_n(g*)=conj(M_n(-g)); -g=g.T; transpose then conjugate',
            numerical_scope='floating-point saved-spectral evaluation, not enclosure or basis error bound'),
        diagnostics=dict(all_outputs_finite=bool(np.all(np.isfinite(operators))
                and np.all(np.isfinite(adjoint_operators))),
            raw_supplied_adjoint_defect_fro_by_probe=np.linalg.norm(
                adjoint_operators - operators.conj().swapaxes(-1, -2), axis=(1, 2)).tolist(),
            moment_evaluations=segments, adjoint_moment_evaluations=0,
            eigensolves=0, inverse_recomputations=0, condition_recomputations=0))
