"""Full, unit-separated representation of the existing causal force filter.

G_ij = sum_p integral K_ip(t) conjugate(K_jp(t)) dt, separately for ionic
and field-momentum forces. Every channel survives: K=U a, a=U^dagger K.
The bilinear output coupling is K.T F = a.T (U.T F), so a physical modal
operator uses U[:,mode], NOT its conjugate. Complex modes need independent
Hilbert-adjoint actions even when the original Cartesian forces are Hermitian.

This is a change of representation, not a passive kernel, a physical loss,
a discarded-noise estimate, or a certificate of the filter's physical error.
Importing performs no packet access, material actions or numerical solve.
"""
from __future__ import annotations

import numpy as np


def _max_abs(value):
    return float(np.max(np.abs(value), initial=0.0))


def _kind_result(kind, coefficients, boundary, gram, canonical_filter):
    channels = coefficients.shape[2]
    # eigh uses the lower triangle; retain the ORIGINAL Gram and its full
    # anti-Hermitian defect rather than averaging, clipping or repairing it.
    if channels:
        eigenvalues, basis = np.linalg.eigh(gram, UPLO='L')
        eigenvalues, basis = eigenvalues[::-1], basis[:, ::-1]
    else:
        eigenvalues = np.empty(0, dtype=float)
        basis = np.empty((0, 0), dtype=complex)
    modal_chebyshev = np.einsum('ia,skip->skap', basis.conj(), coefficients)
    modal_boundary = np.einsum('ia,bip->bap', basis.conj(), boundary)
    defects = {
        'gram_hermiticity_defect': gram-gram.conj().T,
        'basis_unitarity_defect': basis.conj().T@basis-np.eye(channels),
        'eigenpair_defect': gram@basis-basis*eigenvalues,
        'chebyshev_reconstruction_defect': (
            np.einsum('ia,skap->skip', basis, modal_chebyshev)-coefficients),
        'boundary_reconstruction_defect': (
            np.einsum('ia,bap->bip', basis, modal_boundary)-boundary),
    }
    return {
        'basis': basis,
        'gram': gram,
        'eigenvalues': eigenvalues,
        'modal_chebyshev': modal_chebyshev,
        'modal_boundary': modal_boundary,
        'metadata': {
            'kind': kind,
            'force_insertion': ('plus physical ionic force into p-prime' if kind == 'ionic'
                                else 'plus physical field-momentum force into Pi-prime'),
            'channels': channels,
            'modes': channels,
            'nprobes': coefficients.shape[3],
            'empty': channels == 0,
            'q_index': canonical_filter.metadata.get('q_index'),
            'channel_basis': 'original canonical force-coordinate order; no mass/source/cell factors',
            'basis_axes': ['channel', 'mode'],
            'gram_axes': ['channel', 'channel'],
            'modal_chebyshev_axes': ['segment', 'degree', 'mode', 'probe'],
            'modal_boundary_axes': ['boundary', 'mode', 'probe'],
            'modal_coefficients': 'a=U-dagger K; K=U a',
            'physical_modal_force': 'F_mode=sum_i U[i,mode] F_i (no conjugation of U)',
            'eigensolver': 'numpy.linalg.eigh, lower triangle of raw Gram; descending eigenvalues',
            'rank_truncation': False,
            'clipping': False,
            'physical_error_bound': None,
        },
        'diagnostics': {
            **defects,
            'max_abs': {name: _max_abs(value) for name, value in defects.items()},
        },
    }


def decompose_force_channels(canonical_filter):
    """Return full ionic and field channel dictionaries plus shared time axes.

Consumes the actual CanonicalForceFilter arrays AND scalar-time evaluator.
Eight Gauss-Legendre nodes on each accepted physical-time interval integrate
all probe products (degree <=14) without changing the filter or its knots.
Boundary values are transformed independently and are not added as discrete
weights to the continuous Gram. Ionic and field Grams never share an axis.

Each kind contains basis[channel,mode], raw gram[channel,channel], descending
unclipped eigenvalues[mode], modal_chebyshev[segment,degree,mode,probe],
modal_boundary[boundary,mode,probe], metadata and raw diagnostic arrays.
The top-level segment_bounds and boundary_times are copies of the canonical
axes. A zero-channel field has an explicit complete zero-dimensional basis;
a zero-time filter has a zero Gram but still retains every channel and knot.
There is no truncation or fitted numerical acceptance tolerance.
"""
    arrays = canonical_filter.arrays
    bounds = np.asarray(arrays['segment_bounds'])
    times = np.asarray(arrays['boundary_times'])
    if bounds.ndim != 2 or bounds.shape[1] != 2:
        raise ValueError('canonical segment bounds must have shape [segment,2]')
    count = len(bounds)
    if (times.shape != (count+1,) or not np.isfinite(bounds).all()
            or not np.isfinite(times).all() or times[0] != 0.0
            or times[-1] != canonical_filter.final_time):
        raise ValueError('canonical boundaries must span zero through final_time')
    if count and (np.any(np.diff(times) <= 0.0)
                  or not np.array_equal(bounds[:, 0], times[:-1])
                  or not np.array_equal(bounds[:, 1], times[1:])):
        raise ValueError('canonical accepted intervals must be contiguous and ascending')
    coefficients, boundaries, grams = {}, {}, {}
    nprobes = None
    for kind in ('ionic', 'field'):
        coefficient = np.asarray(arrays[kind+'_chebyshev'])
        boundary = np.asarray(arrays[kind+'_boundary'])
        if (coefficient.ndim != 4 or coefficient.shape[:2] != (count, 8)
                or boundary.shape != (count+1, *coefficient.shape[2:])
                or not np.isfinite(coefficient).all() or not np.isfinite(boundary).all()):
            raise ValueError(kind+' must retain canonical degree-seven coefficient and boundary axes')
        if nprobes is not None and coefficient.shape[3] != nprobes:
            raise ValueError('ionic and field terminal probe counts differ')
        nprobes = coefficient.shape[3]
        coefficients[kind], boundaries[kind] = coefficient, boundary
        channels = coefficient.shape[2]
        grams[kind] = np.zeros((channels, channels), dtype=complex)
    nodes, weights = np.polynomial.legendre.leggauss(8)
    for left, right in bounds:
        half = (right-left)/2.0
        for node, weight in zip(nodes, weights):
            time = left+half*(node+1.0)
            ionic, field = canonical_filter.evaluate(time)
            for kind, value in (('ionic', ionic), ('field', field)):
                value = np.asarray(value)
                if value.shape != coefficients[kind].shape[2:] or not np.isfinite(value).all():
                    raise ValueError('canonical evaluator changed '+kind+' force/probe axes')
                grams[kind] += (half*weight)*(value@value.conj().T)
    return {
        **{kind: _kind_result(kind, coefficients[kind], boundaries[kind], grams[kind],
                              canonical_filter) for kind in ('ionic', 'field')},
        'segment_bounds': bounds.copy(),
        'boundary_times': times.copy(),
        'metadata': {
            'schema': 'P254-full-force-channel-representation-v1',
            'canonical_filter': dict(canonical_filter.metadata),
            'gram': 'G_ij=sum_probe integral K_i,probe(t) conjugate(K_j,probe(t)) dt',
            'quadrature': '8-point Gauss-Legendre per accepted interval; degree14 products',
            'time_measure': 'physical dt in inverse Ry; (right-left)/2 Jacobian once',
            'force_kinds': 'separate ionic and field units; no cross-kind Gram or rotation',
            'boundary_policy': 'transform exact stored knots independently; no point masses in Gram',
            'extra_mass_source_cell_factors': False,
            'all_modes_retained': True,
            'scope': 'full output-directed representation only; not passive kernel completion or physical loss',
        },
    }
