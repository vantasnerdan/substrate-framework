"""Bounded, full-support local projection with SciPy complex batch FFTs.

Integration (after authenticating/archiving this NEW producer dependency)::

    from reduced_local_projection import project_edge
    arrays, timings = project_edge(coupling, source_columns, target_columns, gs)

For the outer-q/inner-slot loop, retain ONE preparation while the existing local
cache retains its support/gradient, then pass ``prepared=preparation``. Create a
new preparation when that support changes::

    preparation = prepare_local_projection(coupling)
    arrays, timings = project_edge(coupling, source, target, gs,
                                  prepared=preparation)

``project_local_pairs`` supplies just the local/density/potential part;
``add_nonlocal_projection`` adds the unchanged small projector contractions.
No source-star, opposite-q inference, occupation, or spin weight is introduced.

Derivation: let N=prod(shape), Omega=cell volume, ell=reciprocal_fold and
p(r)=exp(-2 pi i ell.r). The directed dyad is
n_ab(r)=u_ta(r) conj(u_sb(r)) conj(p(r))/Omega. Both endpoint inverse FFTs
reconstruct amplitudes (norm='forward'); the pair FFT has the 1/N forward
normalization. Thus density_G=sqrt(Omega)*FFT(n)_G and the profile-action
coefficient is conj(density_G), with NO further N/Omega/spin factors. Physical
local coefficients are independently evaluated as M.T @ FFT(n) and
M.H @ conj(FFT(n)), where M is the existing full-G gradient including Omega.
All flattening uses the physical F-order convention. Internally C-order
[pair,z,y,x] is exactly the same grid, not a changed Miller/axis convention.

If physical G and G' alias on the discrete grid, FFT(n)_G=FFT(n)_G'. Only THEN
may their gradient rows be summed: A[slot,d]=sum_(G->slot) M[G,d]. This is an
exact discrete Fourier identity, retaining every radial/position/q+G term,
including Nyquist and opposite-edge extensions. It is not a support cutoff.
Forward and reverse still use two separate contractions; no matrix is repaired.

Cost/limitations: this deliberately retains m*m FULL-GRID pair transforms.
It reduces copies/Python assembly overhead, uses bounded in-place SciPy batches,
and reduces both GEMM row counts only where genuine grid aliases exist. There
is no measured speedup claim. Without aliases the arithmetic remains
O(m*m*N log N + m*m*NG*ndof). An all-coordinate real-space field contraction
would replace NG by N in the expensive dense contraction and require ndof
field transforms/storage plus the density features; it is NOT generically a
cheaper route. Arbitrary atomic positions do not permit an exact ordinary FFT
over atoms. Feature-only sparse projection cannot replace the full ionic work.

The scratch budget bounds explicit pair/sampling buffers conservatively (two
pair grids plus one sampled grid); it is NOT an RSS or FFT-library-scratch cap.
Inputs, two m-column wave grids, outputs, existing M and any reusable aliased
M are additional. Pair batches never exceed ``batch_size``. Independent local
roundoff changes from alias summation and pocketfft/BLAS ordering are expected;
no threshold, clipping, real-part projection or Hermitian repair is performed.
"""
from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy import fft
from scipy.linalg.blas import zgemm

from reduced_local_gradient import local_gradient_matrix
from reduced_projection_algebra import projector_overlaps


@dataclass(frozen=True, eq=False)
class PreparedLocalProjection:
    """One immutable q/support binding; do not retain a full bank of these.

    ``support`` and ``original_gradient`` must remain immutable. Identity checks
    reject reuse across a replacement support/gradient even of equal shape.
    """
    support: object
    original_gradient: np.ndarray
    slots: np.ndarray
    gradient: np.ndarray
    preparation_seconds: float
    gradient_build_seconds: float
    alias_count: int


def prepare_local_projection(ionic):
    """Reuse the physical M, additively fold ONLY identical discrete FFT slots."""
    started = perf_counter()
    support = ionic.bridge.density_fft
    original = getattr(ionic, 'local_gradient', None)
    build_seconds = 0.
    if original is None:
        t = perf_counter()
        original, _ = local_gradient_matrix(ionic.inputs, support, ionic.vloc)
        build_seconds = perf_counter()-t
    original = np.asarray(original)
    nl = np.asarray(support.nl)
    if (original.shape != (len(nl), ionic.ndof) or original.dtype != np.complex128
            or nl.ndim != 1 or nl.dtype.kind not in 'iu'
            or np.any(nl < 0) or np.any(nl >= support.nnr)):
        raise ValueError('complex128 full physical gradient and valid F-order slots required')
    unique, inverse = np.unique(nl, return_inverse=True)
    aliases = len(nl)-len(unique)
    if aliases:
        matrix = np.zeros((len(unique), ionic.ndof), dtype=np.complex128, order='F')
        np.add.at(matrix, inverse, original)
        slots = unique
        matrix.setflags(write=False)
    else:
        # Keep the original row order: no need to copy/reorder a large matrix.
        matrix = original
        slots = np.array(nl, copy=True)
    slots.setflags(write=False)
    return PreparedLocalProjection(support, original, slots, matrix,
                                   perf_counter()-started, build_seconds, aliases)


def _waves(mapping, columns, workers):
    """[band,z,y,x] complex inverse FFT; additive placement is alias-safe."""
    nband = columns.shape[1]
    padded = np.zeros((nband, mapping.nnr), dtype=np.complex128)
    # Unlike assignment, this also correctly handles additive endpoint aliases.
    np.add.at(padded.T, mapping.nl, columns)
    grid = padded.reshape((nband, *mapping.physical[::-1]))
    return fft.ifftn(grid, axes=(-3, -2, -1), norm='forward', workers=workers,
                     overwrite_x=True).reshape((nband, mapping.nnr))


def _inputs(ionic, source, target, gs, prepared):
    bridge = ionic.bridge
    shape = tuple(bridge.source_fft.physical)
    source, target, gs = np.asarray(source), np.asarray(target), np.asarray(gs)
    if (len(shape) != 3 or tuple(bridge.target_fft.physical) != shape
            or tuple(bridge.density_fft.physical) != shape
            or bridge.target_fft.omega != bridge.source_fft.omega
            or bridge.density_fft.omega != bridge.source_fft.omega):
        raise ValueError('endpoint and density physical grids/volumes must agree')
    if (source.ndim != 2 or target.ndim != 2 or not source.shape[1]
            or target.shape[1] != source.shape[1]
            or source.shape[0] != len(bridge.source_fft.nl)
            or target.shape[0] != len(bridge.target_fft.nl)):
        raise ValueError('equal nonempty active band counts and full endpoint rows required')
    if gs.ndim != 2 or gs.shape[1] != 3 or gs.dtype.kind not in 'iu':
        raise ValueError('integer Fourier feature rows [nfeature,3] required')
    if ionic.edge.ispin not in (1, 2):
        raise ValueError('explicit physical spin 1 or 2 required')
    if prepared.support is not bridge.density_fft:
        raise ValueError('prepared local projection belongs to a different physical support')
    current_gradient = getattr(ionic, 'local_gradient', None)
    if current_gradient is not None and current_gradient is not prepared.original_gradient:
        raise ValueError('prepared local projection belongs to a different physical gradient')
    if np.shape(bridge.multiplication_phase) != (bridge.source_fft.nnr,):
        raise ValueError('F-order reciprocal-fold multiplication phase required')
    return source, target, gs, shape


def project_local_pairs(ionic, source, target, gs, *, prepared=None,
                        batch_size=8, workspace_mib=128., workers=1):
    """Return the complete local part in the existing project_edge layout.

    ``prepared`` is optional and can be shared by slots with the SAME immutable
    support and M objects. Workers default to one, avoiding hidden thread fanout.
    This function makes no native/material calls. It does not mutate inputs.
    """
    started = perf_counter()
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError('positive integer pair batch_size required')
    if type(workers) is not int or workers < 1:
        raise ValueError('positive explicit FFT worker count required')
    if not np.isfinite(workspace_mib) or workspace_mib <= 0:
        raise ValueError('positive finite explicit pair workspace MiB required')
    preparation_start = perf_counter()
    prepared = prepare_local_projection(ionic) if prepared is None else prepared
    preparation_seconds = perf_counter()-preparation_start
    source, target, gs, shape = _inputs(ionic, source, target, gs, prepared)
    bridge = ionic.bridge
    omega, nnr = bridge.source_fft.omega, bridge.source_fft.nnr
    m, ng, spin = source.shape[1], len(gs), ionic.edge.ispin-1
    # Two grids account for a library choosing not to overwrite its input.
    per_pair_bytes = np.dtype(np.complex128).itemsize*(2*nnr+len(prepared.slots))
    capacity = int(workspace_mib*1024**2)//per_pair_bytes
    if capacity < 1:
        raise MemoryError('pair workspace cannot hold one full-support pair transform')
    batch = min(batch_size, capacity, m*m)
    wave_start = perf_counter()
    source_wave = _waves(bridge.source_fft, source, workers)
    target_wave = _waves(bridge.target_fft, target, workers)
    np.conjugate(source_wave, out=source_wave)
    phase = np.conjugate(bridge.multiplication_phase)/omega
    wave_seconds = perf_counter()-wave_start
    feature_slots = np.ravel_multi_index(tuple(gs[:, a] % shape[a] for a in range(3)),
                                        shape, order='F')
    density = np.zeros((2*ng, m, m), dtype=np.complex128)
    forward = np.empty((ionic.ndof, m, m), dtype=np.complex128)
    reverse = np.empty_like(forward)
    density_flat = density.reshape((2*ng, m*m))
    forward_flat = forward.reshape((ionic.ndof, m*m))
    reverse_flat = reverse.reshape((ionic.ndof, m*m))
    pair = np.empty((batch, nnr), dtype=np.complex128)
    sampled = np.empty((len(prepared.slots), batch), dtype=np.complex128, order='F')
    fft_seconds = gemm_seconds = 0.
    pair_start = perf_counter()
    for start in range(0, m*m, batch):
        stop = min(start+batch, m*m)
        count = stop-start
        fft_start = perf_counter()
        # Assemble contiguous source-band runs without two large advanced-index
        # copies of the wave grids and without a Python loop per grid point.
        offset = 0
        while offset < count:
            a, b = divmod(start+offset, m)
            length = min(m-b, count-offset)
            np.multiply(target_wave[a:a+1], source_wave[b:b+length],
                        out=pair[offset:offset+length])
            offset += length
        pair[:count] *= phase
        spectrum = fft.fftn(pair[:count].reshape((count, *shape[::-1])),
                            axes=(-3, -2, -1), norm='forward',
                            workers=workers, overwrite_x=True).reshape((count, nnr))
        current = sampled[:, :count]
        # Validated indices allow mode='clip', avoiding np.take's mandatory
        # mode='raise' output-buffer protection copy.
        np.take(spectrum, prepared.slots, axis=1, out=current.T, mode='clip')
        features = np.take(spectrum, feature_slots, axis=1).T
        density_flat[spin*ng:(spin+1)*ng, start:stop] = np.sqrt(omega)*features
        fft_seconds += perf_counter()-fft_start
        gemm_start = perf_counter()
        # Independent physical moment and action, NOT conj(reverse) assignment.
        reverse_flat[:, start:stop] = zgemm(1., prepared.gradient, current, trans_a=1)
        np.conjugate(current, out=current)
        forward_flat[:, start:stop] = zgemm(1., prepared.gradient, current, trans_a=2)
        gemm_seconds += perf_counter()-gemm_start
    arrays = dict(ionic=forward, ionic_reverse=reverse, density=density,
                  potential=density.conj().transpose(1, 2, 0))
    timings = dict(wave_fft_seconds=wave_seconds,
                   pair_fft_all_local_seconds=perf_counter()-pair_start,
                   pair_density_fft_seconds=fft_seconds,
                   local_gradient_gemm_seconds=gemm_seconds,
                   local_gradient_prepare_seconds=preparation_seconds,
                   local_gradient_build_seconds=prepared.gradient_build_seconds,
                   prepared_local_build_seconds=prepared.preparation_seconds,
                   pair_batch_size=batch, pair_ffts=m*m, coordinate_full_actions=0,
                   local_density_support=len(prepared.support.nl),
                   local_contraction_rows=len(prepared.slots),
                   local_additive_aliases=prepared.alias_count,
                   all_coordinates=ionic.ndof, fft_workers=workers,
                   pair_explicit_scratch_bound_bytes=per_pair_bytes*batch,
                   endpoint_wave_bytes=source_wave.nbytes+target_wave.nbytes,
                   folded_gradient_extra_bytes=(prepared.gradient.nbytes if prepared.alias_count else 0),
                   algorithm='scipy_pair_fft_full_support_additive_alias_fold',
                   total_seconds=perf_counter()-started)
    return arrays, timings


def add_nonlocal_projection(ionic, source, target, forward, reverse):
    """Add existing independently derived NC moments, reusing cached overlaps."""
    if hasattr(ionic, 'source_overlaps'):
        so, to = ionic.source_overlaps, ionic.target_overlaps
    else:
        so = projector_overlaps(ionic.source_factors, source, ionic.ps)
        to = projector_overlaps(ionic.target_factors, target, ionic.pt)
    for ss, ts in ionic.species:
        sb, sp = so[ss.name]
        tb, tp = to[ts.name]
        d = ss.active_d
        for atom_slot, atom in enumerate(ss.atom_indices):
            for axis in range(3):
                forward[3*(atom-1)+axis] += 1j*(
                    tb[atom_slot].conj().T @ d @ sp[axis, atom_slot]
                    -tp[axis, atom_slot].conj().T @ d @ sb[atom_slot])
                reverse[3*(atom-1)+axis] += 1j*(
                    (sb[atom_slot].conj().T @ d.conj().T @ tp[axis, atom_slot]).T
                    -(sp[axis, atom_slot].conj().T @ d.conj().T @ tb[atom_slot]).T)


def project_edge(ionic, source, target, gs, *, prepared=None, batch_size=8,
                 workspace_mib=128., workers=1):
    """Drop-in full local + nonlocal edge export; existing positional API."""
    started = perf_counter()
    arrays, timings = project_local_pairs(ionic, source, target, gs,
        prepared=prepared, batch_size=batch_size, workspace_mib=workspace_mib, workers=workers)
    nonlocal_start = perf_counter()
    add_nonlocal_projection(ionic, source, target, arrays['ionic'], arrays['ionic_reverse'])
    timings['overlap_nonlocal_seconds'] = perf_counter()-nonlocal_start
    timings['total_seconds'] = perf_counter()-started
    return arrays, timings
