"""Physical finite-q pair FFT and small nonlocal projector contractions."""
from time import perf_counter
import numpy as np
from scipy.linalg.blas import zgemm
from reduced_local_gradient import local_gradient_matrix


def feature_indices(count=27):
    if count not in (7, 27):
        raise ValueError('7 or 27 Fourier features required')
    first = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    rest = [i for i in np.ndindex(3, 3, 3) if tuple(v-1 for v in i) not in first]
    return np.asarray((first+[tuple(v-1 for v in i) for i in rest])[:count], dtype=int)


def q_feature_geometry(iq, family=27):
    """Nested physical-wavevector closure, including self-inverse Nyquist q."""
    indices = np.indices((4, 4, 4)).reshape(3, -1).T
    index = indices[iq]
    iminus = int(np.ravel_multi_index(tuple(-index % 4), (4, 4, 4)))
    q = index/4.; q -= np.floor(q+.5)
    qm = indices[iminus]/4.; qm -= np.floor(qm+.5)
    shift = np.rint(q+qm).astype(int)
    base = feature_indices(family)
    def list_for(i):
        if iq == iminus:
            items = [tuple(g) for g in base]
            for g in base:
                partner = tuple(-g-shift)
                if partner not in items:
                    items.append(partner)
            return np.asarray(items, int)
        return base if i < (iminus if i == iq else iq) else -base-shift
    gs = list_for(iq)
    opposite = list_for(iminus)
    lookup = {tuple(g): j for j, g in enumerate(opposite)}
    partner = np.asarray([lookup[tuple(-g-shift)] for g in gs], int)
    return gs, np.r_[partner, len(gs)+partner]


def profile(shape, g, omega):
    out = np.ones(shape, complex)
    for axis in range(3):
        dims = [1, 1, 1]
        dims[axis] = shape[axis]
        out *= np.exp(2j*np.pi*g[axis]*np.arange(shape[axis])/shape[axis]).reshape(dims)
    return out.reshape(-1, order='F')/np.sqrt(omega)


def feature_projection(field, shape, omega, gs):
    """Integral conj(phi_G) field, preserving F order and complex values."""
    spectrum = np.fft.fftn(np.asarray(field).reshape((*shape, *np.shape(field)[1:]), order='F'),
                           axes=(0, 1, 2), norm='forward')
    return np.sqrt(omega)*spectrum[tuple(gs[:, a] % shape[a] for a in range(3))]


def projector_overlaps(factors, columns, momenta):
    result = {}
    for species in factors.species:
        shape = (len(species.atom_indices), species.active_d.shape[0], columns.shape[1])
        base = (species.beta.conj().T @ columns).reshape(shape)
        derivative = np.asarray([(species.beta.conj().T @ (momenta[:, axis, None]*columns)).reshape(shape)
                                 for axis in range(3)])
        result[species.name] = (base, derivative)
    return result


def project_edge(ionic, source, target, gs):
    """All full-support local and NC coordinates; only m² pair FFTs per edge.

    reverse[i,a,b] is the independently evaluated positive-gradient coefficient
    Tr(g[e_i]^dagger |t,a><s,b|). Forward is <t,a|g[e_i]|s,b>.
    Potential is the exact physical density-pairing dual, not a symmetry repair.
    """
    t0 = perf_counter()
    bridge = ionic.bridge
    omega, shape = bridge.source_fft.omega, bridge.source_fft.physical
    m, ng, spin = source.shape[1], len(gs), ionic.edge.ispin-1
    density = np.zeros((2*ng, m, m), complex)
    forward = np.zeros((ionic.ndof, m, m), complex)
    reverse = np.zeros_like(forward)
    sw = bridge.source_fft.from_g(source)
    tw = bridge.target_fft.from_g(target)
    wave_seconds = perf_counter()-t0
    gradient_start = perf_counter()
    gradient = getattr(ionic, 'local_gradient', None)
    if gradient is None:
        gradient, _ = local_gradient_matrix(ionic.inputs, bridge.density_fft, ionic.vloc)
    gradient_seconds = perf_counter()-gradient_start
    feature_slots = np.ravel_multi_index(tuple(gs[:, axis] % shape[axis] for axis in range(3)), shape, order='F')
    t1 = perf_counter()
    fft_seconds = gemm_seconds = 0.
    batch_size = 8
    for start in range(0, m*m, batch_size):
        pairs = [(i//m, i % m) for i in range(start, min(start+batch_size, m*m))]
        fft_start = perf_counter()
        # Bounded multi-pair FFT: one transform supplies full local support AND
        # selected Fourier features. No atom loop remains inside the pair loop.
        pair_density = np.empty((bridge.source_fft.nnr, len(pairs)), complex, order='F')
        for column, (a, b) in enumerate(pairs):
            pair_density[:, column] = tw[:, a]*sw[:, b].conj()*bridge.multiplication_phase.conj()/omega
        spectrum = np.fft.fftn(pair_density.reshape((*shape, len(pairs)), order='F'),
                               axes=(0, 1, 2), norm='forward')
        flat = spectrum.reshape((bridge.source_fft.nnr, len(pairs)), order='F')
        sampled = np.asfortranarray(flat[bridge.density_fft.nl])
        features = np.sqrt(omega)*flat[feature_slots]
        fft_seconds += perf_counter()-fft_start
        gemm_start = perf_counter()
        # TWO independent evaluations of the physical transpose/action formulas:
        # reverse=M.T*n and forward=M.H*conj(n). This is not a matrix repair.
        local_reverse = zgemm(1., gradient, sampled, trans_a=1)
        local_forward = zgemm(1., gradient, sampled.conj(), trans_a=2)
        gemm_seconds += perf_counter()-gemm_start
        for column, (a, b) in enumerate(pairs):
            density[spin*ng:(spin+1)*ng, a, b] = features[:, column]
            reverse[:, a, b] = local_reverse[:, column]
            forward[:, a, b] = local_forward[:, column]
    pair_seconds = perf_counter()-t1
    t2 = perf_counter()
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
                forward[3*(atom-1)+axis] += 1j*(tb[atom_slot].conj().T @ d @ sp[axis, atom_slot]
                    -tp[axis, atom_slot].conj().T @ d @ sb[atom_slot])
                reverse[3*(atom-1)+axis] += 1j*(
                    (sb[atom_slot].conj().T @ d.conj().T @ tp[axis, atom_slot]).T
                    -(sp[axis, atom_slot].conj().T @ d.conj().T @ tb[atom_slot]).T)
    potential = density.conj().transpose(1, 2, 0)
    timings = dict(wave_fft_seconds=wave_seconds, pair_fft_all_local_seconds=pair_seconds,
                   pair_density_fft_seconds=fft_seconds, local_gradient_gemm_seconds=gemm_seconds,
                   local_gradient_prepare_seconds=gradient_seconds, pair_batch_size=batch_size,
                   overlap_nonlocal_seconds=perf_counter()-t2, total_seconds=perf_counter()-t0,
                   pair_ffts=m*m, coordinate_full_actions=0,
                   local_density_support=int(len(ionic.gq)), all_coordinates=ionic.ndof)
    return dict(ionic=forward, ionic_reverse=reverse, density=density, potential=potential), timings


def project_kernel(background, q_cart, gs):
    t0 = perf_counter()
    shape, omega = background.fft.physical, background.fft.omega
    r = 2*len(gs)
    out = np.empty((r, r), complex)
    for j in range(r):
        field = np.zeros((background.fft.nnr, 2), complex, order='F')
        field[:, j//len(gs)] = profile(shape, gs[j % len(gs)], omega)
        action = background.apply_kernel(field, q_bohr_inv=q_cart)
        out[:, j] = feature_projection(action, shape, omega, gs).T.reshape(-1)
    return out, dict(seconds=perf_counter()-t0, kernel_actions=r)


def projection_witness(ionic, source, target, arrays, field_vertices, field_matrix, gs):
    """Representative existing physical actions, intentionally not all axes."""
    u = np.zeros((ionic.nat, 3)); u[0, 0] = 1
    action = ionic.prepare(u)
    direct = target.conj().T @ action.apply(source)
    reverse = (source.conj().T @ action.adjoint_apply(target)).T
    f = field_vertices.prepare(np.array([1., 0., 0.]))
    direct_field = source.conj().T @ f.apply(source)
    spin, ng = ionic.edge.ispin-1, len(gs)
    j = min(1, ng-1)
    v = profile(ionic.bridge.source_fft.physical, gs[j], ionic.bridge.source_fft.omega)
    direct_v = target.conj().T @ ionic.bridge.local_apply(v, source)
    defect = lambda x, y: float(np.linalg.norm(x-y)/max(1., np.linalg.norm(x), np.linalg.norm(y)))
    return dict(ionic_axis0_relative=defect(direct, arrays['ionic'][0]),
                ionic_reverse_axis0_relative=defect(reverse, arrays['ionic_reverse'][0]),
                field_axis0_relative=defect(direct_field, field_matrix[0]),
                potential_feature_relative=defect(direct_v, arrays['potential'][:, :, spin*ng+j]),
                witness_actions=dict(ionic_forward=1, ionic_adjoint=1, field=1, potential=1))
