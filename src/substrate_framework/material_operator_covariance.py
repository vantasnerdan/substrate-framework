"""Current-operator covariance bounds, not a future nonlinear symmetry claim.

For declared unitary U, this assesses Ht-U Hs U*. Fourier embeddings use the
normalized physical-grid inner product; no PW-square matrix is constructed.
All upper bounds are analytic inequalities evaluated in ordinary floating
point, NOT outward-rounded interval certificates. Probe residuals are only
empirical witnesses. An antiunitary intertwiner reverses the sign of i and
therefore cannot reuse positive-time propagation by the unitary rule.

A zero exact-real component sum proves covariance. A nonzero sum only supports
an explicitly budgeted current-stage approximation; for Hermitian Hs,Ht the
Duhamel propagator error is at most |time| times this sum. Nothing here licenses
future nonlinear stages or omitted displacement/field operator contributions.
"""
from dataclasses import dataclass
import numpy as np


def _finite(value, name):
    value = np.asarray(value)
    if not np.isfinite(value).all():
        raise ValueError(name + ' must be finite')
    return value


def _permutation(value, size, name):
    value = np.asarray(value)
    if (value.shape != (size,) or value.dtype.kind not in 'iu'
            or not np.array_equal(np.sort(value), np.arange(size))):
        raise ValueError(name + ' must be a complete permutation')
    return value


@dataclass(frozen=True)
class UnitaryPWMap:
    """target[t] = phase[t] * source[row_pullback[t]], no conjugation."""
    row_pullback: np.ndarray
    phase: np.ndarray
    source_identity: str
    target_identity: str
    antiunitary: bool = False

    def __post_init__(self):
        if self.antiunitary:
            raise ValueError('antiunitary positive-time reuse is not supported')
        phase = _finite(self.phase, 'phase')
        _permutation(self.row_pullback, phase.size, 'PW row map')
        if (phase.ndim != 1 or not self.source_identity or not self.target_identity
                or np.max(np.abs(np.abs(phase)-1), initial=0.) > 2e-13):
            raise ValueError('declared unit-modulus phase and original identities required')

    def apply(self, vectors):
        vectors = np.asarray(vectors)
        if vectors.ndim not in (1, 2) or vectors.shape[0] != self.phase.size:
            raise ValueError('complete original PW rows required')
        return self.phase.reshape((-1,) + (1,)*(vectors.ndim-1))*vectors[self.row_pullback]

    def adjoint(self, vectors):
        vectors = np.asarray(vectors)
        if vectors.ndim not in (1, 2) or vectors.shape[0] != self.phase.size:
            raise ValueError('complete target PW rows required')
        result = np.empty_like(vectors, dtype=complex)
        result[self.row_pullback] = self.phase.conj().reshape((-1,) + (1,)*(vectors.ndim-1))*vectors
        return result


def fourier_witness(source_fft, target_fft, mapping, *, rotation, translation,
                    reciprocal_shift, source_k_fractional, fft_pullback):
    """Analytic columnwise ||Et U-W Es||F without constructing either embedding.

    W f(x)=exp(-2pi i[g.x+k.f]) f(S.T x-f). In the normalized grid
    inner product, different discrete Fourier frequencies are orthogonal.
    Nonmatching frequencies contribute 2 to the squared defect; matching ones
    contribute |actual_phase-expected_phase|^2. No huge-square subtraction.
    The *actual* supplied FFT permutation must equal the declared affine map.
    """
    shape = np.asarray(source_fft.physical, dtype=int)
    if (tuple(shape) != target_fft.physical or len(set(shape)) != 1
            or source_fft.miller.shape != target_fft.miller.shape
            or source_fft.miller.shape[1] != mapping.phase.size):
        raise ValueError('same cubic physical grid and full PW dimension required')
    s, f = np.asarray(rotation), _finite(translation, 'translation')
    g, k = np.asarray(reciprocal_shift), _finite(source_k_fractional, 'source k')
    if (s.shape != (3, 3) or s.dtype.kind not in 'iu'
            or not np.array_equal(s.T@s, np.eye(3, dtype=int))
            or f.shape != (3,) or g.shape != (3,) or g.dtype.kind not in 'iu' or k.shape != (3,)):
        raise ValueError('integer signed-permutation rotation/fold and finite translation/k required')
    shift = f*shape
    if not np.allclose(shift, np.rint(shift), atol=2e-10, rtol=0):
        raise ValueError('translation is not grid commensurate')
    pull = _permutation(fft_pullback, source_fft.nnr, 'FFT pullback')
    indices = np.indices(tuple(shape)).reshape(3, -1, order='F')
    expected_pull = np.ravel_multi_index(tuple((s.T@indices-np.rint(shift).astype(int)[:, None]) % shape[:, None]), tuple(shape), order='F')
    if not np.array_equal(pull, expected_pull):
        raise ValueError('FFT pullback differs from declared affine map')
    source_m = source_fft.miller[:, mapping.row_pullback]
    expected_m = s@source_m-g[:, None]
    same = np.all((target_fft.miller-expected_m) % shape[:, None] == 0, axis=0)
    # W uses the commensurate integer shift, not a repaired atomic geometry.
    grid_f = np.rint(shift)/shape
    expected_phase = np.exp(-2j*np.pi*(source_m.T@grid_f + k@f))
    squares = np.where(same, np.abs(mapping.phase-expected_phase)**2,
                       np.abs(mapping.phase)**2+1.)
    return dict(embedding_defect_frobenius=float(np.sqrt(np.sum(squares))),
                frequency_mismatch_columns=int(np.count_nonzero(~same)),
                phase_max_defect=float(np.max(np.abs(mapping.phase-expected_phase), initial=0.)),
                phase_modulus_roundoff=float(np.max(np.abs(np.abs(mapping.phase)-1), initial=0.)),
                grid_translation_roundoff=float(np.max(np.abs(grid_f-f))),
                fft_pullback_exact_integer_match=True)


def covariance_bound(source_factors, target_factors, mapping, source_potential,
                     target_potential, *, fft_pullback, fourier, atom_pullback,
                     angular_intertwiners=None):
    """Sum safe kinetic/local/per-atom low-rank bounds in Ry.

    Intertwiners map target projector columns to source columns: Bt=C Q+R,
    C=U Bs. They are declared, never least-squares fitted. Q need not be unitary:
    both its actual residual and Q Dt Q*-Ds are retained. Atom indices in
    PWFactors are one-based; atom_pullback is zero-based target -> source.
    """
    ts = _finite(source_factors.kinetic_ry, 'source kinetic')
    tt = _finite(target_factors.kinetic_ry, 'target kinetic')
    if ts.shape != mapping.phase.shape or tt.shape != ts.shape:
        raise ValueError('complete kinetic diagonals required')
    vs, vt = _finite(source_potential, 'source potential'), _finite(target_potential, 'target potential')
    pull = _permutation(fft_pullback, source_factors.fft.nnr, 'FFT pullback')
    if vs.shape != pull.shape or vt.shape != pull.shape or target_factors.fft.nnr != len(pull):
        raise ValueError('physical scalar potentials required')
    eta = float(fourier['embedding_defect_frobenius'])
    if not np.isfinite(eta) or eta < 0:
        raise ValueError('finite nonnegative Fourier defect required')
    kinetic = float(np.max(np.abs(tt-ts[mapping.row_pullback]), initial=0.))
    local_defect = float(np.max(np.abs(vt-vs[pull]), initial=0.))
    local_map = 2*eta*float(np.max(np.abs(vs), initial=0.))
    source_by_name = {x.name: x for x in source_factors.species}
    target_names = [x.name for x in target_factors.species]
    if len(source_by_name) != len(source_factors.species) or len(set(target_names)) != len(target_names) or set(target_names) != set(source_by_name):
        raise ValueError('exact original species required')
    atoms_s = [a for x in source_factors.species for a in x.atom_indices]
    atoms_t = [a for x in target_factors.species for a in x.atom_indices]
    nat = len(atoms_s)
    if sorted(atoms_s) != list(range(1, nat+1)) or sorted(atoms_t) != list(range(1, nat+1)):
        raise ValueError('every original atom exactly once required')
    atom_map = _permutation(atom_pullback, nat, 'atom pullback')
    intertwiners = {} if angular_intertwiners is None else angular_intertwiners
    if set(intertwiners)-set(target_names):
        raise ValueError('unknown intertwiner species')
    records = []
    for target in target_factors.species:
        source = source_by_name[target.name]
        ds, dt = _finite(source.active_d, 'source D'), _finite(target.active_d, 'target D')
        ns, nt = ds.shape[0], dt.shape[0]
        if ds.shape != (ns, ns) or dt.shape != (nt, nt):
            raise ValueError('square original active D blocks required')
        q = _finite(intertwiners.get(target.name, np.eye(ns, nt)), 'angular intertwiner')
        if q.shape != (ns, nt):
            raise ValueError('intertwiner dimensions differ')
        if (source.beta.shape != (ts.size, len(source.atom_indices)*ns)
                or target.beta.shape != (tt.size, len(target.atom_indices)*nt)):
            raise ValueError('full atom-major projector factors required')
        d_defect = float(np.linalg.norm(q@dt@q.conj().T-ds, 2))
        dnorm = float(np.linalg.norm(dt, 2))
        for j, atom in enumerate(target.atom_indices):
            source_atom = int(atom_map[atom-1])+1
            if source_atom not in source.atom_indices:
                raise ValueError('atom pullback changes species')
            i = source.atom_indices.index(source_atom)
            c = mapping.apply(source.beta[:, i*ns:(i+1)*ns])
            t = _finite(target.beta[:, j*nt:(j+1)*nt], 'target beta')
            a = c@q
            residual = float(np.linalg.norm(t-a))
            projector_term = (float(np.linalg.norm(t))+float(np.linalg.norm(a)))*residual*dnorm
            conjugacy_term = float(np.linalg.norm(c))**2*d_defect
            records.append(dict(species=target.name, target_atom=atom, source_atom=source_atom,
                projector_residual_frobenius=residual, d_conjugacy_defect_norm2=d_defect,
                projector_bound_ry=projector_term, d_conjugacy_bound_ry=conjugacy_term,
                bound_ry=projector_term+conjugacy_term,
                intertwiner='declared' if target.name in intertwiners else 'identity-not-fitted'))
    nonlocal_bound = sum(r['bound_ry'] for r in records)
    return dict(kinetic_bound_ry=kinetic, local_potential_defect_ry=local_defect,
                local_fourier_bound_ry=local_map, local_bound_ry=local_defect+local_map,
                nonlocal_bound_ry=nonlocal_bound, atoms=records, fourier=dict(fourier),
                total_bound_ry=kinetic+local_defect+local_map+nonlocal_bound,
                source_identity=mapping.source_identity, target_identity=mapping.target_identity,
                qualification='analytic upper inequality evaluated in floating point, not interval certified',
                exact_real_zero_implies_covariance=True, numerical_zero_proves_exact_symmetry=False)


def action_witness(source_action, target_action, mapping, vectors, *, block_size=8,
                   cached_source_action=None):
    """All supplied columns in full PW space; maxima are not norm certificates."""
    vectors = _finite(vectors, 'probe vectors')
    if vectors.ndim != 2 or vectors.shape[0] != mapping.phase.size or block_size < 1:
        raise ValueError('complete full-PW columns and positive block size required')
    if cached_source_action is not None and np.shape(cached_source_action) != vectors.shape:
        raise ValueError('cached source action shape differs')
    residuals, relative, norms = [], [], []
    for start in range(0, vectors.shape[1], block_size):
        block = vectors[:, start:start+block_size]
        source = source_action(block) if cached_source_action is None else cached_source_action[:, start:start+block_size]
        residual = target_action(mapping.apply(block))-mapping.apply(source)
        rn, vn = np.linalg.norm(residual, axis=0), np.linalg.norm(block, axis=0)
        if np.any(vn == 0) or not np.isfinite(rn).all():
            raise ValueError('finite residual and nonzero probes required')
        residuals.extend(rn.tolist()); relative.extend((rn/vn).tolist()); norms.extend(vn.tolist())
    return dict(columns=vectors.shape[1], column_norms=norms, residual_norms_ry=residuals,
                relative_residuals_ry=relative, residual_frobenius_ry=float(np.linalg.norm(residuals)),
                empirical_max_relative_ry=max(relative, default=0.),
                qualification='actual operator actions on specified probes, not an operator-norm upper bound')
