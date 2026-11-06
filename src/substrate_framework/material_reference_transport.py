"""Checked full-column source-orbit reuse, not a saved-background symmetry claim.

For the admitted cubic unitary convention u_t(x)=exp(-2pi i[g.x+k_s.f])
u_s(S.T x-f). Density therefore pulls back on the physical Fortran FFT mesh.
The positive ionic gradient at target atom a is S g_s[irt(a)]. The latter
is a covariance prediction at the transformed geometry, not an assertion
that rounded saved positions or the saved PBE background are invariant.
Antiunitary Bloch slots remain unsupported by the existing BlochBasis API.
"""
from dataclasses import dataclass

import numpy as np

from .material_reference_state import SourceState, source_population_indices
from .qe_bloch_frame import transport_bloch_frame


def raw_state_transport(state, *, antiunitary=False):
    """Raw-column algebra only: U K conjugates BOTH S and D, not populations.

    This helper does not admit antiunitary physical slots or spin reversal.
    Unitary transport shares immutable source arrays without a Gram rebuild.
    """
    if not antiunitary:
        return state
    gram, coefficient = state.gram.conj(), state.density_coefficients.conj()
    gram.flags.writeable = coefficient.flags.writeable = False
    return SourceState(gram, state.populations, coefficient, state.admissibility,
                       state.source_identity)


@dataclass(frozen=True)
class SlotTransport:
    slot: object
    fft_pullback: np.ndarray
    atom_pullback: np.ndarray
    cartesian_rotation: np.ndarray
    geometry_residual_bohr: np.ndarray

    def density(self, field):
        """Scalar or spin-channel field; no spin swap or quadrature factor."""
        field = np.asarray(field)
        if field.shape[0] != len(self.fft_pullback):
            raise ValueError('physical FFT field size differs')
        return field[self.fft_pullback]

    def gradient(self, positive_gradient):
        value = np.asarray(positive_gradient)
        if value.size != 3 * len(self.atom_pullback):
            raise ValueError('full original atom-major gradient required')
        return (value.reshape(-1, 3)[self.atom_pullback]
                @ self.cartesian_rotation.T).reshape(-1)

    def onsite(self, blocks):
        blocks = np.asarray(blocks)
        if blocks.shape != (len(self.atom_pullback), 3, 3):
            raise ValueError('all atom onsite Cartesian blocks required')
        rotation = self.cartesian_rotation
        return np.einsum('ij,ajk,lk->ail', rotation,
                         blocks[self.atom_pullback], rotation)

    def diagnostics(self, background=None):
        result = dict(geometry_max_bohr=float(np.max(np.abs(self.geometry_residual_bohr))),
                      h0_transport_supported=False,
                      gradient_scope='covariance at transformed geometry; discrepancies not repaired')
        if background is not None:
            result['background_scalar_max_residuals'] = {
                name: float(np.max(np.abs(self.density(getattr(background, name))
                                         - getattr(background, name))))
                for name in ('rho_updown', 'ionic_ry', 'hartree_ry', 'xc_ry', 'vrs')}
        return result


def slot_transport(grid, slot):
    """Build actual FFT/atom maps from an admitted slot, recording roundoff."""
    grid.require_slot(slot)
    ref, op = grid.reference, slot.operation
    if op.effective_time_reversal is not False:
        raise ValueError('physical antiunitary Bloch transport is not admitted')
    shape = np.asarray(ref.fft_smooth, dtype=np.int64)
    s, f = op.rotation, op.fractional_translation
    if (not np.array_equal(s.T @ s, np.eye(3, dtype=int))
            or len(set(shape)) != 1):
        raise ValueError('checked signed-permutation equal FFT grid required')
    shift = f * shape
    if not np.allclose(shift, np.rint(shift), atol=2e-10, rtol=0):
        raise ValueError('fractional translation is not FFT commensurate')
    indices = np.indices(tuple(shape)).reshape(3, -1, order='F')
    pulled = (s.T @ indices - np.rint(shift).astype(np.int64)[:, None]) % shape[:, None]
    fft_map = np.ravel_multi_index(tuple(pulled), tuple(shape), order='F')
    atoms = np.asarray(op.equivalent_atoms, dtype=np.int64) - 1
    if not np.array_equal(np.sort(atoms), np.arange(ref.nat)):
        raise ValueError('saved equivalent atoms are not a permutation')
    if any(ref.atom_names[a] != ref.atom_names[b] for a, b in enumerate(atoms)):
        raise ValueError('saved atom map changes species')
    cell = ref.cell_bohr
    rotation = cell @ s @ np.linalg.inv(cell)
    x = np.linalg.solve(cell, ref.positions_bohr)
    residual = cell @ (s.T @ x - f[:, None] - x[:, atoms] - slot.atom_image_offsets)
    return SlotTransport(slot, fft_map, atoms, rotation, residual)


def transport_source(grid, slot, source, state, *, source_admission=None):
    """Return checked full raw frame and SAME unitary Gram/D source state."""
    grid.require_slot(slot)
    if (state.source_identity != source.source['sha256']
            or state.gram.shape != (source.nbnd, source.nbnd)):
        raise ValueError('source state and all raw WFC columns differ')
    if slot.operation.effective_time_reversal is not False:
        raise ValueError('antiunitary slots require a separate admitted Bloch basis')
    frame = transport_bloch_frame(grid.reference, slot, source,
                                  source_admission=source_admission)
    return frame, raw_state_transport(state)


@dataclass(frozen=True)
class SourceObservables:
    """Unweighted moments at the original identity representative, all columns."""
    representative: object
    state: SourceState
    density_bohr3: np.ndarray
    positive_gradient_ry_bohr: np.ndarray
    conditional_nonlocal_ry_bohr2: np.ndarray | None = None
    original_mv_nonlocal_ry_bohr2: np.ndarray | None = None
    mv_density_bohr3: np.ndarray | None = None


def aggregate_complete_slots(grid, ledger, sources, *,
                             gradient_policy='direct-unless-exact', direct_observables=None):
    """Consume each source once; retain explicit every-slot weights and spin.

    sources is an iterable of SourceObservables (it may stream checkpoints).
    Default uses covariance only for the original identity map; otherwise
    direct_observables(slot, state) must return gradient and, when requested,
    conditional nonlocal blocks. Explicit 'declared-geometric' permits the measured
    geometry/operator covariance approximation, retaining its residual. Density
    and MV density are exact discrete frame pullbacks in both modes. No H0,
    PBE kernel, held density, ionic force or electron-ion total is inferred.
    MV nonlocal blocks are subtracted ONCE per original source using its XML
    IBZ weight (raw occupations=wg/wk), never transported without proof.
    """
    if gradient_policy not in ('direct-unless-exact', 'declared-geometric'):
        raise ValueError('explicit transport policy required')
    nr, nat = int(np.prod(grid.reference.fft_smooth)), grid.reference.nat
    nspin = len(grid.reference.nbnd_by_spin)
    density = np.zeros((nr, nspin), complex)
    mv_density = np.zeros_like(density)
    gradient = np.zeros(3 * nat, complex)
    onsite = np.zeros((nat, 3, 3), complex)
    original_mv = np.zeros_like(onsite)
    visited, reports, optional_scope = set(), [], None
    for item in sources:
        rep, state = item.representative, item.state
        grid.require_slot(rep)
        if (not np.array_equal(rep.operation.rotation, np.eye(3, dtype=int))
                or np.any(rep.operation.fractional_translation) or np.any(rep.reciprocal_shift)
                or state.source_identity != rep.source_binding['sha256']):
            raise ValueError('moments must belong to the original identity source frame')
        population, pairs = source_population_indices(grid, ledger, rep.ispin, rep.source_ik)
        if not np.array_equal(population, state.populations):
            raise ValueError('source state populations differ from complete slot ledger')
        scope = (item.conditional_nonlocal_ry_bohr2 is not None, item.mv_density_bohr3 is not None)
        if (item.original_mv_nonlocal_ry_bohr2 is not None) != scope[0]:
            raise ValueError('conditional and original MV nonlocal blocks must both be supplied')
        if optional_scope is not None and optional_scope != scope:
            raise ValueError('partial-source optional curvature/MV coverage is forbidden')
        optional_scope = scope
        if scope[0]:
            if (np.shape(item.original_mv_nonlocal_ry_bohr2) != onsite.shape
                    or np.shape(item.conditional_nonlocal_ry_bohr2) != onsite.shape):
                raise ValueError('complete onsite source blocks required')
            point, _, _ = grid.reference.select(rep.source_ik, rep.ispin)
            original_mv += point.weight * item.original_mv_nonlocal_ry_bohr2
        for _, slot in pairs:
            key = (slot.ispin, slot.mesh_index)
            if key in visited:
                raise ValueError('duplicate spin/mesh slot contribution')
            visited.add(key)
            transform = slot_transport(grid, slot)
            report = transform.diagnostics()
            exact_geometry = (not np.any(transform.geometry_residual_bohr)
                              and np.array_equal(slot.operation.rotation, np.eye(3, dtype=int))
                              and not np.any(slot.operation.fractional_translation))
            if exact_geometry or gradient_policy == 'declared-geometric':
                g = transform.gradient(item.positive_gradient_ry_bohr)
                c = transform.onsite(item.conditional_nonlocal_ry_bohr2) if scope[0] else None
                report['evaluation'] = 'geometric-covariance'
            else:
                if direct_observables is None:
                    raise ValueError('rounded geometry requires direct slot evaluation or explicit approximation')
                actual = direct_observables(slot, state)
                g = np.asarray(actual['positive_gradient_ry_bohr'])
                c = np.asarray(actual['conditional_nonlocal_ry_bohr2']) if scope[0] else None
                if g.shape != gradient.shape or (scope[0] and c.shape != onsite.shape):
                    raise ValueError('direct complete-coordinate moments have wrong shape')
                report['evaluation'] = 'direct-at-saved-geometry'
            density[:, slot.ispin-1] += slot.weight * transform.density(item.density_bohr3)
            gradient += slot.weight * g
            if scope[0]:
                onsite += slot.weight * c
            if scope[1]:
                mv_density[:, slot.ispin-1] += slot.weight * transform.density(item.mv_density_bohr3)
            report.update(spin=slot.ispin, mesh_index=list(slot.mesh_index),
                          source_ik=slot.source_ik, source_sha256=state.source_identity,
                          measure=slot.weight)
            reports.append(report)
    expected = {(slot.ispin, slot.mesh_index) for slot in grid.slots}
    if visited != expected:
        raise ValueError('full spin/mesh source coverage is incomplete')
    return dict(density_bohr3=density, positive_gradient_ry_bohr=gradient,
                nonlocal_correction_ry_bohr2=onsite-original_mv if optional_scope[0] else None,
                conditional_nonlocal_ry_bohr2=onsite if optional_scope[0] else None,
                original_mv_nonlocal_ry_bohr2=original_mv if optional_scope[0] else None,
                mv_density_bohr3=mv_density if optional_scope[1] else None,
                slots=reports, gradient_policy=gradient_policy,
                h0_transport_supported=False)
