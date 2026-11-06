"""Fixed-state ordinary-NC onsite electronic curvature, in Ry/bohr**2.

QE7.5 PHonon/PH/dynmat_us.f90:92-123,144-165,209-219 and
PW/src/compute_deff.f90:35-42. The admitted BlochIonicAction supplies the SAME
placed beta, active D, k+G, density G, local radial tables and physical FFT as
the first derivative. No screened, Pulay, NLCC, US/PAW, ionic Ewald, ASR or
mass terms enter. There are no inter-atom electronic explicit terms here.

The four nonlocal terms are B_ij D B† + B D B_ij† + B_i D B_j† + B_j D B_i†,
where B_i=-i p_i B and B_ij=-p_i p_j B. Both second legs and both cross terms
are retained independently, even for non-Hermitian numerical dyads. No real
part or Hermitian projection is used. Six Cartesian pairs suffice because
these explicit partial derivatives commute, NOT because of symmetrization.

The complete material conversion is local[n_cond-full - n_saved-QE] plus
sum_full_slots weight*nonlocal[D_cond] - sum_original_sources XML_weight*
nonlocal[diag(raw-MV)]. Modules/qexsd_init.f90:1247-1249 writes wg/wk.
Saved MV coefficients are numerical reference coefficients, never probabilities.
The orbital-MV local subtraction is a separate diagnostic, not a replacement
for saved QE rho. Ion-ion curvature is unchanged by electronic preparation.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg.blas import zgemm


@dataclass(frozen=True)
class OnsiteCurvature:
    local: np.ndarray
    nonlocal_term: np.ndarray

    @property
    def total(self):
        return self.local + self.nonlocal_term

    def matrix(self):
        """Small atom-major coordinate matrix; offsite explicit terms are zero."""
        blocks = self.total
        result = np.zeros((3*len(blocks), 3*len(blocks)), complex)
        for atom, block in enumerate(blocks):
            result[3*atom:3*atom+3, 3*atom:3*atom+3] = block
        return result


def local_density_curvature(density_bohr3, coupling):
    """C_local[n] for total physical density on the ORIGINAL signed-G support.

    Input can be complex and need not be a preparation. It is one scalar total
    density, so callers sum spin channels themselves. Omega sum_G n(-G)
    [-G_i G_j V_species(G) exp(-i G.R_a)] is the bilinear trace, using Ry,
    physical bohr^-1, radial tables' original 1/Omega and Fortran FFT order.
    No q, alpha or Ewald cutoff occurs in this ELECTRONIC onsite functional.
    """
    ionic, fft = coupling.ionic, coupling.factors.fft
    density = np.asarray(density_bohr3)
    if density.shape != (fft.nnr,) or not np.isfinite(density).all():
        raise ValueError('one full physical total-density field required')
    transposed = np.fft.ifftn(density.reshape(fft.physical, order='F'))
    miller = ionic.inputs.density.miller
    sampled = transposed[tuple(miller[i] % fft.physical[i] for i in range(3))]
    ref, g = ionic.inputs.reference, ionic.g
    out = np.empty((ref.nat, 3, 3), complex)
    moduli = np.linalg.norm(g, axis=1)
    for atom, name in enumerate(ref.atom_names):
        if name not in ionic.vloc:
            ionic.vloc[name] = ionic.inputs.tables[name].vloc(moduli)
        weighted = sampled * ionic.vloc[name] * np.exp(-1j*(g @ ref.positions_bohr[:, atom]))
        out[atom] = -fft.omega * (g.T @ (weighted[:, None] * g))
    return out


def _weighted_frame(psi, coefficient, coupling):
    psi, coefficient = np.asarray(psi), np.asarray(coefficient)
    if (psi.ndim != 2 or psi.shape[0] != coupling.npw
            or coefficient.shape != (psi.shape[1], psi.shape[1])
            or not np.isfinite(psi).all() or not np.isfinite(coefficient).all()):
        raise ValueError('finite full PW raw frame and square raw-column coefficients required')
    return psi, zgemm(1., psi, coefficient)


def dyad_density(left, right, coupling):
    """Physical density of L R†; only FFT columns blocked, no band omission."""
    if left.shape != right.shape or left.shape[0] != coupling.npw:
        raise ValueError('same complete PW dyad factors required')
    fft = coupling.factors.fft
    density = np.zeros(fft.nnr, complex)
    for start in range(0, left.shape[1], coupling.block_size):
        stop = start + coupling.block_size
        l = fft.from_g(left[:, start:stop])
        r = fft.from_g(right[:, start:stop])
        density += np.sum(l*r.conj(), axis=1)/fft.omega
    return density


def nonlocal_dyad_curvature(left, right, coupling, *, diagonal_reference=None):
    """All atom blocks from wide species overlaps, never PW-square or 3N actions.

    Only base/three first moments and one second moment per side/species live
    simultaneously. The raw-frame orthogonal complement has zero rho but is
    NOT projected out of the operator; no inserted saved-band identity appears.
    Optional diagonal_reference returns a second contraction R diag(w) R†,
    reusing the right overlaps. It does not add another PW matrix product.
    """
    if left.shape != right.shape or left.shape[0] != coupling.npw:
        raise ValueError('same complete PW dyad factors required')
    out = np.zeros((coupling.ionic.inputs.reference.nat, 3, 3), complex)
    reference = None
    if diagonal_reference is not None:
        weights = np.asarray(diagonal_reference)
        if weights.shape != (right.shape[1],) or not np.isfinite(weights).all():
            raise ValueError('finite complete raw-column numerical weights required')
        reference = np.zeros_like(out)
    p = coupling.ionic.p
    for species in coupling.factors.species:
        atoms = np.asarray(species.atom_indices, dtype=int) - 1
        shape = (len(atoms), species.active_d.shape[0], left.shape[1])
        beta, d = species.beta, species.active_d

        def overlap(frame, momentum=None):
            value = frame if momentum is None else momentum[:, None]*frame
            return zgemm(1., beta, value, trans_a=2).reshape(shape)

        def apply_d(value):
            return np.einsum('ij,ajb->aib', d, value, optimize=True)

        def pair(bra, ket):
            return np.einsum('aib,aib->a', bra.conj(), ket, optimize=True)

        bl, br = overlap(left), overlap(right)
        dbl = apply_d(bl)
        first_l = [apply_d(overlap(left, p[:, axis])) for axis in range(3)]
        first_r = [overlap(right, p[:, axis]) for axis in range(3)]
        if reference is not None:
            dbr = apply_d(br) * weights
            first_reference = [apply_d(value) * weights for value in first_r]
        for i in range(3):
            for j in range(i, 3):
                second_l = apply_d(overlap(left, p[:, i]*p[:, j]))
                second_r = overlap(right, p[:, i]*p[:, j])
                value = (-pair(br, second_l) - pair(second_r, dbl)
                         + pair(first_r[i], first_l[j]) + pair(first_r[j], first_l[i]))
                out[atoms, i, j] = value
                if i != j:
                    out[atoms, j, i] = value
                if reference is not None:
                    second_reference = apply_d(second_r) * weights
                    reference_value = (-pair(br, second_reference) - pair(second_r, dbr)
                        + pair(first_r[i], first_reference[j]) + pair(first_r[j], first_reference[i]))
                    reference[atoms, i, j] = reference_value
                    if i != j:
                        reference[atoms, j, i] = reference_value
    return out if reference is None else (out, reference)


def frozen_onsite_curvature(psi, coefficient, coupling, *, density_bohr3=None):
    """Unweighted full-column fixed-state expectation; no MV interpretation."""
    psi, left = _weighted_frame(psi, coefficient, coupling)
    density = dyad_density(left, psi, coupling) if density_bohr3 is None else density_bohr3
    return OnsiteCurvature(local_density_curvature(density, coupling),
                           nonlocal_dyad_curvature(left, psi, coupling))


def conditional_minus_mv(psi, state, coupling, raw_mv_coefficients, *, conditional_density=None):
    """One SOURCE's unweighted numerical orbital-MV difference and diagnostics.

    Nonlocal conditional and MV blocks are exported separately: conditional
    blocks use every full-grid slot; original MV blocks use XML IBZ weights
    at original representatives, exactly as dynmat_us uses wg. The per-source
    difference and orbital local difference are diagnostics, not full totals.
    Use common_reference_correction for the actual saved density subtraction.
    Zero/signed raw MV entries survive; no occupation clipping.
    """
    psi = np.asarray(psi)
    mv = np.asarray(raw_mv_coefficients)
    if (state.source_identity != coupling.ionic.basis.slot.source_binding['sha256']
            or mv.shape != (psi.shape[1],) or not np.isfinite(mv).all()
            or not np.array_equal(mv, coupling.ionic.basis.slot.raw_occupations)):
        raise ValueError('exact same-source saved XML numerical MV coefficients required')
    delta = np.array(state.density_coefficients, copy=True)
    delta[np.diag_indices_from(delta)] -= mv
    psi, left = _weighted_frame(psi, delta, coupling)
    nonlocal_correction, mv_nonlocal = nonlocal_dyad_curvature(
        left, psi, coupling, diagonal_reference=mv)
    # A saved conditional density avoids repeating its FFTs. MV is diagonal,
    # so its frame requires no wide raw-column coefficient multiplication.
    mv_density = dyad_density(psi * mv[None, :], psi, coupling)
    if conditional_density is None:
        conditional_density = dyad_density(left, psi, coupling) + mv_density
    orbital_local = local_density_curvature(conditional_density - mv_density, coupling)
    return dict(nonlocal_correction_ry_bohr2=nonlocal_correction,
                conditional_nonlocal_ry_bohr2=nonlocal_correction + mv_nonlocal,
                original_mv_nonlocal_ry_bohr2=mv_nonlocal,
                orbital_local_correction_ry_bohr2=orbital_local,
                orbital_total_correction_ry_bohr2=orbital_local + nonlocal_correction,
                mv_density_bohr3=mv_density,
                conditional_density_bohr3=np.asarray(conditional_density),
                mv_number=np.asarray(np.dot(mv, np.diag(state.gram))),
                conditional_number=np.asarray(np.einsum('ij,ji->', state.density_coefficients, state.gram)))


def common_reference_correction(conditional_density_updown, saved_density_updown,
                                full_nonlocal_correction, coupling):
    """Actual preparation conversion: local saved-rho difference plus NL delta.

    Both densities must already carry the complete spin/full-k measures.
    full_nonlocal_correction is the complete conditional-slot sum minus the
    original XML-weighted source MV sum. This function cannot authenticate
    external coverage; use aggregate_complete_slots
    and the original accepted saved background, not a selected source density.
    Add its matrix to the original MV frozen curvature without any ASR repair.
    The original ionic-q alpha/cutoffs/phases are untouched by this correction.
    """
    conditional, saved = np.asarray(conditional_density_updown), np.asarray(saved_density_updown)
    if (conditional.shape != (coupling.nr, 2) or saved.shape != conditional.shape
            or np.shape(full_nonlocal_correction) != (coupling.ionic.inputs.reference.nat, 3, 3)):
        raise ValueError('complete up/down physical fields and onsite blocks required')
    local = local_density_curvature(np.sum(conditional - saved, axis=1), coupling)
    return OnsiteCurvature(local, np.asarray(full_nonlocal_correction))
