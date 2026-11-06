"""Full raw-column physical UPF field moments without PW-square operators.

A28 full-field-moments-registration specifies the scalar and covariance.
Only fixed-state zero-field coefficients are returned. No H0 transport,
anchored undressed offset, equilibrium, field state or SCF response is inferred.
"""
from types import SimpleNamespace

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_material_field_derivatives import UPFRadial, centered_derivatives
from .material_reference_transport import slot_transport


AXES = np.eye(3)
PAIRS = tuple((i, j) for i in range(3) for j in range(i, 3))


def cartesian_jets(radial, momenta):
    """Ten centered value/first/symmetric-second jets, ONE radial evaluation."""
    values = radial.evaluate(np.einsum('pi,pi->p', momenta, momenta))
    cached = SimpleNamespace(tables=radial.tables, evaluate=lambda s: values)
    jets = [None] * 10
    for index, (i, j) in enumerate(PAIRS):
        b, bi, bj, bij = centered_derivatives(cached, momenta, AXES[i], AXES[j])
        if index == 0:
            jets[0] = b
        if i == j:
            jets[1+i] = bi
        jets[4+index] = bij
    return jets


def ionic_jets(beta, momenta):
    return [beta] + [-1j*momenta[:, i, None]*beta for i in range(3)] + [
        -momenta[:, i, None]*momenta[:, j, None]*beta for i, j in PAIRS]


def overlap_moments(overlaps, active_d, density_coefficients):
    """Ordered complex trace, independent left/right legs; D need not be real."""
    x = np.asarray(overlaps)
    if x.ndim != 3 or len(x) != 10:
        raise ValueError('ten projector-overlap jets required')
    weighted = np.stack([active_d @ a @ density_coefficients for a in x])
    def pair(i, j):
        return np.vdot(x[i], weighted[j])
    gradient = np.array([pair(1+i, 0)+pair(0, 1+i) for i in range(3)])
    curvature = np.empty((3, 3), complex)
    for index, (i, j) in enumerate(PAIRS):
        value = (pair(4+index, 0)+pair(1+i, 1+j)
                 +pair(1+j, 1+i)+pair(0, 4+index))
        curvature[i, j] = curvature[j, i] = value
    return pair(0, 0), gradient, curvature


def require_scalar_channels(upf, active_d, tolerance=2e-12):
    """Ordinary scalar NC d: radial mixing identical for every m within l.

    This explicit block test implies commutation with all real-harmonic O(3)
    representations, including improper signed permutations; no spin action.
    """
    d = np.asarray(active_d)
    expected = np.zeros_like(d)
    for a in range(upf.nbeta):
        rows = np.flatnonzero(upf.channel_radial == a)
        for b in range(upf.nbeta):
            cols = np.flatnonzero(upf.channel_radial == b)
            if upf.angular_l[a] == upf.angular_l[b]:
                if not np.array_equal(upf.channel_lm[rows], upf.channel_lm[cols]):
                    raise ValueError('ordinary QE channel ordering required')
                expected[rows, cols] = d[rows[0], cols[0]]
    residual = float(np.max(np.abs(d-expected)))
    if residual > tolerance*max(1., float(np.max(np.abs(d)))):
        raise ValueError('nonlocal channel matrix is not rotationally scalar')
    return residual


def full_column_moments(hamiltonian, frame, density_coefficients, *, row_block=256):
    """Exact retained radial quadrature, all PW rows and independent raw columns.

    Returns unweighted one-spin/original-source means. Work is O(P*N^2)
    kinetic plus O(P*N*29*sum_atoms(nh)); no local-H FFT or P-square object.
    Delta nonlocal curvature is atom onsite: all off-atom entries are exactly
    zero for the fixed-state one-center nonlocal scalar, not for SCF response.
    """
    from .qe_material_current import MaterialOpticalCurrent
    psi = np.asarray(frame)
    d = np.asarray(density_coefficients)
    p = MaterialOpticalCurrent(hamiltonian).momenta
    if (psi.ndim != 2 or psi.shape[0] != len(p) or d.shape != (psi.shape[1],)*2
            or not np.isfinite(d).all() or type(row_block) is not int or row_block < 1):
        raise ValueError('full frame, finite independent D and positive row block required')
    inputs = hamiltonian.background.inputs
    nat = inputs.reference.positions_bohr.shape[1]
    kinetic_gradient = np.zeros(3, complex)
    number = 0j
    for start in range(0, len(p), row_block):
        stop = min(start+row_block, len(p))
        v = psi[start:stop]
        if not np.isfinite(v).all():
            raise ValueError('nonfinite raw frame')
        diagonal = np.einsum('pi,pi->p', v.conj(), v @ d)
        number += diagonal.sum()
        kinetic_gradient += 2*p[start:stop].T @ diagonal
    gradient = kinetic_gradient.copy()
    contact = 2*number*np.eye(3, dtype=complex)
    delta_energy = 0j
    delta_gradient = np.zeros((nat, 3), complex)
    delta_onsite = np.zeros((nat, 3, 3), complex)
    checks, overlap_columns = {}, 0
    for species in hamiltonian.factors.species:
        tables = inputs.tables[species.name]
        checks[species.name] = require_scalar_channels(tables.upf, species.active_d)
        nh = species.active_d.shape[0]
        atoms = tuple(species.atom_indices)
        radial = UPFRadial(tables, row_block=row_block)
        # The zero jet is shared by physical field and ionic families.
        overlaps = np.zeros((len(atoms), 29, nh, psi.shape[1]), complex)
        overlap_columns += len(atoms)*29*nh
        for start in range(0, len(p), row_block):
            stop = min(start+row_block, len(p))
            block_p = p[start:stop]
            vectors = np.asfortranarray(psi[start:stop])
            centered = cartesian_jets(radial, block_p)
            for slot, atom in enumerate(atoms):
                tau = inputs.reference.positions_bohr[:, atom-1]
                phase = np.exp(-1j*(block_p @ tau))[:, None]
                field = [b*phase for b in centered]
                physical_ionic = ionic_jets(field[0], block_p)
                legacy = species.beta[start:stop, slot*nh:(slot+1)*nh]
                jets = field + physical_ionic[1:] + ionic_jets(legacy, block_p)
                wide = np.asfortranarray(np.concatenate(jets, axis=1))
                overlaps[slot] += zgemm(1., wide, vectors, trans_a=2).reshape(
                    29, nh, psi.shape[1])
        for slot, atom in enumerate(atoms):
            ov = overlaps[slot]
            _, g, c = overlap_moments(ov[:10], species.active_d, d)
            gradient += g
            contact += c
            physical = overlap_moments(np.concatenate((ov[:1], ov[10:19])), species.active_d, d)
            legacy = overlap_moments(ov[19:29], species.active_d, d)
            delta_energy += physical[0]-legacy[0]
            delta_gradient[atom-1] += physical[1]-legacy[1]
            delta_onsite[atom-1] += physical[2]-legacy[2]
    return dict(number=number, gradient=gradient, current_hartree=-gradient/2,
                contact=contact, kinetic_gradient=kinetic_gradient,
                kinetic_contact=2*number*np.eye(3),
                delta_nonlocal_energy_ry=delta_energy,
                delta_nonlocal_positive_gradient_ry_bohr=delta_gradient,
                delta_nonlocal_force_ry_bohr=-delta_gradient,
                delta_nonlocal_onsite_ry_bohr2=delta_onsite,
                rotational_scalar_residuals=checks,
                work=dict(pw_rows=len(p), raw_columns=psi.shape[1],
                          wide_projector_columns=overlap_columns,
                          complex_overlap_multiply_adds=len(p)*psi.shape[1]*overlap_columns,
                          kinetic_complex_multiply_adds=len(p)*psi.shape[1]**2,
                          full_frame_bytes=int(psi.nbytes), local_fft_count=0,
                          pw_square_matrices=0),
                normalization='unweighted original spin/source; no spin degeneracy',
                family='unanchored physical UPF radial quadrature at zero field')


MOMENT_KEYS = ('number', 'gradient', 'current_hartree', 'contact',
               'kinetic_gradient', 'kinetic_contact', 'delta_nonlocal_energy_ry',
               'delta_nonlocal_positive_gradient_ry_bohr', 'delta_nonlocal_force_ry_bohr',
               'delta_nonlocal_onsite_ry_bohr2')


def transport_moments(moments, rotation, atom_pullback):
    """Tensor covariance at the transformed geometry; NOT H0 covariance."""
    r = np.asarray(rotation)
    if r.shape != (3, 3) or not np.allclose(r.T @ r, np.eye(3), atol=2e-12, rtol=0):
        raise ValueError('orthogonal Cartesian rotation required')
    out = {key: np.asarray(moments[key]).copy() for key in MOMENT_KEYS}
    for key in ('gradient', 'current_hartree', 'kinetic_gradient'):
        out[key] = r @ out[key]
    for key in ('contact', 'kinetic_contact'):
        out[key] = r @ out[key] @ r.T
    for key in ('delta_nonlocal_positive_gradient_ry_bohr', 'delta_nonlocal_force_ry_bohr'):
        out[key] = out[key][atom_pullback] @ r.T
    out['delta_nonlocal_onsite_ry_bohr2'] = np.einsum(
        'ij,ajk,lk->ail', r, out['delta_nonlocal_onsite_ry_bohr2'][atom_pullback], r)
    return out


def aggregate_slots(grid, source_moments, *, geometry_policy):
    """Every slot exactly once; source inputs MUST remain unweighted.

    source_moments maps (spin, original_ik) to full_column_moments results.
    Original XML weighted diagnostics are returned separately, never combined.
    """
    if geometry_policy != 'declared-geometric':
        raise ValueError('explicit declared-geometric optical approximation required')
    expected = {(s.ispin, s.source_ik) for s in grid.slots}
    if set(source_moments) != expected:
        raise ValueError('all and only original sources required')
    template = next(iter(source_moments.values()))
    total = {k: np.zeros_like(template[k], dtype=complex) for k in MOMENT_KEYS}
    original = {k: np.zeros_like(template[k], dtype=complex) for k in MOMENT_KEYS}
    for key, item in source_moments.items():
        if item.get('normalization') != 'unweighted original spin/source; no spin degeneracy':
            raise ValueError('weighted or unspecified source normalization forbidden')
        point, _, _ = grid.reference.select(key[1], key[0])
        for name in MOMENT_KEYS:
            original[name] += point.weight*np.asarray(item[name])
    reports, visited = [], set()
    for slot in grid.slots:
        key = (slot.ispin, slot.mesh_index)
        if key in visited:
            raise ValueError('duplicate full-grid slot')
        visited.add(key)
        transform = slot_transport(grid, slot)
        values = transport_moments(source_moments[(slot.ispin, slot.source_ik)],
                                   transform.cartesian_rotation, transform.atom_pullback)
        for name in MOMENT_KEYS:
            total[name] += slot.weight*values[name]
        reports.append(dict(spin=slot.ispin, source_ik=slot.source_ik,
                            mesh_index=list(slot.mesh_index), slot_weight=slot.weight,
                            **transform.diagnostics()))
    return dict(moments=total, original_source_weighted_diagnostic=original,
                slots=reports, geometry_policy=geometry_policy,
                normalization='full spin/mesh mean; slot weight applied exactly once',
                h0_transport_supported=False, field_state=None,
                correction_scope='fixed-state zero-field nonlocal scalar only; no SCF or finite-field ionic join')
