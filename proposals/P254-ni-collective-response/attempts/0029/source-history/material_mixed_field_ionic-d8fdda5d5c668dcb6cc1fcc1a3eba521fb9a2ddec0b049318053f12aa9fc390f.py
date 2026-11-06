"""Direct fixed-state mixed UPF ionic/field coefficient (A28 registration).

b(p,a,u)=exp(-ip.(tau+u))*b_center(p+a); atom common Peierls phases
cancel in b d b^dagger. No density derivative, SCF response or legacy
undressed offset is hidden in the mixed coefficient. All original raw
columns and the independently supplied complex D remain in the trace.
"""
from types import SimpleNamespace

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_material_field_derivatives import UPFRadial, centered_derivatives
from .qe_material_current import MaterialOpticalCurrent
from .material_physical_field_moments import require_scalar_channels
from .material_reference_transport import slot_transport


NORMALIZATION = 'unweighted original spin/source; no spin degeneracy'


def first_cartesian_jets(radial, momenta):
    """Four centered jets, one retained radial quadrature evaluation."""
    values = radial.evaluate(np.einsum('pi,pi->p', momenta, momenta))
    cached = SimpleNamespace(tables=radial.tables, evaluate=lambda s: values)
    jets = []
    for axis in np.eye(3):
        b, ba, _, _ = centered_derivatives(cached, momenta, axis)
        if not jets:
            jets.append(b)
        jets.append(ba)
    return jets


def mixed_projector_jets(beta, field_first, momenta):
    """Order: b, b_a[3], b_u[3], b_ua[3,3] (ionic-major)."""
    p = np.asarray(momenta)
    b = np.asarray(beta)
    a = np.asarray(field_first)
    if b.ndim != 2 or p.shape != (len(b), 3) or a.shape != (3, *b.shape):
        raise ValueError('three field jets and matching PW momenta required')
    return [b, *a, *[-1j*p[:, i, None]*b for i in range(3)],
            *[-1j*p[:, i, None]*a[j] for i in range(3) for j in range(3)]]


def mixed_overlap_coefficient(overlaps, active_d, density_coefficients):
    """Four independent ordered complex terms; no real projection/repair."""
    x = np.asarray(overlaps)
    d = np.asarray(active_d)
    rho = np.asarray(density_coefficients)
    if (x.ndim != 3 or len(x) != 16 or d.shape != (x.shape[1],)*2
            or rho.shape != (x.shape[2],)*2):
        raise ValueError('16 complete overlap jets and matching d,D required')
    weighted = np.stack([d @ item @ rho for item in x])
    pair = lambda i, j: np.vdot(x[i], weighted[j])
    return np.asarray([[pair(7+3*i+j, 0)+pair(4+i, 1+j)
                        +pair(1+j, 4+i)+pair(0, 7+3*i+j)
                        for j in range(3)] for i in range(3)])


def full_column_mixed(hamiltonian, frame, density_coefficients, *, row_block=256):
    """NEW 16-jet overlaps; existing scalar moment receipts cannot supply them."""
    psi = np.asarray(frame)
    d = np.asarray(density_coefficients)
    p = MaterialOpticalCurrent(hamiltonian).momenta
    if (psi.ndim != 2 or psi.shape[0] != len(p) or d.shape != (psi.shape[1],)*2
            or not np.isfinite(d).all() or type(row_block) is not int or row_block < 1):
        raise ValueError('full frame, independent finite D and positive row block required')
    inputs = hamiltonian.background.inputs
    mixed = np.zeros((inputs.reference.positions_bohr.shape[1], 3, 3), complex)
    checks, columns = {}, 0
    for species in hamiltonian.factors.species:
        tables = inputs.tables[species.name]
        checks[species.name] = require_scalar_channels(tables.upf, species.active_d)
        nh, atoms = species.active_d.shape[0], tuple(species.atom_indices)
        radial = UPFRadial(tables, row_block=row_block)
        overlaps = np.zeros((len(atoms), 16, nh, psi.shape[1]), complex)
        columns += len(atoms)*16*nh
        for start in range(0, len(p), row_block):
            stop = min(start+row_block, len(p))
            block_p = p[start:stop]
            vectors = np.asfortranarray(psi[start:stop])
            if not np.isfinite(vectors).all():
                raise ValueError('nonfinite raw frame')
            centered = first_cartesian_jets(radial, block_p)
            for slot, atom in enumerate(atoms):
                tau = inputs.reference.positions_bohr[:, atom-1]
                phase = np.exp(-1j*(block_p @ tau))[:, None]
                placed = [b*phase for b in centered]
                jets = mixed_projector_jets(placed[0], placed[1:], block_p)
                wide = np.asfortranarray(np.concatenate(jets, axis=1))
                overlaps[slot] += zgemm(1., wide, vectors, trans_a=2).reshape(
                    16, nh, psi.shape[1])
        for slot, atom in enumerate(atoms):
            mixed[atom-1] = mixed_overlap_coefficient(overlaps[slot], species.active_d, d)
    return dict(mixed_ionic_field=mixed, normalization=NORMALIZATION,
                rotational_scalar_residuals=checks,
                work=dict(pw_rows=len(p), raw_columns=psi.shape[1],
                          wide_projector_columns=columns,
                          complex_overlap_multiply_adds=len(p)*psi.shape[1]*columns,
                          overlap_workspace_upper_bytes=16*columns*psi.shape[1],
                          full_frame_bytes=int(psi.nbytes), local_fft_count=0,
                          pw_square_matrices=0, saved_overlap_reuse=False))


def aggregate_mixed(grid, sources, *, geometry_policy):
    """Original sources to all independent slots; not a claim of H0 covariance."""
    if geometry_policy != 'declared-geometric':
        raise ValueError('explicit declared-geometric approximation required')
    expected = {(s.ispin, s.source_ik) for s in grid.slots}
    if set(sources) != expected:
        raise ValueError('all and only original mixed sources required')
    nat = grid.reference.positions_bohr.shape[1]
    total = np.zeros((nat, 3, 3), complex)
    visited, diagnostics = set(), []
    for slot in grid.slots:
        key = (slot.ispin, slot.mesh_index)
        if key in visited:
            raise ValueError('duplicate full-grid slot')
        visited.add(key)
        source = sources[(slot.ispin, slot.source_ik)]
        if source.get('normalization') != NORMALIZATION:
            raise ValueError('source weights/spin degeneracy must not be pre-applied')
        value = np.asarray(source['mixed_ionic_field'])
        if value.shape != (nat, 3, 3) or not np.isfinite(value).all():
            raise ValueError('full finite atom/ionic-axis/field-axis coefficient required')
        transform = slot_transport(grid, slot)
        r = transform.cartesian_rotation
        total += slot.weight*np.einsum('ij,ajk,lk->ail', r,
                    value[transform.atom_pullback], r)
        diagnostics.append(dict(spin=slot.ispin, mesh_index=list(slot.mesh_index),
                                source_ik=slot.source_ik, **transform.diagnostics()))
    return dict(mixed_ionic_field=total, slots=diagnostics,
                normalization='full spin/mesh mean; slot weight applied exactly once',
                geometry_policy=geometry_policy, h0_transport_supported=False)


class UPFPreparedVertices:
    """Full-row actions for DeltaV0, Delta_g and physical L at zero field.

    A runtime sector must bind its ACTUAL Hamiltonian/basis, not a transported
    H0 claim. No fixed contact or mixed term is inserted in electron H: those
    terms are prepared scalar coefficients, whose rho derivatives vanish.
    """
    def __init__(self, hamiltonian):
        self.hamiltonian = hamiltonian
        self.momenta = MaterialOpticalCurrent(hamiltonian).momenta

    def evaluate(self, frame, density_coefficients, displacement, shift):
        psi = np.asarray(frame)
        d = np.asarray(density_coefficients)
        p = self.momenta
        inputs = self.hamiltonian.background.inputs
        nat = inputs.reference.positions_bohr.shape[1]
        u = np.asarray(displacement).reshape(nat, 3)
        a = np.asarray(shift)
        if (psi.ndim != 2 or len(psi) != len(p) or d.shape != (psi.shape[1],)*2
                or a.shape != (3,) or not all(np.isfinite(x).all() for x in (psi,d,u,a))):
            raise ValueError('finite full frame,D,u and Cartesian a required')
        diagonal = np.einsum('pi,pi->p', psi.conj(), psi @ d)
        field_gradient = 2*p.T @ diagonal
        delta_gradient = np.zeros((nat, 3), complex)
        delta_energy = 0j
        action = 2*(p @ a)[:, None]*psi
        for species in self.hamiltonian.factors.species:
            tables = inputs.tables[species.name]
            centered = first_cartesian_jets(UPFRadial(tables), p)
            nh = species.active_d.shape[0]
            for slot, atom in enumerate(species.atom_indices):
                phase = np.exp(-1j*(p @ inputs.reference.positions_bohr[:, atom-1]))[:, None]
                b, *ba = [v*phase for v in centered]
                legacy = species.beta[:, slot*nh:(slot+1)*nh]
                def value_and_ionic(beta):
                    x = beta.conj().T @ psi
                    dx = species.active_d @ x
                    wd = dx @ d
                    value = np.vdot(x, wd)
                    g = np.empty(3, complex)
                    bu = -1j*(p @ u[atom-1])[:, None]*beta
                    xu = bu.conj().T @ psi
                    acted = beta @ dx + bu @ dx + beta @ (species.active_d @ xu)
                    for i in range(3):
                        xi = (-1j*p[:, i, None]*beta).conj().T @ psi
                        g[i] = np.vdot(xi, wd)+np.vdot(x, species.active_d @ xi @ d)
                    return value, g, acted, x, dx
                val, grad, acted, x, dx = value_and_ionic(b)
                oldval, oldgrad, oldacted, _, _ = value_and_ionic(legacy)
                delta_energy += val-oldval
                delta_gradient[atom-1] += grad-oldgrad
                action += acted-oldacted
                for j in range(3):
                    xj = ba[j].conj().T @ psi
                    dxj = species.active_d @ xj
                    field_gradient[j] += np.vdot(xj, dx @ d)+np.vdot(x, dxj @ d)
                    action += a[j]*(ba[j] @ dx+b @ dxj)
        return dict(action=action, delta_energy=delta_energy,
                    delta_gradient=delta_gradient.reshape(-1), field_gradient=field_gradient)
