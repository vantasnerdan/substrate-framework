"""Bare finite-q ordinary-NC derivative and its positive-gradient transpose.

QE compute_dvloc/phq_init/dvqpsi_us_only: local -i(q+G).u V(q+G),
plus BOTH i Bt D Bs^dagger ps.u and -i pt.u Bt D Bs^dagger. Placed
projectors already contain the atomic phase; no additional q phase is added.
This is unscreened, Ry/bohr, all actual atoms, full-PW on both legs.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_bloch_finiteq import BlochFiniteQBridge
from .qe_material_hamiltonian import build_pw_factors_for_basis
from .qe_nc_radial import build_radial_tables


def finiteq_radial_tables(inputs, edge):
    """Extend only numerical q tables to every requested |G+q| and |k+G|.

    Original UPF objects/bytes, cell volume and species D remain unchanged.
    The caller may reuse the returned mapping over edges with a smaller range.
    """
    ref = inputs.reference
    gq = ref.reciprocal_bohr_inv @ inputs.density.miller + edge.q_cart_bohr_inv[:, None]
    required = float(np.max(np.linalg.norm(gq, axis=0), initial=0.))
    for basis in (edge.source_basis, edge.target_basis):
        p = basis.reciprocal_bohr_inv @ basis.miller + basis.k_cart_bohr_inv[:, None]
        required = max(required, float(np.max(np.linalg.norm(p, axis=0), initial=0.)))
    tables = {}
    for name, table in inputs.tables.items():
        tables[name] = (table if table.qmax_bohr_inv >= required else
                        build_radial_tables(table.upf, table.omega_bohr3, required))
    return tables


@dataclass(frozen=True)
class FiniteQDirectionalAction:
    owner: object
    direction: np.ndarray
    field: np.ndarray

    def apply(self, vectors):
        """Complex-linear forward derivative; no conjugation of direction."""
        return self.owner._apply(self.direction, self.field, vectors, adjoint=False)

    def adjoint_apply(self, vectors):
        """Independent Hilbert adjoint using the SAME forward local support."""
        return self.owner._apply(self.direction, self.field, vectors, adjoint=True)


class BlochFiniteQIonic:
    """All-coordinate matrix-free coupling; factors and inputs must stay immutable.

    moments returns X_ts density and grad_a=weight Tr(g_ts[e_a]^dagger X_ts).
    Thus sum_a conj(u_a) grad_a = weight Tr(g_ts[u]^dagger X_ts).
    Equivalently sum_a u_a grad_a pairs X_ts with the opposite-edge derivative
    g_st[u] whose local support is the exact adjoint of this forward support.
    This is a positive gradient, NOT its negative (force). No implicit real part.
    """

    def __init__(self, inputs, edge, *, source_factors=None, target_factors=None,
                 tables=None, block_size=8):
        self.bridge = BlochFiniteQBridge(inputs, edge, block_size=block_size)
        self.inputs, self.edge = inputs, edge
        self.source_factors = (build_pw_factors_for_basis(inputs, edge.source_basis)
                               if source_factors is None else source_factors)
        self.target_factors = (build_pw_factors_for_basis(inputs, edge.target_basis)
                               if target_factors is None else target_factors)
        self.nat = inputs.reference.nat
        self.ndof = 3*self.nat
        for factors, basis in ((self.source_factors, edge.source_basis),
                               (self.target_factors, edge.target_basis)):
            self._check_factors(factors, basis)
        source = {s.name: s for s in self.source_factors.species}
        self.species = tuple((source[t.name], t) for t in self.target_factors.species)
        self.ps = (edge.source_basis.reciprocal_bohr_inv @ edge.source_basis.miller
                   + edge.source_basis.k_cart_bohr_inv[:, None]).T
        self.pt = (edge.target_basis.reciprocal_bohr_inv @ edge.target_basis.miller
                   + edge.target_basis.k_cart_bohr_inv[:, None]).T
        self.gq = self.bridge.density_fft.wavevectors.T
        self.tables = finiteq_radial_tables(inputs, edge) if tables is None else tables
        moduli = np.linalg.norm(self.gq, axis=1)
        self.vloc = {}
        for name in dict.fromkeys(inputs.reference.atom_names):
            table = self.tables[name]
            original = inputs.tables[name]
            if (table.upf is not original.upf or table.omega_bohr3 != original.omega_bohr3
                    or table.qmax_bohr_inv < float(np.max(moduli, initial=0.))):
                raise ValueError('local q tables must retain original UPF/volume and cover |G+q|')
            self.vloc[name] = table.vloc(moduli)

    def _check_factors(self, factors, basis):
        ref = self.inputs.reference
        if (not factors.sources or factors.sources[-1].get('identity') != basis.identity
                or factors.sources[:-1] != self.inputs.sources
                or not np.array_equal(factors.fft.miller, basis.miller)
                or not np.array_equal(factors.k_bohr_inv, basis.k_cart_bohr_inv)
                or not np.array_equal(factors.fft.reciprocal_bohr_inv, basis.reciprocal_bohr_inv)
                or factors.fft.physical != self.bridge.source_fft.physical
                or factors.fft.omega != self.bridge.source_fft.omega):
            raise ValueError('NC factors must match the bound endpoint/geometry/physical grid')
        covered, names = [], set()
        for species in factors.species:
            if species.name in names or species.name not in self.inputs.tables:
                raise ValueError('one projector factor group per original species required')
            names.add(species.name)
            table = self.inputs.tables[species.name]
            atoms = tuple(a+1 for a, name in enumerate(ref.atom_names) if name == species.name)
            nh = table.upf.nh
            if (not atoms or species.atom_indices != atoms
                    or species.beta.shape != (basis.miller.shape[1], len(atoms)*nh)
                    or table.omega_bohr3 != factors.fft.omega
                    or not np.array_equal(species.active_d, table.upf.expanded_d)
                    or not np.array_equal(species.active_d, species.active_d.conj().T)):
                raise ValueError('NC factors require original atom order and full Hermitian species D')
            covered.extend(atoms)
        if sorted(covered) != list(range(1, self.nat+1)):
            raise ValueError('all actual atoms must be present exactly once')

    def prepare(self, direction):
        """Dimensionless complex direction[nat,3]; amplitude is in bohr."""
        u = np.array(direction, dtype=np.complex128, copy=True)
        if u.shape != (self.nat, 3) or not np.isfinite(u).all():
            raise ValueError('finite full atom-order direction[nat,3] required')
        local = np.zeros(self.gq.shape[0], dtype=np.complex128)
        ref = self.inputs.reference
        for atom in np.flatnonzero(np.any(u != 0, axis=1)):
            local += (-1j*(self.gq @ u[atom])*self.vloc[ref.atom_names[atom]]
                      * np.exp(-1j*(self.gq @ ref.positions_bohr[:, atom])))
        field = self.bridge.density_fft.from_g(local)
        u.setflags(write=False)
        field.setflags(write=False)
        return FiniteQDirectionalAction(self, u, field)

    def _apply(self, direction, field, vectors, *, adjoint):
        bridge = self.bridge
        rows = self.edge.target_dimension if adjoint else self.edge.source_dimension
        vectors = bridge._matrix(vectors, rows)
        outrows = self.edge.source_dimension if adjoint else self.edge.target_dimension
        out = np.empty((outrows, vectors.shape[1]), dtype=np.complex128)
        u = direction.conj() if adjoint else direction
        pin, pout = (self.pt, self.ps) if adjoint else (self.ps, self.pt)
        local_action = bridge.local_adjoint_apply if adjoint else bridge.local_apply
        for start in range(0, vectors.shape[1], bridge.block_size):
            sl = slice(start, start+bridge.block_size)
            block = vectors[:, sl]
            result = local_action(field, block)
            for source, target in self.species:
                incoming, outgoing = (target, source) if adjoint else (source, target)
                atoms = np.asarray(incoming.atom_indices, dtype=int)-1
                nh, width = incoming.active_d.shape[0], block.shape[1]
                shape = (len(atoms), nh, width)
                d = incoming.active_d.conj().T if adjoint else incoming.active_d
                base = zgemm(1., incoming.beta, block, trans_a=2).reshape(shape)
                dbase = np.einsum('ij,ajb->aib', d, base, optimize=True)
                for axis in range(3):
                    ua = u[atoms, axis, None, None]
                    if not np.any(ua):
                        continue
                    derivative = zgemm(1., incoming.beta, pin[:, axis, None]*block,
                                       trans_a=2).reshape(shape)
                    weighted = np.einsum('ij,ajb->aib', d, derivative, optimize=True)*ua
                    result += 1j*zgemm(1., outgoing.beta, weighted.reshape(-1, width))
                    weighted = dbase*ua
                    result -= 1j*pout[:, axis, None]*zgemm(
                        1., outgoing.beta, weighted.reshape(-1, width))
            out[:, sl] = result
        return out

    def _nonlocal_gradient_block(self, left, right):
        """Tr(g_ts[e_a]^dagger L_t R_s^dagger), all species/axes."""
        result = np.zeros((self.nat, 3), dtype=np.complex128)
        width = left.shape[1]
        for source, target in self.species:
            atoms = np.asarray(source.atom_indices, dtype=int)-1
            shape = (len(atoms), source.active_d.shape[0], width)
            bl = zgemm(1., target.beta, left, trans_a=2).reshape(shape)
            br = zgemm(1., source.beta, right, trans_a=2).reshape(shape)
            d = source.active_d.conj().T
            dbl = np.einsum('ij,ajb->aib', d, bl, optimize=True)
            for axis in range(3):
                dl = zgemm(1., target.beta, self.pt[:, axis, None]*left, trans_a=2).reshape(shape)
                dr = zgemm(1., source.beta, self.ps[:, axis, None]*right, trans_a=2).reshape(shape)
                ddl = np.einsum('ij,ajb->aib', d, dl, optimize=True)
                result[atoms, axis] = 1j*(np.einsum('aib,aib->a', br.conj(), ddl)
                                                - np.einsum('aib,aib->a', dr.conj(), dbl))
        return result

    def moments(self, left, right, coefficient=None, *, weight):
        """Return (density[nr,2], positive_gradient[3*nat]) for weight L_t C R_s†.

        No PW-square matrix, coordinate derivative matrix, or global species-D
        matrix is formed. Both transient factor legs are bounded column blocks.
        """
        bridge = self.bridge
        weight = bridge._weight(weight)
        density = np.zeros((bridge.source_fft.nnr, 2), dtype=np.complex128)
        channel = density[:, self.edge.ispin-1]
        gradient = np.zeros((self.nat, 3), dtype=np.complex128)
        for lblock, rblock in bridge.iter_dyad_blocks(left, right, coefficient):
            lwave = bridge.target_fft.from_g(lblock)
            rwave = bridge.source_fft.from_g(rblock)
            channel += np.sum(lwave*rwave.conj(), axis=1)
            gradient += self._nonlocal_gradient_block(lblock, rblock)
        channel *= weight*bridge.multiplication_phase.conj()/bridge.source_fft.omega
        gradient *= weight
        # Exact support-pinned reverse transpose: Omega sum_G n_q(G) conj(v_q,a(G)).
        sampled = bridge.density_fft.to_g(channel)
        ref = self.inputs.reference
        for atom, name in enumerate(ref.atom_names):
            weighted = (sampled*self.vloc[name].conj()
                        * np.exp(1j*(self.gq @ ref.positions_bohr[:, atom])))
            gradient[atom] += 1j*bridge.source_fft.omega*(weighted @ self.gq)
        return density, gradient.reshape(-1)
