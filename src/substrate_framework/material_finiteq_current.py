"""Rectangular physical current for one complex Bloch-Fourier vector potential.

Declared model: kinetic (p+A)^2 in Ry and the ordinary-NC scalar kernel
V_A(r,r')=exp[-i integral_straight_(r'->r) A.dl] V_NL(r,r').
A(r)=amplitude*a*exp(i q.r); amplitude is bohr^-1. prepare(a) returns
L_Ry=dH/damplitude (Ry bohr), and current methods return J_Ha=-L_Ry/2.
This is an explicit path choice, not unique all-electron transverse accuracy.
No local-potential variation, cell/spin/k weight, A25 offset, contact term,
photon mode, boundary condition, occupation or populated response is added.

The line integral uses a caller-controlled Gauss-Legendre order (not an error
certificate). Same-endpoint q=0 uses the existing analytic uniform vertex.
Only the full output columns and streamed endpoint/projector/column blocks
are allocated; there is no PW-square matrix. Inputs/factors remain borrowed.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .material_physical_operators import _require_physical, PhysicalFieldVertices
from .material_finiteq_projectors import placed_blocks, forward_block, adjoint_block
from .qe_bloch_finiteq import BlochFiniteQBridge
from .qe_bloch_finiteq_ionic import BlochFiniteQIonic


MODEL = 'ordinary-NC-straight-line-Peierls-physical-current-v1'


def _positive(value, name):
    if type(value) is not int or value < 1:
        raise ValueError('positive integer '+name+' required')
    return value


@dataclass(frozen=True)
class FiniteQCurrentAction:
    owner: object
    direction: np.ndarray

    def apply(self, vectors):
        """L_Ry[a] sourcePW x columns -> targetPW x columns."""
        return self.owner._apply(self.direction, vectors, adjoint=False)

    def adjoint_apply(self, vectors):
        """True L_Ry[a]^dagger target -> source, including conjugate(a)."""
        return self.owner._apply(self.direction, vectors, adjoint=True)

    def current_apply(self, vectors):
        """Unweighted electronic J_Ha[a]=-L_Ry[a]/2."""
        return -.5*self.apply(vectors)

    def current_adjoint_apply(self, vectors):
        return -.5*self.adjoint_apply(vectors)


class MaterialFiniteQCurrent:
    """Bind an actual checked edge and direct-physical endpoint Hamiltonians.

    block_size in 1..16 bounds column products. atom_block_size bounds placed
    projector columns by atom_block_size*max_species_nh; radial rows are already
    streamed by UPFRadial. Line nodes and species are streamed, never cached.
    Existing physical endpoint factor banks are authenticated, not corrected.
    """
    def __init__(self, edge, source_hamiltonian, target_hamiltonian, *,
                 line_order=12, block_size=2, atom_block_size=1):
        self.inputs = source_hamiltonian.background.inputs
        self.edge = edge
        self.bridge = BlochFiniteQBridge(self.inputs, edge, block_size=block_size)
        self.atom_block_size = _positive(atom_block_size, 'atom_block_size')
        self.line_order = _positive(line_order, 'line_order')
        self.source_hamiltonian = source_hamiltonian
        self.target_hamiltonian = target_hamiltonian
        if (target_hamiltonian.background.inputs is not self.inputs
                or source_hamiltonian.ispin != edge.ispin
                or target_hamiltonian.ispin != edge.ispin):
            raise ValueError('same actual inputs and explicit edge spin required')
        # Reuse the owning endpoint factor/UPF/atom/geometry validator, without
        # constructing ionic radial tables or allocating an ionic action.
        self.nat = self.inputs.reference.nat
        for h, basis in ((source_hamiltonian, edge.source_basis),
                         (target_hamiltonian, edge.target_basis)):
            _require_physical(h.factors)
            BlochFiniteQIonic._check_factors(self, h.factors, basis)
        self.ps = self._momenta(edge.source_basis)
        self.pt = self._momenta(edge.target_basis)
        self.species = source_hamiltonian.factors.species
        nodes, weights = leggauss(self.line_order)
        self.nodes, self.weights = (nodes+1)/2, weights/2
        self.uniform = (not np.any(edge.q_cart_bohr_inv)
                        and edge.source_basis.identity == edge.target_basis.identity)
        self.uniform_vertices = (PhysicalFieldVertices(source_hamiltonian,
            projector_block_size=block_size) if self.uniform else None)
        # An integer join avoids FFT cyclic aliases and any PW-square matrix.
        lookup = {tuple(m): i for i, m in enumerate(edge.target_basis.miller.T)}
        pairs = [(lookup[tuple(m-edge.reciprocal_fold)], j)
                 for j, m in enumerate(edge.source_basis.miller.T)
                 if tuple(m-edge.reciprocal_fold) in lookup]
        self.target_rows = np.asarray([i for i, _ in pairs], dtype=np.intp)
        self.source_rows = np.asarray([j for _, j in pairs], dtype=np.intp)

    @staticmethod
    def _momenta(basis):
        return (basis.reciprocal_bohr_inv @ basis.miller+basis.k_cart_bohr_inv[:, None]).T

    def prepare(self, polarization):
        """Complex-linear Cartesian polarization (3,); never normalize/project it."""
        a = np.array(polarization, dtype=complex, copy=True)
        if a.shape != (3,) or not np.isfinite(a).all():
            raise ValueError('finite complex Cartesian polarization (3,) required')
        a.setflags(write=False)
        return FiniteQCurrentAction(self, a)

    def _uniform_apply(self, direction, vectors, adjoint):
        a = direction.conj() if adjoint else direction
        result = self.uniform_vertices.prepare(a.real).apply(vectors)
        if np.any(a.imag):
            result += 1j*self.uniform_vertices.prepare(a.imag).apply(vectors)
        return result

    def _apply(self, direction, vectors, *, adjoint):
        rows = self.edge.target_dimension if adjoint else self.edge.source_dimension
        x = self.bridge._matrix(vectors, rows)
        if self.uniform:
            return self._uniform_apply(direction, x, adjoint)
        outrows = self.edge.source_dimension if adjoint else self.edge.target_dimension
        out = np.zeros((outrows, x.shape[1]), complex)
        tr, sr = self.target_rows, self.source_rows
        kinetic = (self.pt[tr]+self.ps[sr]) @ direction
        incoming, outgoing = (tr, sr) if adjoint else (sr, tr)
        if adjoint:
            kinetic = kinetic.conj()
        size = self.bridge.block_size
        for start in range(0, x.shape[1], size):
            sl = slice(start, start+size)
            out[outgoing, sl] = kinetic[:, None]*x[incoming, sl]
        action = adjoint_block if adjoint else forward_block
        for t, weight in zip(self.nodes, self.weights):
            for species in self.species:
                for banks in placed_blocks(self.inputs, species, self.pt, self.ps,
                        self.edge.q_cart_bohr_inv, direction, t, self.atom_block_size):
                    for start in range(0, x.shape[1], size):
                        sl = slice(start, start+size)
                        out[:, sl] += weight*action(*banks, x[:, sl])
        return out
