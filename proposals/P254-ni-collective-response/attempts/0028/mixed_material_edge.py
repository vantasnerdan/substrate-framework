"""Endpoint-local physical banks and complex-linear mixed material actions.

Only actual physical H0, finite-q ionic/local vertices and Gamma L are used.
The mean action is the completed moving history's PhysicalFrozenAction, not
an eigenvalue, stationary frame or a reconstructed mean density.
"""

import numpy as np

from substrate_framework.material_physical_operators import (
    physical_hamiltonian, physical_coupling, physical_finiteq,
    PhysicalFieldVertices, PhysicalFrozenAction)



class ComplexFieldAction:
    """Complex-linear extension of the real physical Cartesian L vertex."""
    def __init__(self, vertex, direction):
        a = np.asarray(direction)
        if a.shape != (3,) or not np.isfinite(a).all():
            raise ValueError('finite complex Cartesian field direction required')
        self.real = vertex.prepare(a.real) if np.any(a.real != 0) else None
        self.imag = vertex.prepare(a.imag) if np.any(a.imag != 0) else None

    def add(self, out, vectors):
        if self.real is not None:
            out += self.real.apply(vectors)
        if self.imag is not None:
            out += 1j*self.imag.apply(vectors)
        return out


class NativeEdgeOperators:
    """One source plus its incident endpoint banks, never derivative payloads."""
    def __init__(self, owner, index, stage):
        self.owner, self.index, self.stage = owner, index, stage
        self.h0, self.gamma, self.vertices = {}, {}, {}
        self.couplings, self.mean_actions, self.perturbations = {}, {}, {}

    def bare(self, index):
        if index not in self.h0:
            self.owner.budget.check()
            basis = self.owner.baseline.bases[index]
            self.h0[index] = physical_hamiltonian(self.owner.case.background,
                basis, basis.slot.ispin, self.owner.block_size)
        return self.h0[index]

    def gamma_coupling(self, index):
        if index not in self.gamma:
            h0 = self.bare(index)
            self.gamma[index] = physical_coupling(self.owner.case.background.inputs,
                self.owner.baseline.bases[index], h0.factors, self.owner.block_size)
        return self.gamma[index]

    def vertex(self, index):
        if index not in self.vertices:
            self.vertices[index] = PhysicalFieldVertices(self.bare(index),
                projector_block_size=self.owner.block_size)
        return self.vertices[index]

    def coupling(self, direction, source):
        key = direction, source
        if key not in self.couplings:
            link = self.owner.links[direction][source]
            self.couplings[key] = physical_finiteq(self.owner.case.background.inputs,
                link.coupling.edge, self.bare(source).factors,
                self.bare(link.target).factors, self.owner.block_size)
        return self.couplings[key]


    def mean(self, index, vectors):
        if not np.any(vectors):
            return np.zeros_like(vectors, dtype=np.complex128)
        if index not in self.mean_actions:
            spin = self.owner.baseline.bases[index].slot.ispin
            self.mean_actions[index] = PhysicalFrozenAction(self.bare(index),
                self.gamma_coupling(index), potential=self.stage.potential[:, spin-1],
                direction=self.stage.u, shift=self.stage.shift,
                fft_block_size=self.owner.block_size,
                projector_block_size=self.owner.block_size)
        self.owner.budget.action(vectors)
        return self.mean_actions[index].apply(vectors)

    def perturbation(self, direction, source, vectors, *, potential, small):
        link = self.owner.links[direction][source]
        shape = (self.owner.baseline.bases[link.target].miller.shape[1], vectors.shape[1])
        gamma = self.owner.histories[direction].response.gamma
        if (not np.any(vectors) or
                (not np.any(potential) and not np.any(small.u)
                 and (not gamma or not np.any(small.coordinate)))):
            return np.zeros(shape, dtype=np.complex128)
        key = direction, source
        if key not in self.perturbations:
            coupling = self.coupling(direction, source)
            ionic = coupling.prepare(np.asarray(small.u).reshape(-1, 3))
            field = None
            if self.owner.histories[direction].response.gamma:
                cartesian = np.asarray(small.coordinate) @ self.owner.case.field.polarizations
                field = ComplexFieldAction(self.vertex(source), cartesian)
            self.perturbations[key] = coupling, ionic, field
        coupling, ionic, field = self.perturbations[key]
        self.owner.budget.action(vectors)
        out = ionic.apply(vectors)
        out += coupling.bridge.local_apply(potential[:, coupling.edge.ispin-1], vectors)
        if field is not None:
            field.add(out, vectors)
        return out
