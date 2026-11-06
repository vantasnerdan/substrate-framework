"""Full-PW H0 flow and directly assembled moving-mean residual actions.

These operations retain the existing physical H0, original PW bases and raw D.
No Hfull-H0 subtraction, saved-band projection or mean extrapolation is used.
The supplied history/budget owns the actual authenticated physical inputs and
wall/action/memory limits. Exponential receipts are floating error estimates,
not interval certificates or a material long-window accuracy claim.
"""
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from joint_column_mean import ColumnLimits
from moving_exponential_h0 import PhysicalH0Exponential
from substrate_framework.material_joint_response import RealBlochEdge
from substrate_framework.material_physical_operators import (
    physical_hamiltonian, physical_finiteq, PhysicalFieldVertices)
from substrate_framework.streamed_joint_response import StreamedJointResponse
from moving_exponential_response import FeedbackContext, SmallState


class DirectMeanResidual:
    """Apply Vmean + g(Umean) + L(Amean), with independent adjoints."""
    def __init__(self, coupling, potential, direction, field):
        self.bridge = coupling.bridge
        self.potential = potential
        self.ionic = (coupling.prepare(np.asarray(direction).reshape(-1, 3))
                      if np.any(direction) else None)
        self.field = field

    def apply(self, vectors):
        if not np.any(vectors):
            return vectors
        out = self.bridge.local_apply(self.potential, vectors)
        if self.ionic is not None:
            out += self.ionic.apply(vectors)
        if self.field is not None:
            out += self.field.apply(vectors)
        return out

    def adjoint_apply(self, vectors):
        if not np.any(vectors):
            return vectors
        out = self.bridge.local_adjoint_apply(self.potential, vectors)
        if self.ionic is not None:
            out += self.ionic.adjoint_apply(vectors)
        if self.field is not None:
            out += self.field.adjoint_apply(vectors)
        return out


class MaterialH0Flow:
    """Source-local exact-H0 exponential approximation under explicit limits."""
    def __init__(self, history, *, limits, receipt_directory):
        if not isinstance(limits, ColumnLimits):
            raise TypeError('explicit original-column numerical limits required')
        self.history, self.limits = history, limits
        self.receipt_directory = Path(receipt_directory)
        self.receipt_directory.mkdir(parents=True, exist_ok=True)
        self.calls = self.actions = self.columns = 0

    def _apply(self, index, vectors, interval, *, adjoint):
        self.history.budget.check()
        self.calls += 1
        if interval == 0 or not np.any(vectors):
            return vectors
        basis = self.history.bases[index]
        h0 = physical_hamiltonian(self.history.background, basis,
                                   basis.slot.ispin, self.history.block_size)
        operation = PhysicalH0Exponential(h0, budget=self.history.budget,
            limits=self.limits, identity=self.history.final.identity+':'+basis.identity,
            receipt_directory=self.receipt_directory, fft_block_size=self.history.block_size)
        try:
            return operation.propagate(vectors, interval, adjoint=adjoint)
        finally:
            self.actions += operation.actions
            self.columns += operation.columns

    def propagate(self, index, vectors, interval):
        return self._apply(index, vectors, interval, adjoint=False)

    def propagate_adjoint(self, index, vectors, interval):
        return self._apply(index, vectors, interval, adjoint=True)

    def free_pair(self, edge_index, pair, interval):
        edge = self.history.response.edges[edge_index]
        left, right = pair
        return (self.propagate(edge.target, left, interval),
                self.propagate_adjoint(edge.source, right.conj().T, interval).conj().T)


@contextmanager
def _at_time(history, time):
    """Right-hand stage at an internal boundary, independent of replay order."""
    active = history.active
    try:
        history.clear_interval()
        if time == history.final.time:
            history.active = len(history.records)-1
        else:
            history.set_interval(time, 0.)
        yield
    finally:
        history.active = active


@contextmanager
def material_residual_edge(history, index, time):
    """Open one actual moving edge, without relying on mutable active interval."""
    history.budget.check()
    with _at_time(history, time):
        _, _, stage = history._stage(time)
        edge = history.response.edges[index]
        source, target = edge.source, edge.target
        hs = physical_hamiltonian(history.background, history.bases[source],
                                  history.bases[source].slot.ispin, history.block_size)
        ht = hs if source == target else physical_hamiltonian(
            history.background, history.bases[target],
            history.bases[target].slot.ispin, history.block_size)

        def residual(basis, h0):
            geometry = RealBlochEdge.create(history.background.inputs.reference,
                                            basis, basis, np.zeros(3))
            coupling = physical_finiteq(history.background.inputs, geometry,
                                        h0.factors, h0.factors, history.block_size)
            field = (PhysicalFieldVertices(h0).prepare(stage.shift)
                     if np.any(stage.shift) else None)
            return DirectMeanResidual(coupling, stage.potential[:, basis.slot.ispin-1],
                                      stage.u, field)

        rs = residual(history.bases[source], hs)
        rt = rs if source == target else residual(history.bases[target], ht)
        coupling = physical_finiteq(history.background.inputs, edge.coupling.edge,
                                    hs.factors, ht.factors, history.block_size)
        ps = history.frame(source, time)
        pt = ps if source == target else history.frame(target, time)
        yield ps, pt, rs, rt, coupling, PhysicalFieldVertices(hs) if history.response.gamma else None


class MaterialMidpointModel:
    """Complete source-local forward feedback for the existing physical scalar.

    Reuse the settled streamed canonical map without allocating RK4 scratch.
    Feedback completes every source before applying K; neither moments nor
    forward drives multiply an already embedded raw D a second time.
    """
    dynamics = 'screened-joint'

    def __init__(self, history, flow):
        if flow.history is not history:
            raise ValueError('H0 flow and moving response must share physical history')
        self.history, self.flow = history, flow
        self.response = self.r = history.response
        self.free_pair = flow.free_pair

    def prepare_interval(self, time, dt):
        active = self.history.active
        try:
            self.history.clear_interval()
            self.history.set_interval(time, dt)
        finally:
            self.history.active = active

    def context(self, time, state):
        r = self.r
        with _at_time(self.history, time):
            self.history._index(time)
        density = np.zeros(r.density_shape, complex)
        gradient, ell = np.zeros(r.ndof, complex), np.zeros(3, complex)
        for index, edge in enumerate(r.edges):
            self.history.budget.check()
            left, right = state.edge(index)
            if not np.any(left) and not np.any(right):
                continue
            weight = r.sectors[edge.source].weight
            with material_residual_edge(self.history, index, time) as data:
                ps, pt, _, _, coupling, vertices = data
                n1, g1 = coupling.moments(left, ps, weight=weight)
                n2, g2 = coupling.moments(pt, right.conj().T, weight=weight)
                density += n1+n2
                gradient += g1+g2
                if r.gamma:
                    for j, axis in enumerate(np.eye(3)):
                        action = vertices.prepare(axis)
                        ell[j] += weight*(np.vdot(action.apply(ps), left)
                            +np.einsum('ij,ji->', right, action.adjoint_apply(pt)))
        potential = np.asarray(r.kernel(density, q_bohr_inv=r.q_bohr_inv), complex)
        if potential.shape != r.density_shape or not np.isfinite(potential).all():
            raise ValueError('complete finite physical Kq density potential required')
        small = SmallState(*(np.zeros_like(getattr(state, name), dtype=complex)
                             for name in ('u', 'p', 'coordinate', 'momentum')))
        StreamedJointResponse.small_rhs(self, state, small, gradient, ell)
        return FeedbackContext(float(time), potential, small)

    def edge_rhs(self, context, index, pair, small):
        r = self.r
        edge = r.edges[index]
        source, target = r.sectors[edge.source], r.sectors[edge.target]
        left, right = pair
        if (left.shape != (target.shape[0], source.shape[1])
                or right.shape != (target.shape[1], source.shape[0])):
            raise ValueError('original rectangular independent full-PW legs required')
        with material_residual_edge(self.history, index, context.time) as data:
            ps, pt, residual_source, residual_target, coupling, vertices = data
            ionic = coupling.prepare(small.u.reshape(-1, 3))
            potential = context.potential[:, source.density_component]
            drive = coupling.bridge.local_apply(potential, ps)+ionic.apply(ps)
            adjoint = coupling.bridge.local_adjoint_apply(potential, pt)+ionic.adjoint_apply(pt)
            if r.gamma:
                field = r.sqrt_n*r.polarizations.T @ small.coordinate
                for coefficient, axis in zip(field, np.eye(3)):
                    if coefficient != 0:
                        action = vertices.prepare(axis)
                        drive += coefficient*action.apply(ps)
                        adjoint += coefficient.conjugate()*action.adjoint_apply(pt)
            return (-1j*(residual_target.apply(left)+drive @ source.d),
                    1j*(residual_source.adjoint_apply(right.conj().T).conj().T
                        +target.d @ adjoint.conj().T))
