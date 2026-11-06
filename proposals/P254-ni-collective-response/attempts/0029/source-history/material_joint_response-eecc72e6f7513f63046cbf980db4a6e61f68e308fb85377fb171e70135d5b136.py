"""Explicit physical finite-q and completed joint-mean history adapters.

RealBlochEdge is a real extension of the checked geometry contract: both basis
kinds validate their own source/scope and every mesh/spin/fold invariant is
checked. It does not pretend RealBlochBasis is BlochBasis, conjugate a positive
-time state, or change the saved physical factor provenance.
"""
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .qe_bloch_edge import BlochEdge
from .qe_bloch_frame import BlochBasis
from .qe_real_bloch_grid import RealBlochBasis
from .qe_bloch_grid import qe_fold, readonly
from .material_physical_operators import (physical_hamiltonian, physical_coupling,
    physical_finiteq, PhysicalFieldVertices, PhysicalFrozenAction)
from .joint_bloch_response import JointBlochResponse, JointSector, JointStage, CartesianFieldActions
from .streamed_joint_state import DiskJointMean
from .streamed_joint_reduction import JointStage as MeanStage
from .streamed_material_state import authenticate, DiskMean, file_binding
from .streamed_material_evolution import exponential_action
from .qe_bloch_finiteq import BlochFiniteQBridge
from .paired_finiteq_moments import paired_moments_and_drive as paired_matter_terms
from .paired_field_moments import (CartesianPhysicalFieldAxes,
    paired_field_moments_and_drive as paired_field_terms)


def _detached_basis(basis):
    """Retain metadata without retaining its entire native WFC byte buffer."""
    metadata = {}
    for name in ('reciprocal_bohr_inv', 'source_miller'):
        value = getattr(basis, name).copy(order='K')
        value.setflags(write=False)
        metadata[name] = value
    return replace(basis, **metadata)


class RealBlochEdge(BlochEdge):
    @classmethod
    def create(cls, reference, source_basis, target_basis, q_fractional):
        for basis in (source_basis, target_basis):
            if not isinstance(basis, (BlochBasis, RealBlochBasis)):
                raise TypeError('independently admitted ordinary or real-collinear basis required')
            basis.validate(reference)
        source, target = source_basis.slot, target_basis.slot
        mesh = tuple(reference.monkhorst_pack['nk%d' % i] for i in (1, 2, 3))
        shift = np.array([reference.monkhorst_pack['k%d' % i] for i in (1, 2, 3)])
        if min(mesh) < 1 or not np.isin(shift, (0, 1)).all() or source.ispin != target.ispin or source.ispin not in (1, 2):
            raise ValueError('complete original same-spin mesh required')
        for slot in (source, target):
            index = np.asarray(slot.mesh_index)
            if (index.shape != (3,) or index.dtype.kind not in 'iu' or np.any(index < 0) or np.any(index >= mesh)
                    or not np.allclose(slot.k_fractional, qe_fold((index+shift/2)/mesh), atol=2e-10, rtol=0)):
                raise ValueError('basis slot does not match original declared mesh')
        q = np.asarray(q_fractional, float)
        if q.shape != (3,) or not np.isfinite(q).all():
            raise ValueError('finite reciprocal fractional q required')
        step = np.rint(q*np.asarray(mesh))
        if not np.allclose(q*np.asarray(mesh), step, atol=2e-10, rtol=0):
            raise ValueError('commensurate q required')
        index = tuple(int(v) for v in step.astype(np.int64) % mesh)
        if not np.array_equal((np.asarray(source.mesh_index)+index) % mesh, target.mesh_index):
            raise ValueError('target is not source+q')
        delta = target.k_fractional-source.k_fractional-q
        fold = np.rint(delta).astype(np.int64)
        if not np.allclose(delta, fold, atol=2e-10, rtol=0):
            raise ValueError('noninteger reciprocal fold')
        return cls(source_basis, target_basis, readonly(q), readonly(reference.reciprocal_bohr_inv @ q),
                   readonly(fold), index, mesh, source.ispin)


@dataclass(frozen=True)
class ResponseLink:
    source: int
    target: int
    coupling: object


class LazyDirection:
    def __init__(self, coupling, direction):
        self.coupling, self.direction = coupling, np.array(direction, complex)

    def apply(self, vectors):
        with self.coupling.open() as coupling:
            return coupling.prepare(self.direction).apply(vectors)

    def adjoint_apply(self, vectors):
        with self.coupling.open() as coupling:
            return coupling.prepare(self.direction).adjoint_apply(vectors)


class LazyCoupling:
    def __init__(self, owner, edge):
        self.owner, self.edge, self.inputs = owner, edge, owner.background.inputs

    @contextmanager
    def open(self):
        a = physical_hamiltonian(self.owner.background, self.edge.source_basis, self.edge.ispin, self.owner.block_size)
        b = a if self.edge.source_basis is self.edge.target_basis else physical_hamiltonian(
            self.owner.background, self.edge.target_basis, self.edge.ispin, self.owner.block_size)
        yield physical_finiteq(self.inputs, self.edge, a.factors, b.factors, self.owner.block_size)

    def prepare(self, direction):
        return LazyDirection(self, direction)


class LazyFieldAxis:
    def __init__(self, owner, index, axis):
        self.owner, self.index, self.axis = owner, index, axis

    def apply(self, vectors):
        h = physical_hamiltonian(self.owner.background, self.owner.bases[self.index],
                                 self.owner.bases[self.index].slot.ispin, self.owner.block_size)
        return PhysicalFieldVertices(h).prepare(np.eye(3)[self.axis]).apply(vectors)

    def adjoint_apply(self, vectors):
        return self.apply(vectors)


class InitialFrames:
    def __init__(self, state):
        self.state = state

    def __len__(self):
        return len(self.state.slots)

    def __getitem__(self, index):
        return self.state.arrays(self.state.slots[index])[0]


class MaterialJointHistory:
    def __init__(self, final_binding, *, model_identity, background, grid, ordered,
                 slot_identity, budget, block_size):
        self.final = DiskJointMean(authenticate(final_binding), expected_identity=model_identity)
        if not self.final.history or self.final.time <= 0:
            raise ValueError('completed actually moving physical joint history required')
        self.background, self.budget, self.block_size = background, budget, block_size
        self.records, self.priors = self.final.history, []
        self.cached_stage = self.cached_index = None
        last = 0.
        expected = tuple(slot_identity(s) for s in ordered)
        if tuple(s.identity for s in self.final.slots) != expected:
            raise ValueError('final boundary explicit slot coverage differs')
        authenticated = set()
        for record in self.records:
            if record['kind'] != 'actual-frozen-midpoint':
                raise ValueError('only completed actual physical midpoint histories admitted')
            prior = DiskJointMean(authenticate(record['prior_boundary_binding']), expected_identity=model_identity, verify=False)
            dt = record['dt_ry_inverse']
            if (prior.time != last or record['start_time_ry_inverse'] != last or dt <= 0
                    or record['time_ry_inverse'] != last+dt/2
                    or tuple(s.identity for s in prior.slots) != expected
                    or prior.origin != self.final.origin):
                raise ValueError('broken completed mean history/source chain')
            if self.priors:
                if prior.history != self.records[:len(self.priors)]:
                    raise ValueError('history prefix differs from completed prior boundary')
            for a, b in zip(prior.slots, self.final.slots):
                if a.metadata() != b.metadata() or a.coefficient != b.coefficient:
                    raise ValueError('original raw basis/gauge/D changed along history')
                for binding in (a.frame, a.coefficient):
                    key = (binding['path'], binding['sha256'])
                    if key not in authenticated:
                        budget.check()
                        authenticate(binding)
                        authenticated.add(key)
            self.priors.append(prior)
            last += dt
        if last != self.final.time or self.priors[0].time != 0:
            raise ValueError('history does not cover original preparation through final boundary')
        self.initial = self.priors[0]
        original = DiskMean(authenticate(self.initial.origin), verify=False)
        if (original.time != 0 or len(original.slots) != len(self.initial.slots)
                or any(a.metadata() != b.metadata() or a.frame != b.frame or a.coefficient != b.coefficient
                       for a, b in zip(original.slots, self.initial.slots))):
            raise ValueError('joint t0 is not the literal original electronic preparation')
        self.bases = []
        for slot, record in zip(ordered, self.initial.slots):
            frame = grid.load_frame(slot)
            if (frame.basis.identity != record.basis_identity or frame.gauge_identity != record.gauge_identity
                    or frame.coefficients.shape != record.shape or slot.source_binding['sha256'] != record.source_identity):
                raise ValueError('original basis/gauge/frame source differs')
            self.bases.append(_detached_basis(frame.basis))
            del frame
        self.active = None
        self.frame_work = None
        self.frame_records = {}

    def stage_identity(self, time):
        index = self._index(time)
        saved = self._saved_frames(time)
        boundary = self.priors[index] if saved is None else saved
        return dict(time=float(time), interval=index, stage=self.records[index],
                    baseline_boundary=file_binding(boundary.path), model=self.final.identity,
                    origin=self.initial.origin, block_size=self.block_size)

    def _stage(self, time):
        index = self._index(time)
        if self.cached_index != index:
            self.cached_stage = MeanStage.restore(self.records[index])
            self.cached_index = index
        return index, self.priors[index], self.cached_stage

    def _saved_frames(self, time):
        index, prior, _ = self._stage(time)
        endpoint = self.priors[index+1] if index+1 < len(self.priors) else self.final
        return prior if time == prior.time else endpoint if time == endpoint.time else None

    def prepare_time(self, time):
        """One actual exponential per slot/time, not two per incident edge/pass."""
        if self.frame_work is None or self._saved_frames(time) is not None:
            return
        identity = self.stage_identity(time)
        records = self.frame_work.execute('actual-baseline-frames', identity, len(self.bases),
                                          lambda index: dict(frame=self._frame(index, time)))
        self.frame_records[(self._index(time), float(time))] = records

    def _frame(self, index, time):
        _, prior, stage = self._stage(time)
        saved = self._saved_frames(time)
        if saved is not None:
            return saved.arrays(saved.slots[index])[0]
        h0 = physical_hamiltonian(self.background, self.bases[index],
                                  self.bases[index].slot.ispin, self.block_size)
        coupling = physical_coupling(self.background.inputs, self.bases[index], h0.factors, self.block_size)
        action = PhysicalFrozenAction(h0, coupling, potential=stage.potential[:, self.bases[index].slot.ispin-1],
                                      direction=stage.u, shift=stage.shift)
        return exponential_action(self.counted(action.apply), prior.arrays(prior.slots[index])[0], time-prior.time)

    def frame(self, index, time):
        saved = self._saved_frames(time)
        if saved is not None:
            return saved.arrays(saved.slots[index])[0]
        if self.frame_work is None:
            return self._frame(index, time)
        records = self.frame_records.get((self._index(time), float(time)))
        if records is None:
            raise ValueError('complete authenticated per-time baseline frame barrier required')
        binding = records[index]['arrays']['frame']
        self.frame_work.verify(binding)
        value = np.load(binding['path'], mmap_mode='r', allow_pickle=False)
        if value.shape != self.initial.slots[index].shape or value.dtype != np.complex128:
            raise ValueError('cached actual baseline frame dimensions differ')
        return value

    @contextmanager
    def feedback(self, index, time):
        """Local potential only: no H0, ionic g, or projector factor rebuild."""
        link = self.response.edges[index]
        bridge = BlochFiniteQBridge(self.background.inputs, link.coupling.edge, block_size=self.block_size)
        yield self.frame(link.source, time), self.frame(link.target, time), bridge

    def clear_interval(self):
        self.active = None

    def set_interval(self, time, dt):
        for i, record in enumerate(self.records):
            lo, hi = record['start_time_ry_inverse'], record['start_time_ry_inverse']+record['dt_ry_inverse']
            if lo <= time < hi and time+dt <= hi+4*np.spacing(hi):
                self.active = i
                return
        raise ValueError('response step must not cross completed mean stage boundary')

    def _index(self, time):
        if not np.isfinite(time) or not 0 <= time <= self.final.time:
            raise ValueError('time not covered by completed moving history')
        if self.active is not None:
            record = self.records[self.active]
            if record['start_time_ry_inverse'] <= time <= record['start_time_ry_inverse']+record['dt_ry_inverse']+4*np.spacing(self.final.time):
                return self.active
        return next(i for i, record in enumerate(self.records)
                    if time <= record['start_time_ry_inverse']+record['dt_ry_inverse'])

    def counted(self, action):
        def apply(vectors):
            self.budget.action(vectors)
            return action(vectors)
        return apply

    def bind_response(self, *, q, coefficients, common_curvature, masses, charges,
                      polarizations, field_mass, ncells, density_shape, scalar_identity, coverage_identity):
        # Geometry, not source-star symmetry, selects the actual explicit target.
        lookup = {(b.slot.ispin, tuple(b.slot.mesh_index)): i for i, b in enumerate(self.bases)}
        ref = self.background.inputs.reference
        mesh = np.array([ref.monkhorst_pack['nk%d' % i] for i in (1, 2, 3)])
        step = np.rint(np.asarray(q)*mesh).astype(int)
        edges = []
        for i, basis in enumerate(self.bases):
            j = lookup[(basis.slot.ispin, tuple((np.asarray(basis.slot.mesh_index)+step) % mesh))]
            edge = RealBlochEdge.create(ref, basis, self.bases[j], q)
            edges.append(ResponseLink(i, j, LazyCoupling(self, edge)))
        sectors = [JointSector(slot.shape, self.initial.arrays(slot)[1], slot.weight, slot.component,
                               slot.basis_identity) for slot in self.initial.slots]
        fields = tuple(CartesianFieldActions(tuple(LazyFieldAxis(self, i, axis) for axis in range(3)))
                       for i in range(len(sectors))) if np.all(np.asarray(q) == 0) else ()
        self.response = JointBlochResponse(sectors, edges, ncells=ncells, q_fractional=q,
            q_bohr_inv=ref.reciprocal_bohr_inv @ q, density_shape=density_shape,
            kernel=self.background.apply_kernel, common_curvature=common_curvature,
            coefficients=coefficients, masses=masses, charges=charges, polarizations=polarizations,
            field_mass=field_mass, field_vertices=fields, scalar_identity=scalar_identity,
            coverage_identity=coverage_identity)
        return self.response

    @staticmethod
    def paired_moments_and_drive(coupling, left, right, left2, right2, direction, *, weight):
        return paired_matter_terms(coupling, left, right, left2, right2, direction, weight=weight)

    @staticmethod
    def paired_field_moments_and_drive(axes, left, right, left2, right2, field, *, weight):
        return paired_field_terms(axes.vertex, left, right, left2, right2, field, weight=weight)

    @contextmanager
    def __call__(self, index, time):
        i, prior, stage = self._stage(time)
        link = self.response.edges[index]
        source, target = link.source, link.target
        edge = link.coupling.edge
        h0s = physical_hamiltonian(self.background, self.bases[source], edge.ispin, self.block_size)
        h0t = h0s if source == target else physical_hamiltonian(self.background, self.bases[target], edge.ispin, self.block_size)
        cs = physical_coupling(self.background.inputs, self.bases[source], h0s.factors, self.block_size)
        ct = cs if source == target else physical_coupling(self.background.inputs, self.bases[target], h0t.factors, self.block_size)
        hs = PhysicalFrozenAction(h0s, cs, potential=stage.potential[:, edge.ispin-1], direction=stage.u, shift=stage.shift)
        ht = hs if source == target else PhysicalFrozenAction(h0t, ct, potential=stage.potential[:, edge.ispin-1], direction=stage.u, shift=stage.shift)
        sa, sad = self.counted(hs.apply), self.counted(hs.adjoint_apply)
        ta, tad = self.counted(ht.apply), self.counted(ht.adjoint_apply)
        ps = self.frame(source, time)
        pt = ps if source == target else self.frame(target, time)
        coupling = physical_finiteq(self.background.inputs, edge, h0s.factors, h0t.factors, self.block_size)
        axes = CartesianPhysicalFieldAxes(PhysicalFieldVertices(h0s)) if self.response.gamma else ()
        yield ps, pt, sa, sad, ta, tad, coupling, axes

    def initial_stage(self):
        # Initial H callbacks are not invoked by product preparation; nevertheless
        # they are actual physical actions, not identity/zero placeholders.
        def action(index, adjoint):
            def apply(vectors):
                outgoing = next(i for i, e in enumerate(self.response.edges) if e.source == index)
                with self(outgoing, 0.) as data:
                    return data[3 if adjoint else 2](vectors)
            return apply
        count = len(self.initial.slots)
        return JointStage(InitialFrames(self.initial), tuple(action(i, False) for i in range(count)),
                          tuple(action(i, True) for i in range(count)), 0., self.initial.identity)

    def time_grid(self, final_time, subdivisions):
        if type(subdivisions) is not int or subdivisions < 1 or not 0 <= final_time <= self.final.time:
            raise ValueError('positive subdivisions and covered nonnegative time required')
        times = [0.]
        for record in self.records:
            lo = record['start_time_ry_inverse']
            hi = min(lo+record['dt_ry_inverse'], final_time)
            if hi > lo:
                times.extend(np.linspace(lo, hi, subdivisions+1)[1:])
            if hi == final_time:
                break
        return np.asarray(times)
