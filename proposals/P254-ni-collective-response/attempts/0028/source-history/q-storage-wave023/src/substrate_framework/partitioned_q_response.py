"""Bounded pinned-edge execution of the existing RK4/discrete-adjoint graph.

Only edge scheduling/factor lifetime changes: moments and non-feedback actions
share a physical edge construction, then a complete barrier precedes K/K† and
a bridge-only feedback pass. Immutable content receipts survive invocation loss.
"""
from contextlib import contextmanager
from pathlib import Path
import os

import numpy as np

from .streamed_joint_response import StreamedJointResponse, resource_plan
from .partitioned_joint_slots import array_identity
from .partitioned_q_work import QWork, arrays


class PartitionedQResponse(StreamedJointResponse):
    def __init__(self, *args, cpus, work_directory, **kwargs):
        super().__init__(*args, **kwargs)
        if not hasattr(self.edge_stage, 'stage_identity'):
            raise ValueError('actual immutable moving-stage identity is required')
        self.edge_work = QWork(Path(work_directory)/'edges', self.budget, cpus, self.identity)
        self.frames = QWork(Path(work_directory)/'frames', self.budget, cpus, self.identity)
        if hasattr(self.edge_stage, 'frame_work'):
            self.edge_stage.frame_work = self.frames

    def prepare_time(self, time):
        if hasattr(self.edge_stage, 'prepare_time'):
            self.edge_stage.prepare_time(time)

    def input_identity(self, time, x, dual):
        r = self.r
        return dict(stage=self.edge_stage.stage_identity(time), dual=dual,
            dynamics=self.dynamics, q=r.q_fractional.tolist(), sqrt_n=float(r.sqrt_n),
            input=[array_identity(value) for value in x.values()],
            sectors=[dict(d=array_identity(s.d), weight=float(s.weight), component=s.density_component)
                     for s in r.sectors])

    def independent_edge(self, index, time, x, dual):
        r, e = self.r, self.r.edges[index]
        s, t = r.sectors[e.source], r.sectors[e.target]
        coupled = self.dynamics == 'screened-joint'
        n, g, ell = np.zeros(r.density_shape, complex), np.zeros(r.ndof, complex), np.zeros(3, complex)
        with self.edge_stage(index, time) as (ps, pt, hs, hsa, ht, hta, c, axes):
            if dual:
                left, right = 1j*x.ket[index] @ s.d, ps
                left2, right2 = -1j*pt @ t.d, x.bra[index].conj().T
                weight = 1.
                direction = -x.p
                field = -r.sqrt_n*r.polarizations.T @ x.momentum if r.gamma else np.zeros(3)
            else:
                left, right, left2, right2 = x.ket[index], ps, pt, x.bra[index].conj().T
                weight = s.weight
                direction = x.u
                field = r.sqrt_n*r.polarizations.T @ x.coordinate if r.gamma else np.zeros(3)
            if coupled:
                n1, g1 = c.moments(left, right, weight=weight)
                n2, g2 = c.moments(left2, right2, weight=weight)
                n, g = n1+n2, g1+g2
                if r.gamma:
                    ell = weight*(self._field_moments(axes, left, right)
                                  +self._field_moments(axes, left2, right2))
                op = c.prepare(direction.reshape(-1, 3))
                drive, adj = op.apply(ps), op.adjoint_apply(pt)
                for coefficient, axis in zip(field, axes):
                    drive += coefficient*axis.apply(ps)
                    adj += coefficient.conjugate()*axis.adjoint_apply(pt)
            else:
                drive = np.zeros_like(x.ket[index])
                adj = np.zeros_like(x.bra[index].conj().T)
            if dual:
                ket = 1j*hta(x.ket[index])+s.weight*drive
                bra = -1j*hs(x.bra[index].conj().T).conj().T+s.weight*adj.conj().T
            else:
                ket = -1j*(ht(x.ket[index])+drive @ s.d)
                bra = 1j*(hsa(x.bra[index].conj().T).conj().T+t.d @ adj.conj().T)
        return dict(ket=ket, bra=bra, n=n, g=g, ell=ell)

    @contextmanager
    def feedback_bridge(self, index, time):
        if hasattr(self.edge_stage, 'feedback'):
            with self.edge_stage.feedback(index, time) as data:
                yield data
        else:
            # Finite consumers expose their exact original bridge through stage.
            with self.edge_stage(index, time) as data:
                yield data[0], data[1], data[6].bridge

    def feedback_edge(self, index, time, potential, records, dual):
        r, e = self.r, self.r.edges[index]
        s, t = r.sectors[e.source], r.sectors[e.target]
        old = arrays(records[index])
        with self.feedback_bridge(index, time) as (ps, pt, bridge):
            v = potential[:, s.density_component]
            drive, adj = bridge.local_apply(v, ps), bridge.local_adjoint_apply(v, pt)
            if dual:
                ket, bra = old['ket']+s.weight*drive, old['bra']+s.weight*adj.conj().T
            else:
                ket, bra = old['ket']-1j*drive @ s.d, old['bra']+1j*t.d @ adj.conj().T
        return dict(ket=ket, bra=bra)

    def rhs(self, time, x, out, *, dual=False):
        self.prepare_time(time)
        identity = self.input_identity(time, x, dual)
        records = self.edge_work.execute('edge-free-moments', identity, len(self.r.edges),
                    lambda index: self.independent_edge(index, time, x, dual))
        # Deterministic original edge order; never completion-order summation.
        n = np.zeros(self.r.density_shape, complex)
        g, ell = np.zeros(self.r.ndof, complex), np.zeros(3, complex)
        for record in records:
            value = arrays(record)
            n += value['n']
            g += value['g']
            ell += value['ell']
        if self.dynamics == 'screened-joint':
            kernel_identity = dict(input=identity, complete_density=array_identity(n),
                                   completed=[record['arrays'] for record in records])
            kernel_records = self.edge_work.execute('complete-global-kernel', kernel_identity, 1,
                lambda index: dict(potential=(self.ka if dual else self.r.kernel)(n, q_bohr_inv=self.r.q_bohr_inv)))
            potential = arrays(kernel_records[0])['potential']
            feedback_identity = dict(input=identity, potential=array_identity(potential),
                                     completed=[record['arrays'] for record in records])
            records = self.edge_work.execute('edge-local-feedback', feedback_identity, len(self.r.edges),
                    lambda index: self.feedback_edge(index, time, potential, records, dual))
        self.budget.check_processes(os.getpid(), [])
        for index, record in enumerate(records):
            value = arrays(record)
            out.ket[index][:], out.bra[index][:] = value['ket'], value['bra']
        return self.small_rhs(x, out, g, ell, dual=dual)

    def terminal_probe(self, time, probe, out):
        probe.validate(self.r)
        self.prepare_time(time)
        if not hasattr(self.edge_stage, 'frame'):
            return super().terminal_probe(time, probe, out)
        # Probe construction needs baseline frames, not an unused H0/g stage.
        for index, edge in enumerate(self.r.edges):
            self.budget.check()
            ps = self.edge_stage.frame(edge.source, time)
            pt = self.edge_stage.frame(edge.target, time)
            weight = self.r.sectors[edge.source].weight
            out.ket[index][:] = weight*probe.actions[index](ps)
            out.bra[index][:] = weight*probe.adjoints[index](pt).conj().T
        for name in ('u', 'p', 'coordinate', 'momentum'):
            getattr(out, name)[:] = getattr(probe, name)
        return out


def partition_plan(response, times, second_times, workers):
    plan = resource_plan(response)
    stage_times = set(float(t) for ts in (times, second_times) for t in ts)
    stage_times.update(float((a+b)/2) for ts in (times, second_times) for a,b in zip(ts[:-1], ts[1:]))
    frames = sum(16*int(np.prod(s.shape)) for s in response.sectors)
    moments = len(response.edges)*16*(int(np.prod(response.density_shape))+response.ndof+3)
    plan.update(workers=workers, actual_frame_bytes_per_time=frames,
        retained_frame_payload_upper_bound=len(stage_times)*frames,
        retained_payload_per_distinct_rhs=2*plan['one_q_leg_bytes']+moments+16*int(np.prod(response.density_shape)),
        rhs_per_forward_or_dual_interval=4,
        storage_policy='Retain authenticated frames and completed edge work; hard aggregate disk cap includes failures and existing files. No eviction. Headers, journals and immutable RK4 output checkpoints add storage.',
        recovery='Rerun identical request; regenerate cheap RK4 combinations, reuse authenticated completed edge work. Never restore mutable scratch.',
        memory_policy='Parent plus pinned child RSS watchdog; no all-q allocation, no PW-square matrices.')
    return plan
