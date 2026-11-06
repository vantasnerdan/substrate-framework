"""Bounded pinned-edge execution of the existing RK4/discrete-adjoint graph.

Only edge scheduling/factor lifetime changes: moments and non-feedback actions
share a physical edge construction, then a complete barrier precedes K/K† and
a bridge-only feedback pass. Immutable content receipts survive invocation loss.
"""
from contextlib import contextmanager
from pathlib import Path
import json
import os
from uuid import uuid4

import numpy as np

from .streamed_joint_response import StreamedJointResponse, MappedResponse, NAMES, resource_plan
from .streamed_joint_covariance import StreamedTwoSidedCovariance
from .partitioned_joint_slots import digest
from .partitioned_q_work import QWork, arrays
from .partitioned_q_storage import (QBudget, durable_publish, retire_payloads, sync_file,
                                   bounded_array_identity, release_mapped_pages)
from .streamed_material_state import authenticate, file_binding


def rk4_stage_times(time, dt, dual=False):
    return (time+dt, time+dt/2, time+dt/2, time) if dual else (time, time+dt/2, time+dt/2, time+dt)


def frame_times(times):
    if not len(times):
        raise ValueError('nonempty frame-use grid required')
    result = {float(times[-1])}
    for a, b in zip(times[:-1], times[1:]):
        result.update(float(t) for t in rk4_stage_times(a, b-a))
    return result


class PartitionedQResponse(StreamedJointResponse):
    def __init__(self, *args, cpus, work_directory, **kwargs):
        if not isinstance(kwargs.get('budget'), QBudget):
            raise ValueError('partitioned q requires live-storage QBudget accounting')
        super().__init__(*args, scratch_states=0, **kwargs)
        if not hasattr(self.edge_stage, 'stage_identity'):
            raise ValueError('actual immutable moving-stage identity is required')
        self.edge_work = QWork(Path(work_directory)/'edges', self.budget, cpus, self.identity)
        self.frames = QWork(Path(work_directory)/'frames', self.budget, cpus, self.identity)
        if hasattr(self.edge_stage, 'frame_work'):
            self.edge_stage.frame_work = self.frames
        self.storage = Path(work_directory).resolve()
        self.last_rhs = None
        self._owned_initials = set()
        self.frame_tail_times = None

    def prepare_time(self, time):
        if hasattr(self.edge_stage, 'prepare_time'):
            self.edge_stage.prepare_time(time)

    def retire_frames_except(self, needed_times, barrier):
        """Retire only with the complete remaining schedule, including BOTH legs."""
        needed_times = set(float(t) for t in needed_times)
        for directory in self.frames.directory.glob('*/generation-*'):
            marker = directory/'retire.json'
            if marker.exists():
                self.frames.retire(directory, json.loads(marker.read_text())['barrier'])
            elif (directory/'request.json').exists():
                request = json.loads((directory/'request.json').read_text())
                if float(request['input']['time']) not in needed_times:
                    self.frames.retire(directory, barrier)

    def input_identity(self, time, x, dual):
        r = self.r
        return dict(stage=self.edge_stage.stage_identity(time), dual=dual,
            dynamics=self.dynamics, q=r.q_fractional.tolist(), sqrt_n=float(r.sqrt_n),
            input=[bounded_array_identity(value) for value in x.values()],
            sectors=[dict(d=bounded_array_identity(s.d), weight=float(s.weight), component=s.density_component)
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
                paired = getattr(self.edge_stage, 'paired_moments_and_drive', None)
                if paired is None:
                    n1, g1 = c.moments(left, right, weight=weight)
                    n2, g2 = c.moments(left2, right2, weight=weight)
                    n, g = n1+n2, g1+g2
                    op = c.prepare(direction.reshape(-1, 3))
                    drive, adj = op.apply(ps), op.adjoint_apply(pt)
                else:
                    n, g, drive, adj = paired(c, left, right, left2, right2,
                                              direction.reshape(-1, 3), weight=weight)
                if r.gamma:
                    paired_field = getattr(self.edge_stage, 'paired_field_moments_and_drive', None)
                    if paired_field is None:
                        ell = weight*(self._field_moments(axes, left, right)
                                      +self._field_moments(axes, left2, right2))
                        for coefficient, axis in zip(field, axes):
                            drive += coefficient*axis.apply(ps)
                            adj += coefficient.conjugate()*axis.adjoint_apply(pt)
                    else:
                        ell, field_drive, field_adj = paired_field(
                            axes, left, right, left2, right2, field, weight=weight)
                        drive += field_drive
                        adj += field_adj
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

    def rhs(self, time, x, out=None, *, dual=False):
        self.prepare_time(time)
        identity = self.input_identity(time, x, dual)
        root = self.storage/'rhs'/digest(identity)
        root.mkdir(parents=True, exist_ok=True)
        active = [p for p in root.glob('generation-*') if not (p/'retired.json').exists()]
        if len(active) > 1:
            raise ValueError('concurrent q RHS generations are not supported')
        directory = active[0] if active else root/('generation-'+uuid4().hex)
        directory.mkdir(exist_ok=True)
        path = directory/'complete.json'
        if path.exists():
            record = json.loads(path.read_text())
            if record['identity'] != identity:
                raise ValueError('q RHS checkpoint identity differs')
        else:
            records = self.edge_work.execute('edge-free-moments', identity, len(self.r.edges),
                        lambda index: self.independent_edge(index, time, x, dual))
            free = [str(self.edge_work.last[1])]
            # Deterministic original edge order; never completion-order summation.
            n = np.zeros(self.r.density_shape, complex)
            g, ell = np.zeros(self.r.ndof, complex), np.zeros(3, complex)
            for item in records:
                value = arrays(item)
                n += value['n']
                g += value['g']
                ell += value['ell']
            if self.dynamics == 'screened-joint':
                kernel_identity = dict(input=identity, complete_density=bounded_array_identity(n),
                                       completed=[item['arrays'] for item in records])
                kernel_records = self.edge_work.execute('complete-global-kernel', kernel_identity, 1,
                    lambda index: dict(potential=(self.ka if dual else self.r.kernel)(n, q_bohr_inv=self.r.q_bohr_inv)))
                free.append(str(self.edge_work.last[1]))
                potential = arrays(kernel_records[0])['potential']
                feedback_identity = dict(input=identity, potential=bounded_array_identity(potential),
                                         completed=[item['arrays'] for item in records])
                records = self.edge_work.execute('edge-local-feedback', feedback_identity, len(self.r.edges),
                            lambda index: self.feedback_edge(index, time, potential, records, dual))
                del potential
            else:
                free = []
            # The final worker arrays themselves are immutable RHS progress. No
            # second Q-sized RHS copy or reusable mapped scratch is allocated.
            result = self._rhs_view(records)
            self.small_rhs(x, result, g, ell, dual=dual)
            record = dict(schema='partitioned-q-rhs-v2', identity=identity, records=records,
                small={name: [[float(z.real), float(z.imag)] for z in getattr(result, name)] for name in NAMES},
                free=free, final=str(self.edge_work.last[1]))
            durable_publish(path, record, self.budget)
            del result, value
        binding = file_binding(path)
        result = self._rhs_view(record['records'])
        for name in NAMES:
            setattr(result, name, np.asarray([complex(*z) for z in record['small'][name]]))
        # The durable RHS manifest is the barrier, not the mutable caller output.
        for retired in record['free']:
            self.edge_work.retire(retired, binding)
        self.budget.check_processes(os.getpid(), [])
        self.last_rhs = binding
        if out is not None:
            return out.combine((1, result))
        return result

    def _rhs_view(self, records):
        result = object.__new__(MappedResponse)
        result.response, result.path = self.r, None
        result.budget = self.budget
        result.ket, result.bra = [], []
        for record in records:
            for binding in record['arrays'].values():
                self.edge_work.verify(binding)
            value = arrays(record)
            result.ket.append(value['ket'])
            result.bra.append(value['bra'])
        for name, size in zip(NAMES, (self.r.ndof, self.r.ndof, self.r.nfield, self.r.nfield)):
            setattr(result, name, np.zeros(size, complex))
        return result

    def _retire_rhs(self, binding, barrier):
        record = json.loads(authenticate(binding).read_text())
        for directory in record['free']:
            self.edge_work.retire(directory, binding)
        self.edge_work.retire(record['final'], barrier)
        marker = Path(binding['path']).parent/'retired.json'
        if not marker.exists():
            durable_publish(marker, dict(barrier=barrier), self.budget)

    def _retire_state(self, binding, barrier):
        record = json.loads(authenticate(binding).read_text())
        owner = Path(binding['path']).parent
        if not (owner/'retire.json').exists():
            for index, leg in enumerate(record['legs']):
                if not Path(leg['path']).exists():
                    receipt = json.loads((owner/('retire-edge-%04d.json' % (index//2))).read_text())
                    authenticate(receipt['barrier'])
                    if leg not in receipt['payloads']:
                        raise ValueError('missing state leg has no durable edge replacement')
        else:
            barrier = json.loads((owner/'retire.json').read_text())['barrier']
        return retire_payloads(owner/'retire.json', barrier,
            [leg['path'] for leg in record['legs'] if Path(leg['path']).exists()], owner, self.budget)

    def _stage_cleanup(self, record, barrier):
        for binding in record['obsolete_states']:
            self._retire_state(binding, barrier)
        self._retire_rhs(record['rhs'], barrier)

    @staticmethod
    def _close_state(state):
        for value in state.ket+state.bra:
            value._mmap.close()

    def step(self, x, out, time, dt, *, dual=False):
        request = dict(identity=self.identity, input=[bounded_array_identity(v) for v in x.values()],
            time=float(time), dt=float(dt), dual=dual, dynamics=self.dynamics,
            provenance=getattr(self, 'provenance', None))
        result, _ = self._bounded_step(x, time, dt, self.storage/'single-steps'/digest(request), dual=dual)
        return out.combine((1, result))

    def _load_rhs(self, binding):
        record = json.loads(authenticate(binding).read_text())
        result = self._rhs_view(record['records'])
        for name in NAMES:
            setattr(result, name, np.asarray([complex(*z) for z in record['small'][name]]))
        for directory in record['free']:
            self.edge_work.retire(directory, binding)
        return result

    def _combine_stage(self, x, rhs, previous, stage, dt, dual, directory, rhs_binding):
        """Commit each edge's new pair before retiring that edge's old pair.

        An interrupted stage is a completed prefix plus an untouched remainder.
        Its durable RHS remains whole until the all-edge stage checkpoint exists.
        No mixed prefix is ever exposed as an adoptable full state.
        """
        old = [json.loads(authenticate(b).read_text()) for b in previous]
        accumulator = old[0] if stage else None
        out_path, input_path = directory/'accumulator', directory/'input'
        out_path.mkdir(parents=True, exist_ok=True)
        if stage < 3:
            input_path.mkdir(exist_ok=True)
        coefficient = 1 if dual else (dt/6 if stage in (0, 3) else dt/3)
        next_x = (dt/6 if stage == 2 else dt/3) if dual else 1
        next_rhs = (dt if stage == 0 else dt/2) if dual else (dt if stage == 2 else dt/2)
        for edge in range(len(self.r.edges)):
            path = directory/('edge-%04d.json' % edge)
            if path.exists():
                record = json.loads(path.read_text())
                if record['rhs'] != rhs_binding or record['edge'] != edge:
                    raise ValueError('stage edge completion identity differs')
                for binding in record['arrays']:
                    authenticate(binding)
            else:
                paths = [p/('leg-%04d.npy' % index)
                         for p in ([out_path, input_path] if stage < 3 else [out_path])
                         for index in (2*edge, 2*edge+1)]
                abandoned = [p for p in paths if p.exists()]
                if abandoned:
                    retire_payloads(directory/('abandoned-'+uuid4().hex+'.json'),
                        rhs_binding, abandoned, directory, self.budget)
                bindings = []
                for index, original, derivative in ((2*edge, x.ket[edge], rhs.ket[edge]),
                                                    (2*edge+1, x.bra[edge], rhs.bra[edge])):
                    prior = original if stage == 0 else np.load(
                        authenticate(accumulator['legs'][index]), mmap_mode='r', allow_pickle=False)
                    bindings.append(self.budget.array(out_path/('leg-%04d.npy' % index),
                        sum(a*v for a, v in ((1, prior), (coefficient, derivative)))))
                    if stage:
                        prior._mmap.close()
                    if stage < 3:
                        bindings.append(self.budget.array(input_path/('leg-%04d.npy' % index),
                            sum(a*v for a, v in ((next_x, original), (next_rhs, derivative)))))
                for binding in bindings:
                    sync_file(binding['path'])
                record = dict(edge=edge, rhs=rhs_binding, arrays=bindings)
                durable_publish(path, record, self.budget)
            barrier = file_binding(path)
            for binding, state in zip(previous, old):
                owner = Path(binding['path']).parent
                retire_payloads(owner/('retire-edge-%04d.json' % edge), barrier,
                    [leg['path'] for leg in state['legs'][2*edge:2*edge+2]], owner, self.budget)
        out = MappedResponse(self.r, out_path, self.budget)
        next_input = MappedResponse(self.r, input_path, self.budget) if stage < 3 else None
        for name in NAMES:
            original, derivative = getattr(x, name), getattr(rhs, name)
            prior = original if stage == 0 else np.asarray([complex(*z) for z in accumulator['small'][name]])
            setattr(out, name, sum(a*v for a, v in ((1, prior), (coefficient, derivative))))
            if next_input is not None:
                setattr(next_input, name, sum(a*v for a, v in ((next_x, original), (next_rhs, derivative))))
        return out, next_input

    def _publish_progress(self, state, time, kind, barrier):
        checkpoint = state.path/'checkpoint.json'
        provenance = getattr(self, 'provenance', None)
        if checkpoint.exists():
            binding = file_binding(checkpoint)
            restored, record = MappedResponse.restore(self.r, binding, self.budget, self.identity)
            same = (record['time'] == time and record['kind'] == kind and record['provenance'] == provenance
                    and all(np.array_equal(getattr(restored, name), getattr(state, name)) for name in NAMES))
            self._close_state(restored)
            if not same:
                raise ValueError('interrupted stage checkpoint metadata differs')
            return binding
        temporary = checkpoint.with_name(checkpoint.name+'.incomplete')
        if temporary.exists():
            retire_payloads(state.path/('abandoned-json-'+uuid4().hex+'.json'),
                barrier, [temporary], state.path, self.budget)
        return state.publish(time=time, identity=self.identity, budget=self.budget, kind=kind, provenance=provenance)

    def _bounded_step(self, x, time, dt, directory, *, dual=False, future_times=None):
        """Exact existing RK4 graph with per-edge durable accumulator replacement."""
        directory.mkdir(parents=True, exist_ok=True)
        request_path = directory/'request.json'
        request = dict(input=[bounded_array_identity(v) for v in x.values()], time=float(time), dt=float(dt),
                       dual=dual, identity=self.identity, provenance=getattr(self, 'provenance', None))
        if not request_path.exists():
            durable_publish(request_path, request, self.budget)
        elif json.loads(request_path.read_text()) != request:
            raise ValueError('RK4 stage request differs')
        if hasattr(self.edge_stage, 'set_interval'):
            self.edge_stage.set_interval(time, dt)
        completed = sorted(directory.glob('stage-*.json'))
        start, previous = 0, []
        if completed:
            record = json.loads(completed[-1].read_text())
            if ([path.name for path in completed] != ['stage-%d.json' % i for i in range(len(completed))]
                    or len(completed) > 4 or record['stage'] != len(completed)-1):
                raise ValueError('modified or noncontiguous RK4 stage progress')
            self._stage_cleanup(record, file_binding(completed[-1]))
            start = record['stage']+1
            if future_times is not None:
                self.retire_frames_except(set(future_times)|set(rk4_stage_times(time, dt, dual)[start:]),
                                          file_binding(completed[-1]))
            previous = [record['accumulator']]+([record['next_input']] if record['next_input'] else [])
            if start == 4:
                result, _ = MappedResponse.restore(self.r, record['accumulator'], self.budget, self.identity)
                return result, record['accumulator']
        if start == 0 and dual:
            seed_path = directory/'seed.json'
            if seed_path.exists():
                previous = [json.loads(seed_path.read_text())['checkpoint']]
            else:
                # Failed seed writes are regenerable from the actual input.
                garbage = list(directory.glob('seed-*/*.npy'))
                if garbage:
                    retire_payloads(directory/('seed-abandoned-'+uuid4().hex+'.json'),
                        file_binding(request_path), garbage, directory, self.budget)
                seed = MappedResponse(self.r, directory/('seed-'+uuid4().hex), self.budget, create=True)
                seed.combine((dt/6, x))
                binding = seed.publish(time=time+dt, identity=self.identity, budget=self.budget,
                    kind='rk4-stage-input', provenance=getattr(self, 'provenance', None))
                durable_publish(seed_path, dict(checkpoint=binding), self.budget)
                self._close_state(seed)
                previous = [binding]
        for stage in range(start, 4):
            stage_time = rk4_stage_times(time, dt, dual)[stage]
            rhs_path = directory/('rhs-stage-%d.json' % stage)
            if rhs_path.exists():
                rhs_binding = json.loads(rhs_path.read_text())['rhs']
                rhs = self._load_rhs(rhs_binding)
            else:
                stage_input = x if not previous else MappedResponse.restore(
                    self.r, previous[-1], self.budget, self.identity)[0]
                for value in x.values():
                    release_mapped_pages(value)
                rhs = self.rhs(stage_time, stage_input, dual=dual)
                rhs_binding = self.last_rhs
                durable_publish(rhs_path, dict(rhs=rhs_binding), self.budget)
                if stage_input is not x:
                    self._close_state(stage_input)
            out, next_input = self._combine_stage(x, rhs, previous, stage, dt, dual,
                directory/('combine-%d' % stage), rhs_binding)
            next_binding = None
            if next_input is not None:
                next_binding = self._publish_progress(next_input, stage_time, 'rk4-stage-input', rhs_binding)
                self._close_state(next_input)
            out_binding = self._publish_progress(out, time if dual else time+dt,
                ('dual' if dual else 'forward') if stage == 3 else 'rk4-accumulator', rhs_binding)
            record = dict(stage=stage, accumulator=out_binding, next_input=next_binding,
                          obsolete_states=previous, rhs=rhs_binding)
            barrier = durable_publish(directory/('stage-%d.json' % stage), record, self.budget)
            self._close_state(rhs)
            self._stage_cleanup(record, barrier)
            if future_times is not None:
                self.retire_frames_except(set(future_times)|set(rk4_stage_times(time, dt, dual)[stage+1:]), barrier)
            previous = [out_binding]+([next_binding] if next_binding else [])
            if stage < 3:
                self._close_state(out)
        return out, out_binding

    def own_initial(self, state):
        """Explicitly transfer a newly created invocation seed, never a restart."""
        if not (state.path.is_relative_to(self.directory.parent) or state.path.is_relative_to(self.storage)):
            raise ValueError('initial ownership is limited to this invocation or q storage')
        self._owned_initials.add(state.path)

    def _initial_request(self, key, time, kind):
        return dict(identity=self.identity, key=key, time=float(time), kind=kind,
                    dynamics=self.dynamics, provenance=getattr(self, 'provenance', None))

    def initial_state(self, key, initialize, *, time, kind):
        """Reuse an authenticated seed before allocating another full Q state."""
        request = self._initial_request(key, time, kind)
        root = self.storage/'initial'/digest(request)
        root.mkdir(parents=True, exist_ok=True)
        manifest = root/'request.json'
        if not manifest.exists():
            durable_publish(manifest, request, self.budget)
        elif json.loads(manifest.read_text()) != request:
            raise ValueError('q initial request differs')
        active = []
        for directory in root.glob('generation-*'):
            checkpoint = directory/'checkpoint.json'
            if (directory/'retire.json').exists():
                self._retire_state(file_binding(checkpoint), file_binding(manifest))
            elif checkpoint.exists():
                active.append(checkpoint)
            else:
                garbage = list(directory.glob('*.npy'))+list(directory.glob('*.incomplete'))
                if garbage:
                    retire_payloads(directory/('abandoned-'+uuid4().hex+'.json'),
                        file_binding(manifest), garbage, directory, self.budget)
        if len(active) > 1:
            raise ValueError('concurrent initial-state writers are not supported')
        if active:
            state, record = MappedResponse.restore(self.r, file_binding(active[0]), self.budget, self.identity)
            if record['kind'] != kind or record['time'] != time or record['provenance'] != request['provenance']:
                raise ValueError('stored q seed metadata differs')
        else:
            state = MappedResponse(self.r, root/('generation-'+uuid4().hex), self.budget, create=True)
            initialize(state)
            state.publish(time=time, identity=self.identity, budget=self.budget, kind=kind,
                          provenance=request['provenance'])
        state.seed_key = digest(request)
        self.own_initial(state)
        return state

    def _latest_interval(self, completed, times, dual):
        record = json.loads(completed[-1].read_text())
        if (len(completed) >= len(times) or record['ordinal'] != len(completed)-1
                or [path.name for path in completed] != ['interval-%06d.json' % i for i in range(len(completed))]):
            raise ValueError('modified or noncontiguous RK4 interval progress')
        binding, begin = record['checkpoint'], len(completed)
        current, checkpoint = MappedResponse.restore(self.r, binding, self.budget, self.identity)
        target = times[len(times)-1-begin if dual else begin]
        if checkpoint['time'] != target or checkpoint['kind'] != ('dual' if dual else 'forward'):
            raise ValueError('completed RK4 interval checkpoint time/direction differs')
        return current, binding, begin, record

    def resume_propagation(self, times, key, *, time, kind, dual=False):
        """Resume the actual source-bound request before reconstructing its seed."""
        seed_key = digest(self._initial_request(key, time, kind))
        matches = []
        for path in (self.storage/'propagation').glob('*/request.json'):
            request = json.loads(path.read_text())
            if path.parent.name != digest(request):
                raise ValueError('modified q propagation manifest')
            if (request.get('seed_key') == seed_key and request['times'] == list(times)
                    and request['identity'] == self.identity and request['dual'] == dual
                    and request['provenance'] == getattr(self, 'provenance', None)
                    and request['dynamics'] == self.dynamics):
                matches.append(path.parent)
        if not matches:
            return None
        if len(matches) != 1:
            raise ValueError('ambiguous source-bound q initial request')
        root = matches[0]
        completed = sorted(root.glob('interval-*.json'))
        if completed:
            current, binding, begin, record = self._latest_interval(completed, times, dual)
            if record['obsolete'] is not None:
                self._retire_state(record['obsolete'], file_binding(completed[-1]))
        else:
            candidates = [json.loads(path.read_text())['checkpoint'] for path in root.glob('seed-*.json')]
            candidates = [b for b in candidates if not (Path(b['path']).parent/'retire.json').exists()]
            if len(candidates) != 1:
                # Request creation may have been interrupted before seed adoption;
                # initial_state reuses its already-published seed without a copy.
                if not candidates:
                    return None
                raise ValueError('ambiguous active q seed')
            binding, begin = candidates[0], 0
            current, _ = MappedResponse.restore(self.r, binding, self.budget, self.identity)
        return self._continue_propagation(root, np.asarray(times, float), current, binding, begin, dual)

    def propagate(self, initial, times, directory, *, dual=False):
        times = np.asarray(times, float)
        if len(times) < 2 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError('strictly increasing explicit history-covered integration grid required')
        request = dict(identity=self.identity, initial=[bounded_array_identity(v) for v in initial.values()],
            seed_key=getattr(initial, 'seed_key', None), times=times.tolist(), dual=dual,
            dynamics=self.dynamics, provenance=getattr(self, 'provenance', None))
        root = self.storage/'propagation'/digest(request)
        root.mkdir(parents=True, exist_ok=True)
        if not (root/'request.json').exists():
            durable_publish(root/'request.json', request, self.budget)
        elif json.loads((root/'request.json').read_text()) != request:
            raise ValueError('propagation request differs')
        completed = sorted(root.glob('interval-*.json'))
        seed = initial.path/'checkpoint.json'
        owned_seed = None
        if seed.exists() and initial.path in self._owned_initials:
            if json.loads(seed.read_text())['kind'] in ('terminal-cotangent', 'coherent-initial'):
                owned_seed = file_binding(seed)
        if owned_seed is not None:
            seed_receipt = root/('seed-'+owned_seed['sha256']+'.json')
            if not seed_receipt.exists():
                durable_publish(seed_receipt, dict(checkpoint=owned_seed), self.budget)
            for receipt in root.glob('seed-*.json'):
                previous_seed = json.loads(receipt.read_text())['checkpoint']
                if previous_seed != owned_seed:
                    self._retire_state(previous_seed, owned_seed)
        current, current_binding, begin = initial, owned_seed, 0
        if completed:
            current, current_binding, begin, last = self._latest_interval(completed, times, dual)
            if last['obsolete'] is not None:
                self._retire_state(last['obsolete'], file_binding(completed[-1]))
            if owned_seed is not None:
                self._close_state(initial)
                self._retire_state(owned_seed, file_binding(completed[-1]))
        return self._continue_propagation(root, times, current, current_binding, begin, dual)

    def _continue_propagation(self, root, times, current, current_binding, begin, dual):
        indices = list(range(len(times)-2, -1, -1) if dual else range(len(times)-1))
        for ordinal in range(begin, len(indices)):
            index = indices[ordinal]
            future = None
            if self.frame_tail_times is not None:
                future = set(self.frame_tail_times)
                for later in indices[ordinal+1:]:
                    future.update(rk4_stage_times(times[later], times[later+1]-times[later], dual))
            final, binding = self._bounded_step(current, times[index], times[index+1]-times[index],
                root/('step-%06d' % ordinal), dual=dual, future_times=future)
            record = dict(ordinal=ordinal, checkpoint=binding, obsolete=current_binding)
            barrier = durable_publish(root/('interval-%06d.json' % ordinal), record, self.budget)
            if current_binding is not None:
                self._close_state(current)
            del current
            if current_binding is not None:
                self._retire_state(current_binding, barrier)
            current, current_binding = final, binding
        if self.frame_tail_times is not None:
            self.retire_frames_except(self.frame_tail_times, current_binding)
        return current, [current_binding]

    def terminal_probe(self, time, probe, out):
        probe.validate(self.r)
        self.own_initial(out)
        self.budget.record_write(sum(value.nbytes for value in out.ket+out.bra))
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


class PartitionedTwoSidedCovariance(StreamedTwoSidedCovariance):
    """The unchanged two-time contraction with bounded seed/continuation ownership."""
    def pair(self, left, right, times, second_times, directory, *, left_checkpoint=None, right_checkpoint=None):
        previous = self.executor.frame_tail_times
        self._pair_tail = [frame_times(second_times), set()]
        try:
            result = super().pair(left, right, times, second_times, directory,
                                 left_checkpoint=left_checkpoint, right_checkpoint=right_checkpoint)
            self.executor.retire_frames_except((), result['right_checkpoints'][-1])
            return result
        finally:
            self.executor.frame_tail_times = previous
            del self._pair_tail

    def pullback(self, probe, times, directory, *, checkpoint=None):
        if hasattr(self, '_pair_tail'):
            self.executor.frame_tail_times = self._pair_tail.pop(0)
        if checkpoint is not None:
            return super().pullback(probe, times, directory, checkpoint=checkpoint)
        times = np.asarray(times, float)
        if not len(times) or times[0] != 0 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError('pullback grid must start at original t0')
        executor = self.executor
        probe.validate(executor.r)
        previous = getattr(executor, 'provenance', None)
        executor.provenance = dict(previous or {}, pullback=dict(probe=probe.identity, units=probe.units,
            terminal_time=float(times[-1]), grid=times.tolist()))
        key = dict(probe=probe.identity, units=probe.units,
                   small={name: bounded_array_identity(getattr(probe, name)) for name in NAMES})
        try:
            if len(times) > 1:
                resumed = executor.resume_propagation(times, key, time=times[-1], kind='terminal-cotangent', dual=True)
                if resumed is not None:
                    return resumed
            if hasattr(executor.edge_stage, 'clear_interval'):
                executor.edge_stage.clear_interval()
            terminal = executor.initial_state(key, lambda state:executor.terminal_probe(times[-1], probe, state),
                time=times[-1], kind='terminal-cotangent')
            if len(times) == 1:
                return terminal, [file_binding(terminal.path/'checkpoint.json')]
            return executor.propagate(terminal, times, Path(directory)/'pullback', dual=True)
        finally:
            executor.provenance = previous


def partition_plan(response, times, second_times, workers, *, mode='pair', aliased_frame_times=(), dynamics='screened-joint'):
    """Storage ownership, not a claim that a proposed time grid has real history."""
    if mode not in ('forward', 'pair', 'all') or not 1 <= workers <= 5:
        raise ValueError('explicit execution mode and allocated worker count required')
    if dynamics not in ('screened-joint', 'bare-electronic'):
        raise ValueError('known q dynamics required for accounting')
    grids = [np.asarray(ts, float) for ts in (times, second_times)]
    if any(ts.ndim != 1 or not len(ts) or not np.isfinite(ts).all()
           or np.any(np.diff(ts) <= 0) for ts in grids):
        raise ValueError('finite increasing proposed grids required')
    plan = resource_plan(response)
    q, edge = plan['one_q_leg_bytes'], plan['largest_edge_leg_bytes']
    # Same event ordering as production: forward/readout, then both actual
    # pullbacks; terminal frames survive until the first durable reverse stage.
    from collections import Counter
    legs = ([(grids[0], False)] if mode in ('forward', 'all') else [])+(
            [(grids[0], True), (grids[1], True)] if mode in ('pair', 'all') else [])
    events = []
    for leg, (grid, dual) in enumerate(legs):
        if dual:
            events.append((float(grid[-1]), leg, False))
        indices = range(len(grid)-2, -1, -1) if dual else range(len(grid)-1)
        for index in indices:
            events.extend((float(t), leg, True) for t in rk4_stage_times(grid[index], grid[index+1]-grid[index], dual))
        if not dual:
            events.append((float(grid[-1]), leg, True))
    stage_times = {t for t, _, _ in events}
    aliases = set(float(t) for t in aliased_frame_times)
    cached_times = stage_times-aliases
    frames = sum(16*int(np.prod(s.shape)) for s in response.sectors)
    count = len(response.edges)
    moments = count*16*(int(np.prod(response.density_shape))+response.ndof+3)
    coupled = dynamics == 'screened-joint'
    potential = 16*int(np.prod(response.density_shape)) if coupled else 0
    active = max(5*q+moments+potential, 4*q+2*edge) if coupled else 4*q+moments+2*edge
    retained = {'forward': 0, 'pair': 1, 'all': 2}[mode]
    remaining = Counter(t for t, _, _ in events)
    live_frames, peak_frames, peak_offset, peak_event = set(), 0, 0, None
    for time, leg, retire_after in events:
        remaining[time] -= 1
        if time not in aliases:
            live_frames.add(time)
        peak_frames = max(peak_frames, len(live_frames))
        offset = leg*q+len(live_frames)*frames
        if offset > peak_offset:
            peak_offset = offset
            peak_event = dict(time=time, prior_output_states=leg, live_frame_times=sorted(live_frames))
        if retire_after:
            live_frames = {t for t in live_frames if remaining[t] > 0}
    frame_payload = peak_frames*frames
    frame_writes = len(cached_times)*frames
    forward_intervals = len(grids[0])-1 if mode in ('forward', 'all') else 0
    dual_intervals = sum(len(ts)-1 for ts in grids) if mode in ('pair', 'all') else 0
    intervals = forward_intervals+dual_intervals
    seed_states = {'forward': 1, 'pair': 2, 'all': 3}[mode]
    payload_allocations = ((15 if coupled else 11)*intervals+dual_intervals+seed_states)*q+4*intervals*(moments+potential)+frame_writes
    payload_writes = payload_allocations+(dual_intervals+(2 if mode in ('pair', 'all') else 0))*q
    arrays_written = ((42*count+4) if coupled else 34*count)*intervals+2*count*(dual_intervals+seed_states)+len(cached_times)*len(response.sectors)
    plan.update(mode=mode, dynamics=dynamics, workers=workers, scratch_states=0, scratch_leg_bytes=0,
        policy='per-edge immutable RK4 accumulator/input replacement; no reusable full-state scratch',
        one_q_state_bytes=q, largest_edge_two_leg_bytes=edge,
        actual_frame_bytes_per_time=frames, cached_stage_times=sorted(cached_times),
        genuine_boundary_alias_times=sorted(stage_times-set(cached_times)),
        retained_frame_payload_upper_bound=frame_payload,
        cumulative_frame_payload_writes_upper_bound=frame_writes, peak_frame_and_output_event=peak_event,
        independent_moments_bytes=moments, complete_potential_bytes=potential,
        active_stage_peak_payload_bytes=active, retained_prior_output_states=retained,
        peak_live_payload_upper_bound=active+peak_offset,
        resume_seed_overlap_payload_bytes=0,
        resume_peak_live_payload_upper_bound=active+peak_offset,
        final_live_payload_upper_bound=seed_states*q,
        cumulative_payload_writes_upper_bound=payload_writes,
        cumulative_payload_allocations_upper_bound=payload_allocations,
        successful_array_header_reservation_upper_bound=512*arrays_written,
        peak_expression='existing unrelated bytes + max_over_exact_schedule(prior outputs*Q + live frames) + '+(
            'max(5Q+M+P, 4Q+2E)' if coupled else '(4Q+M+2E)')+
            ' + live headers/receipts/reservation slack; production resume adopts existing seed/continuation without a Q copy',
        frame_policy='After a durable stage, retire only times absent from ALL remaining stages/readout and BOTH pair legs. Unscheduled standalone calls retain frames. Genuine material boundary frames alias existing history; no assumed future endpoint aliases.',
        checkpoint_policy='Only current/final propagation payload retained; obsolete stage/interval payload retired after durable replacement. Compact checkpoints/attempt/cost/retirement receipts remain. Returned checkpoint list contains the current final checkpoint only.',
        rhs_per_forward_or_dual_interval=4,
        storage_policy='QBudget caps live reservations separately from cumulative write/action cost. Failed uncommitted arrays reclaimed before retry; completed RHS/edge progress retained until durable successor.',
        recovery='Authenticate latest interval, then whole stage or per-edge completed prefix plus immutable RHS. No mutable scratch adoption. One writer per work directory.',
        limitations='Payload expressions exclude JSON/filesystem overhead and failed-attempt reservation slack, which are measured/charged. Shared times needed by the second leg may still dominate disk on long grids. Unrelated or explicitly external checkpoints are never retired.',
        memory_policy='Parent plus pinned child RSS watchdog; no all-q allocation, no PW-square matrices.')
    return plan


def initial_partition_plan(binding, *, q, mesh, density_shape, times, second_times, workers, mode='pair', dynamics='screened-joint'):
    """Plan from actual t0 boundary metadata without constructing moving history.

    ``mesh`` and ``density_shape`` must be the actual case's reference metadata;
    proposed times are labeled proposals, never positive-time evidence.
    No frame arrays, Hamiltonians or mean evolution are loaded/executed.
    """
    from types import SimpleNamespace
    initial = json.loads(authenticate(binding).read_text())
    if (initial.get('schema') != 'streamed-joint-mean-v1' or initial.get('complete') is not True
            or initial.get('time_ry_inverse') != 0 or initial.get('history')):
        raise ValueError('complete actual initial joint metadata required, not invented moving history')
    mesh, q = np.asarray(mesh), np.asarray(q, float)
    if (mesh.shape != (3,) or not np.issubdtype(mesh.dtype, np.integer) or np.any(mesh <= 0)
            or q.shape != (3,) or not np.isfinite(q).all()
            or not np.allclose(q*mesh, np.rint(q*mesh), atol=1e-12, rtol=0)):
        raise ValueError('actual integer mesh and commensurate q required')
    slots = initial['slots']
    labels = []
    for slot in slots:
        spin, indices = slot['identity'].split(':')
        label = (int(spin), tuple(int(i) for i in indices.split(',')))
        if len(label[1]) != 3 or any(i < 0 or i >= n for i, n in zip(label[1], mesh)):
            raise ValueError('actual explicit slot mesh identity differs')
        if len(slot['shape']) != 2 or any(type(n) is not int or n <= 0 for n in slot['shape']):
            raise ValueError('positive original rectangular slot shape required')
        labels.append(label)
    lookup = {label: i for i, label in enumerate(labels)}
    if len(lookup) != len(slots) or len(slots) != len({s for s, _ in labels})*int(np.prod(mesh)):
        raise ValueError('complete explicit spin/mesh slot coverage required')
    step = np.rint(q*mesh).astype(int)
    edges = [SimpleNamespace(source=i, target=lookup[(spin, tuple((np.asarray(k)+step) % mesh))])
             for i, (spin, k) in enumerate(labels)]
    if len(density_shape) != 2 or any(type(n) is not int or n <= 0 for n in density_shape):
        raise ValueError('actual positive density grid/component dimensions required')
    response = SimpleNamespace(sectors=[SimpleNamespace(shape=s['shape']) for s in slots],
        edges=edges, ndof=len(initial['u_bohr']), density_shape=tuple(density_shape))
    plan = partition_plan(response, times, second_times, workers, mode=mode, dynamics=dynamics, aliased_frame_times=(0.,))
    if times[0] != 0 or second_times[0] != 0:
        raise ValueError('proposed initial-case grids must start at actual t0')
    plan.update(initial_metadata=binding, model_identity=initial['model_identity'], q=q.tolist(),
        mesh=mesh.tolist(), density_shape=list(density_shape), numerical_work_executed=False,
        history_available=False, proposed_time_grid=list(times), proposed_second_time_grid=list(second_times),
        execution_requires='Genuine completed positive-time history and current source-bound finite oracle; this static plan does not authorize propagation.')
    return plan
