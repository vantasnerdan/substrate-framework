"""Exact discrete q pullbacks retaining global contexts, never full Q states.

The unchanged edge equations are replayed from the original terminal probe.
This exchanges additional charged physical actions for bounded disk payloads;
it is not a cheaper mean trajectory or a mapped-state restart format.
"""
from pathlib import Path
import json

import numpy as np

from .partitioned_q_response import rk4_stage_times
from .partitioned_q_storage import (bounded_array_identity, durable_publish,
                                   release_mapped_pages)
from .partitioned_q_work import QWork
from .source_local_q_context import complete_context
from .source_local_q_state import (small_state, encode_small, decode_small,
    combine_small, next_coefficients, edge_state, terminal_edge, replay_interval)
from .streamed_joint_response import NAMES
from .streamed_material_state import authenticate, file_binding
from .source_local_q_edge import dual_moments


def pullback_grid(times):
    times = np.asarray(times, float)
    if (times.ndim != 1 or not len(times) or times[0] != 0
            or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)):
        raise ValueError('finite increasing pullback grid must start at original t0')
    return times


def prepare_pullback_frames(executor, grids):
    times = {float(grid[-1]) for grid in grids}
    for grid in grids:
        for a, b in zip(grid[:-1], grid[1:]):
            times.update(rk4_stage_times(a, b-a, True))
    if hasattr(executor.edge_stage, 'clear_interval'):
        executor.edge_stage.clear_interval()
    for time in sorted(times):
        executor.prepare_time(time)
    return sorted(times)


class _Leg:
    def __init__(self, owner, index):
        self.owner, self.index = owner, index

    def __getitem__(self, edge):
        return self.owner.edge(edge)[self.index]


class SourceLocalPullback:
    """Lazy original-time cotangent; at most one reconstructed edge is cached."""
    def __init__(self, executor, probe, times, intervals, small, binding):
        self.executor, self.probe, self.times = executor, probe, times
        self.intervals, self.binding = intervals, binding
        self.ket, self.bra = _Leg(self, 0), _Leg(self, 1)
        for name in NAMES:
            setattr(self, name, getattr(small, name))
        self._last = None
        self._pair = None

    def edge(self, index):
        if self._last != index:
            self._last, self._pair = None, None
            pair = terminal_edge(self.executor, self.probe, self.times[-1], index)
            for interval in self.intervals:
                pair = replay_interval(self.executor, index, pair, interval)
            self._last, self._pair = index, pair
            for interval in self.intervals:
                for context in interval['contexts']:
                    release_mapped_pages(context['potential'])
        return self._pair

    def release(self):
        self._last, self._pair = None, None


class SourceLocalQResponse:
    def __init__(self, executor, directory):
        self.executor = executor
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.work = QWork(self.directory/'work', executor.budget,
                          executor.edge_work.cpus, executor.identity)

    def _stage_identities(self, times):
        result = []
        for a, b in zip(times[:-1], times[1:]):
            if hasattr(self.executor.edge_stage, 'set_interval'):
                self.executor.edge_stage.set_interval(a, b-a)
            result.append([self.executor.edge_stage.stage_identity(time)
                           for time in rk4_stage_times(a, b-a, True)])
        return result

    def _input(self, probe, times, completed, current, stage, index):
        pair = terminal_edge(self.executor, probe, times[-1], index)
        for interval in completed:
            pair = replay_interval(self.executor, index, pair, interval)
        return replay_interval(self.executor, index, pair, current, prefix=stage)

    def _moments(self, probe, times, completed, current, stage, small, index):
        pair = self._input(probe, times, completed, current, stage, index)
        time = rk4_stage_times(current['time'], current['dt'], True)[stage]
        return dual_moments(
            self.executor, index, time, edge_state(index, pair, small))


    def pullback(self, probe, times):
        executor, r = self.executor, self.executor.r
        times = pullback_grid(times)
        probe.validate(r)
        stages = self._stage_identities(times)
        request = dict(schema='source-local-q-request-v1', identity=executor.identity,
            provenance=getattr(executor, 'provenance', None), dynamics=executor.dynamics,
            q=r.q_fractional.tolist(), times=times.tolist(), stage_identities=stages,
            probe=dict(identity=probe.identity, units=probe.units, small=encode_small(probe)),
            sectors=[dict(d=bounded_array_identity(s.d), weight=float(s.weight),
                          component=s.density_component) for s in r.sectors])
        path = self.directory/'request.json'
        if path.exists():
            if json.loads(path.read_text()) != request:
                raise ValueError('source-local pullback source/probe/grid differs')
        else:
            durable_publish(path, request, executor.budget)
        request_binding = file_binding(path)
        prepare_pullback_frames(executor, (times,))
        completed = []
        terminal = small_state(r, {name: np.array(getattr(probe, name), copy=True) for name in NAMES})
        for index in range(len(times)-2, -1, -1):
            time, dt = float(times[index]), float(times[index+1]-times[index])
            current = dict(time=time, dt=dt, contexts=[])
            small, out = combine_small(r, (dt/6, terminal)), terminal
            for stage, stage_time in enumerate(rk4_stage_times(time, dt, True)):
                identity = dict(pullback=request_binding, interval=index, stage=stage,
                    stage_time=stage_time, input_small=encode_small(small),
                    previous=[ctx['binding'] for interval in completed+[current]
                              for ctx in interval['contexts']])
                context = complete_context(executor, self.work, self.directory/('interval-%04d' % index)/
                    ('stage-%d.json' % stage), identity,
                    lambda edge: self._moments(probe, times, completed, current, stage, small, edge), small)
                current['contexts'].append(context)
                derivative = small_state(r, decode_small(context['record']['derivative_small']))
                out = combine_small(r, (1, out), (1, derivative))
                if stage < 3:
                    a, b = next_coefficients(stage, dt)
                    small = combine_small(r, (a, terminal), (b, derivative))
            terminal = out
            completed.append(current)
        record = dict(schema='source-local-q-pullback-v1', complete=True, request=request_binding,
            time=0., terminal_time=float(times[-1]), small=encode_small(terminal),
            intervals=[dict(time=interval['time'], dt=interval['dt'],
                contexts=[ctx['binding'] for ctx in interval['contexts']]) for interval in completed],
            storage='global kernel contexts and canonical arrays; electronic legs reconstructed, not mapped')
        path = self.directory/'complete.json'
        if path.exists():
            if json.loads(path.read_text()) != record:
                raise ValueError('source-local complete pullback differs')
        else:
            durable_publish(path, record, executor.budget)
        return SourceLocalPullback(executor, probe, times, completed, terminal, file_binding(path))

    def restore(self, probe, times, binding):
        record = json.loads(authenticate(binding).read_text())
        if record.get('schema') != 'source-local-q-pullback-v1' or record.get('complete') is not True:
            raise ValueError('a complete source-local pullback is required, not a mapped state')
        if Path(binding['path']).resolve() != self.directory/'complete.json':
            raise ValueError('source-local pullback belongs to a different directory')
        authenticate(record['request'])
        for interval in record['intervals']:
            for context in interval['contexts']:
                value = json.loads(authenticate(context).read_text())
                if value['kernel'] is not None:
                    for array in value['kernel']['arrays'].values():
                        authenticate(array)
        return self.pullback(probe, times)
