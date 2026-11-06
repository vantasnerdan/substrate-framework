"""One-edge cotangents and the existing reverse-RK4 arithmetic order.

No full response state is allocated.  Small canonical variables remain global;
rectangular electronic legs belong to one explicit directed edge only.
"""
from types import SimpleNamespace

import numpy as np

from .streamed_joint_response import NAMES
from .source_local_q_edge import dual_free


def small_state(response, values=None):
    sizes = (response.ndof, response.ndof, response.nfield, response.nfield)
    return SimpleNamespace(**{
        name: np.zeros(size, complex) if values is None else values[name]
        for name, size in zip(NAMES, sizes)})


def encode_small(state):
    return {name: [[float(z.real), float(z.imag)] for z in getattr(state, name)]
            for name in NAMES}


def decode_small(record):
    return {name: np.asarray([complex(*z) for z in record[name]]) for name in NAMES}


def combine_small(response, *terms):
    return small_state(response, {
        name: sum(a*getattr(value, name) for a, value in terms) for name in NAMES})


def next_coefficients(stage, dt):
    """Coefficients of terminal seed and current derivative, in reverse order."""
    if stage not in (0, 1, 2):
        raise ValueError('only the first three reverse stages have a successor')
    return (dt/6 if stage == 2 else dt/3), (dt if stage == 0 else dt/2)


def combine_pair(*terms):
    return tuple(sum(a*value[index] for a, value in terms) for index in (0, 1))


def edge_state(index, pair, small):
    return SimpleNamespace(ket={index: pair[0]}, bra={index: pair[1]},
                           **{name: getattr(small, name) for name in NAMES})


def terminal_edge(executor, probe, time, index):
    r, stage = executor.r, executor.edge_stage
    edge = r.edges[index]
    if hasattr(stage, 'clear_interval'):
        stage.clear_interval()
    if hasattr(stage, 'frame'):
        ps, pt = stage.frame(edge.source, time), stage.frame(edge.target, time)
        weight = r.sectors[edge.source].weight
        return (weight*probe.actions[index](ps),
                weight*probe.adjoints[index](pt).conj().T)
    with stage(index, time) as (ps, pt, *unused):
        weight = r.sectors[edge.source].weight
        return (weight*probe.actions[index](ps),
                weight*probe.adjoints[index](pt).conj().T)


def dual_feedback(executor, index, time, pair, potential):
    if executor.dynamics != 'screened-joint':
        return pair
    source = executor.r.sectors[executor.r.edges[index].source]
    with executor.feedback_bridge(index, time) as (ps, pt, bridge):
        v = potential[:, source.density_component]
        drive = bridge.local_apply(v, ps)
        adjoint = bridge.local_adjoint_apply(v, pt)
        return (pair[0]+source.weight*drive,
                pair[1]+source.weight*adjoint.conj().T)


def replay_interval(executor, index, terminal, interval, *, prefix=None):
    """Replay committed contexts; prefix returns that stage's input, not output."""
    time, dt = interval['time'], interval['dt']
    if hasattr(executor.edge_stage, 'set_interval'):
        executor.edge_stage.set_interval(time, dt)
    current = combine_pair((dt/6, terminal))
    out = terminal
    for stage, context in enumerate(interval['contexts']):
        if prefix == stage:
            return current
        small = small_state(executor.r, decode_small(context['record']['input_small']))
        free = dual_free(executor, index, context['record']['stage_time'],
                         edge_state(index, current, small))
        derivative = dual_feedback(executor, index, context['record']['stage_time'],
                                   free, context['potential'])
        out = combine_pair((1, out), (1, derivative))
        if stage < 3:
            a, b = next_coefficients(stage, dt)
            current = combine_pair((a, terminal), (b, derivative))
    if prefix == len(interval['contexts']):
        return current
    if len(interval['contexts']) != 4:
        raise ValueError('a final cotangent requires all four committed contexts')
    return out
