"""Readonly committed-state selected observable, one genuine QWork scalar per edge.

The checkpoint's historical producer and the currently executing consumer are
separate authorities. No tangent propagation, PW-square density, D insertion,
real projection, cell multiplier, or reordered floating-point reduction occurs.
"""
from hashlib import sha256
import inspect
import json
import marshal
from pathlib import Path

import numpy as np

from substrate_framework.partitioned_q_work import QWork
from substrate_framework.partitioned_q_storage import bounded_array_identity, release_mapped_pages
from substrate_framework.streamed_joint_response import MappedResponse, NAMES
from substrate_framework.streamed_material_state import authenticate, file_binding


OPERATION = 'committed-selected-weighted-complex-readout-v1'


def _readonly(value):
    result = np.asarray(value).view()
    result.flags.writeable = False
    return result


def _callback_binding(callback):
    """Bind executable callback bytes; construction authority is bound separately."""
    function = getattr(callback, '__func__', callback)
    code = getattr(function, '__code__', None)
    source = inspect.getsourcefile(function if code is not None else type(function))
    result = dict(module=getattr(function, '__module__', type(function).__module__),
                  name=getattr(function, '__qualname__', type(function).__qualname__))
    if source is not None and Path(source).is_file():
        result['source'] = file_binding(source)
    if code is not None:
        result['code_sha256'] = sha256(marshal.dumps(code)).hexdigest()
    return result


def _committed_state(executor, state, time, binding):
    if not isinstance(state, MappedResponse) or state.response is not executor.r:
        raise TypeError('restored MappedResponse for this exact response required')
    record = json.loads(authenticate(binding).read_text())
    if (record.get('schema') != 'streamed-joint-response-v1'
            or record.get('complete') is not True or record.get('kind') != 'forward'
            or record.get('time') != float(time)
            or record.get('q') != executor.r.q_fractional.tolist()
            or not isinstance(record.get('identity'), str) or not record['identity']):
        raise ValueError('complete physical forward checkpoint at the exact selected time/q required')
    if Path(binding['path']).resolve().parent != state.path:
        raise ValueError('state is not mapped from the explicitly authenticated checkpoint')
    if len(record['legs']) != 2*len(executor.r.edges):
        raise ValueError('checkpoint must contain every original independent ket/bra leg')
    for index, leg_binding in enumerate(record['legs']):
        executor.budget.check()
        path = authenticate(leg_binding)
        value = (state.ket if index % 2 == 0 else state.bra)[index//2]
        expected = state.path/('leg-%04d.npy' % index)
        if (Path(path).resolve() != expected or not isinstance(value, np.memmap)
                or Path(value.filename).resolve() != expected
                or value.flags.writeable or value.mode != 'r'
                or value.dtype != np.complex128):
            raise ValueError('state leg is not the original authenticated readonly checkpoint map')
        source = executor.r.edges[index//2]
        ns, rs = executor.r.sectors[source.source].shape
        nt, rt = executor.r.sectors[source.target].shape
        if value.shape != ((nt, rs) if index % 2 == 0 else (rt, ns)):
            raise ValueError('committed original rectangular leg shape differs')
    for name in NAMES:
        expected = np.asarray([complex(*z) for z in record['small'][name]], dtype=np.complex128)
        actual = getattr(state, name)
        if actual.shape != expected.shape or not np.array_equal(actual, expected):
            raise ValueError('state canonical component differs from authenticated checkpoint: '+name)
    return record


def _edge_scalar(executor, state, probe, time, index):
    """Contract one source, then release borrowed resident pages, not its data."""
    executor.budget.check()
    edge = executor.r.edges[index]
    source_frame = target_frame = None
    try:
        source_frame = executor.edge_stage.frame(edge.source, time)
        target_frame = executor.edge_stage.frame(edge.target, time)
        term = executor.r.sectors[edge.source].weight*(
            np.vdot(probe.actions[index](_readonly(source_frame)), state.ket[index])
            + np.einsum('ij,ji->', state.bra[index],
                        probe.adjoints[index](_readonly(target_frame))))
        return dict(value=np.asarray([term], dtype=np.complex128))
    finally:
        # ndarray views hide the memmap type; retain the original frame handles.
        for value in (state.ket[index], state.bra[index], source_frame, target_frame):
            release_mapped_pages(value)


def partitioned_selected_readout(executor, state, probe, time, work, *, state_binding) -> complex:
    """Evaluate the original weighted complex formula in original index sum order.

    ``state_binding`` authenticates a complete forward endpoint, not mutable RK4
    scratch. Required ``executor.readout_source_bindings`` are authenticated
    before cache lookup. Required ``executor.readout_operator_binding`` binds
    exact probe-construction inputs: a label alone is not an operator identity.
    Neither authority relabels the historical producer as the current consumer.
    Completed scalar receipts and original worker costs remain in ``work.last``.
    """
    if not isinstance(work, QWork) or work.budget is not executor.budget:
        raise TypeError('genuine pinned QWork sharing the executor QBudget required')
    if not np.isfinite(time) or time < 0:
        raise ValueError('explicit finite nonnegative committed readout time required')
    probe.validate(executor.r)
    record = _committed_state(executor, state, time, state_binding)
    sources = list(getattr(executor, 'readout_source_bindings', ()))
    construction = getattr(executor, 'readout_operator_binding', None)
    if not sources or not isinstance(construction, dict) or not construction:
        raise ValueError('authenticated current readout_source_bindings and exact nonempty readout_operator_binding required')
    for source in sources:
        authenticate(source)
    if hasattr(executor.edge_stage, 'clear_interval'):
        executor.edge_stage.clear_interval()
    if hasattr(executor, 'prepare_time'):
        executor.prepare_time(time)
    r = executor.r
    identity = dict(checkpoint=dict(state_binding), producer_identity=record['identity'],
        consumer_identity=executor.identity, work_identity=work.identity,
        sources=sources, implementation=file_binding(__file__), time=float(time),
        stage=executor.edge_stage.stage_identity(time), q=r.q_fractional.tolist(),
        scalar_identity=r.scalar_identity, coverage_identity=r.coverage_identity,
        sqrt_n=float(r.sqrt_n),
        sectors=[dict(shape=list(s.shape), d=bounded_array_identity(s.d),
                      weight=float(s.weight), component=s.density_component) for s in r.sectors],
        edges=[dict(source=e.source, target=e.target) for e in r.edges],
        probe=dict(identity=probe.identity, units=probe.units,
            canonical={name: bounded_array_identity(getattr(probe, name)) for name in NAMES},
            actions=[_callback_binding(a) for a in probe.actions],
            adjoints=[_callback_binding(a) for a in probe.adjoints],
            construction=construction))
    # The small sum is computed ONCE, before the same sequential += edge order
    # as the original implementation. The manager gathers only scalar arrays.
    value = sum(np.vdot(getattr(probe, name), getattr(state, name)) for name in NAMES)

    receipts = work.execute(OPERATION, identity, len(r.edges),
        lambda index: _edge_scalar(executor, state, probe, time, index))
    for index, receipt in enumerate(receipts):
        executor.budget.check()
        if receipt['index'] != index or set(receipt['arrays']) != {'value'}:
            raise ValueError('one actual complex scalar receipt per original edge required')
        scalar = np.load(authenticate(receipt['arrays']['value']), allow_pickle=False)
        if scalar.shape != (1,) or scalar.dtype != np.complex128 or not np.isfinite(scalar).all():
            raise ValueError('edge readout must return only finite complex128 value shape (1,)')
        value += scalar[0]
    if not np.isfinite(value):
        raise ValueError('nonfinite selected complex readout')
    return complex(value)
