"""Bounded full-row classical RK4; completed boundaries only, never a repair.

Covariance labels denote INITIAL origins. Streaming retains one boundary, one
stage, one RHS and one weighted increment (no full trajectory or rank projection).
Checkpoints restore completed boundaries, not an interrupted stage or a process.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import time
import numpy as np

from .prepared_factor_dynamics import FactorMean, FactorTangent
from .prepared_covariance_source import OrderedSource

PARTS = ('electronic', 'ionic', 'cross')


@dataclass(frozen=True)
class TimeLimits:
    max_steps: int
    max_mean_calls: int
    max_tangent_calls: int
    max_seconds: float

    def __post_init__(self):
        for value in (self.max_steps, self.max_mean_calls, self.max_tangent_calls):
            if type(value) is not int or value < 0:
                raise ValueError('nonnegative integer limits required')
        if not np.isfinite(self.max_seconds) or self.max_seconds <= 0:
            raise ValueError('positive finite elapsed-time limit required')


@dataclass(frozen=True)
class RealityTolerance:
    atol: float
    rtol: float

    def __post_init__(self):
        if not all(np.isfinite(x) and x >= 0 for x in (self.atol, self.rtol)):
            raise ValueError('explicit nonnegative finite reality allowances required')


@dataclass(frozen=True)
class TimeBoundary:
    mean: FactorMean
    sources: tuple[OrderedSource, ...]
    time: float
    completed_steps: int


class TimeLimitExceeded(RuntimeError):
    pass


class MeanRealityError(ValueError):
    pass


def _combine(base, increment, scale):
    if isinstance(base, FactorMean):
        return FactorMean(tuple(a + scale*b for a, b in zip(base.frames, increment.frames)),
                          base.u + scale*increment.u, base.p + scale*increment.p)
    return FactorTangent(tuple(a + scale*b for a, b in zip(base.ket, increment.ket)),
                         tuple(a + scale*b for a, b in zip(base.bra, increment.bra)),
                         base.u + scale*increment.u, base.p + scale*increment.p)


def _scale(value, scale):
    if isinstance(value, FactorMean):
        return FactorMean(tuple(scale*x for x in value.frames), scale*value.u, scale*value.p)
    return FactorTangent(tuple(scale*x for x in value.ket), tuple(scale*x for x in value.bra),
                         scale*value.u, scale*value.p)


def _source_map(sources, function):
    return tuple(replace(s, **{p: None if getattr(s, p) is None else function(i, p, getattr(s, p))
                              for p in PARTS}) for i, s in enumerate(sources))


def _finite(value):
    arrays = value.frames if isinstance(value, FactorMean) else value.ket + value.bra
    return all(np.all(np.isfinite(x)) for x in (*arrays, value.u, value.p))


def _encode(array):
    a = np.asarray(array)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def _decode(record):
    return np.asarray(record['real']) + 1j*np.asarray(record['imag'])


def _pack(value):
    fields = ('frames',) if isinstance(value, FactorMean) else ('ket', 'bra')
    return {**{f: [_encode(x) for x in getattr(value, f)] for f in fields},
            'u': _encode(value.u), 'p': _encode(value.p)}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


class PreparedTimeEvolution:
    """Single-owner RK4 stream. Identity binds model/preparation/coordinates/sources.

    The owner supplies nonempty source-byte/model/preparation/coordinate/budget
    identities; restore requires exactly the same identity, limits and tolerance.
    Calls are counted immediately before invocation, including calls that raise.
    All four stages are reserved atomically; failed reservations spend nothing,
    failed stages retain their reservations/executions and the last boundary.
    Do not mutate arrays exposed in a boundary while the stream is in use.
    """
    def __init__(self, model, mean, sources, *, identity, limits, reality, initial_time=0.):
        required = ('model', 'preparation', 'coordinates', 'sources', 'budget')
        if any(not isinstance(identity.get(k), str) or not identity[k].strip() for k in required):
            raise ValueError('explicit model/preparation/coordinates/sources/budget identities required')
        if identity['coordinates'] != model.reference.coordinate_identity:
            raise ValueError('coordinate identity mismatch')
        if not np.isfinite(initial_time):
            raise ValueError('finite initial time required')
        model.mean(mean.frames, mean.u, mean.p)
        sources = tuple(sources)
        for source in sources:
            if source.ordering not in ('AB', 'BA') or not source.preparation_authority or not source.right_identity:
                raise ValueError('ordered source authority required')
            if source.electronic is None:
                raise ValueError('electronic source cannot be unknown')
            for part in PARTS:
                x = getattr(source, part)
                if x is not None:
                    model.tangent(x.ket, x.bra, x.u, x.p)
        self.model, self.identity = model, dict(identity)
        self.limits, self.reality = limits, reality
        self.boundary = TimeBoundary(mean, sources, float(initial_time), 0)
        self.counts = dict(reserved_mean=0, executed_mean=0, reserved_tangent=0, executed_tangent=0)
        self._carried_seconds = 0.
        self._opened_at = time.monotonic()
        self.reality_diagnostics = dict(max_imaginary=0., max_allowance_ratio=0., stages_checked=0,
                                        last_raw_u_imag=[], last_raw_p_imag=[])
        self.failed = False

    @property
    def elapsed_seconds(self):
        """Cumulative live-owner wall time, including between yielded steps."""
        return self._carried_seconds + time.monotonic() - self._opened_at

    def _deadline(self):
        if self.elapsed_seconds > self.limits.max_seconds:
            raise TimeLimitExceeded('elapsed time exhausted; completed boundary retained')

    def _reserve(self):
        active = sum(getattr(s, p) is not None for s in self.boundary.sources for p in PARTS)
        c, l = self.counts, self.limits
        if (self.boundary.completed_steps >= l.max_steps or self.elapsed_seconds >= l.max_seconds
                or c['reserved_mean'] + 4 > l.max_mean_calls
                or c['reserved_tangent'] + 4*active > l.max_tangent_calls):
            raise TimeLimitExceeded('entire RK4 step cannot be reserved')
        c['reserved_mean'] += 4
        c['reserved_tangent'] += 4*active

    def _real_rhs(self, raw):
        d = self.reality_diagnostics
        d['stages_checked'] += 1
        for name in ('u', 'p'):
            a = np.asarray(getattr(raw, name))
            residual = np.abs(a.imag)
            allowance = self.reality.atol + self.reality.rtol*np.abs(a.real)
            ratio = np.divide(residual, allowance, out=np.where(residual == 0, 0., np.inf), where=allowance > 0)
            d['last_raw_' + name + '_imag'] = [float(x) if np.isfinite(x) else str(x) for x in a.imag]
            peak = float(np.max(residual, initial=0))
            ratio_peak = float(np.max(ratio, initial=0))
            for key, value in (('max_imaginary', peak), ('max_allowance_ratio', ratio_peak)):
                if not np.isfinite(value):
                    d[key] = None  # raw nonfinite/zero-allowance excess retained above
                elif d[key] is not None:
                    d[key] = max(d[key], value)
            if not np.all(np.isfinite(a)) or np.any(residual > allowance):
                raise MeanRealityError('classical mean derivative exceeds registered real-codomain allowance')
        return FactorMean(raw.frames, raw.u.real, raw.p.real)

    def _rhs(self, mean, sources):
        self._deadline()
        self.counts['executed_mean'] += 1
        evaluation = self.model.evaluate_mean(mean)
        dm = self._real_rhs(evaluation.rhs)

        def tangent(i, part, value):
            self._deadline()
            self.counts['executed_tangent'] += 1
            return self.model.tangent_rhs(evaluation, value)

        ds = _source_map(sources, tangent)
        if not _finite(dm) or any(not _finite(getattr(s, p)) for s in ds for p in PARTS if getattr(s, p) is not None):
            raise ValueError('nonfinite stage RHS; boundary retained')
        return dm, ds

    def step(self, dt):
        if self.failed:
            raise RuntimeError('failed stream requires explicit completed-boundary restoration')
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('positive finite step required')
        self._reserve()
        base = self.boundary
        mean, sources = base.mean, base.sources
        try:
            for stage, (weight, offset) in enumerate(((1/6, .5), (1/3, .5), (1/3, 1.), (1/6, 0.))):
                dm, ds = self._rhs(mean, sources)
                if stage == 0:
                    am = _scale(dm, weight)
                    accum = _source_map(ds, lambda i, p, x: _scale(x, weight))
                else:
                    am = _combine(am, dm, weight)
                    accum = _source_map(accum, lambda i, p, x: _combine(x, getattr(ds[i], p), weight))
                if stage < 3:
                    mean = _combine(base.mean, dm, dt*offset)
                    sources = _source_map(base.sources, lambda i, p, x: _combine(x, getattr(ds[i], p), dt*offset))
            mean = _combine(base.mean, am, dt)
            sources = _source_map(base.sources, lambda i, p, x: _combine(x, getattr(accum[i], p), dt))
            self._deadline()
            if not _finite(mean) or any(not _finite(getattr(s, p)) for s in sources for p in PARTS if getattr(s, p) is not None):
                raise ValueError('nonfinite completed candidate; boundary retained')
            self.boundary = TimeBoundary(mean, sources, base.time + dt, base.completed_steps + 1)
            return self.boundary
        except Exception:
            self.failed = True
            raise

    def stream(self, *, dt, steps):
        if type(steps) is not int or steps < 0:
            raise ValueError('nonnegative integer step count required')
        for _ in range(steps):
            yield self.step(dt)

    def checkpoint(self):
        """JSON-compatible complete boundary and cumulative resource consumption."""
        b = self.boundary
        sources = [{**{p: None if getattr(s, p) is None else _pack(getattr(s, p)) for p in PARTS},
                    **{k: getattr(s, k) for k in ('ordering', 'right_identity', 'preparation_authority')}} for s in b.sources]
        payload = dict(schema='prepared-time-boundary/1', identity=self.identity, limits=asdict(self.limits),
                       reality=asdict(self.reality), time=b.time, completed_steps=b.completed_steps,
                       mean=_pack(b.mean), sources=sources, counts=dict(self.counts),
                       elapsed_seconds=self.elapsed_seconds, reality_diagnostics=self.reality_diagnostics,
                       failed_attempt=self.failed)
        # Round-trip detaches mutable arrays/lists from the live owner.
        payload = json.loads(json.dumps(payload, allow_nan=False))
        return {'payload': payload, 'sha256': _digest(payload)}

    @classmethod
    def restore(cls, model, checkpoint, *, identity, limits, reality):
        p = checkpoint['payload']
        if checkpoint['sha256'] != _digest(p) or p['schema'] != 'prepared-time-boundary/1':
            raise ValueError('checkpoint content identity mismatch')
        if p['identity'] != identity or p['limits'] != asdict(limits) or p['reality'] != asdict(reality):
            raise ValueError('restore requires identical source/model/budget contract')
        m = p['mean']
        mean = model.mean(tuple(_decode(x) for x in m['frames']), _decode(m['u']), _decode(m['p']))
        def unpack(x):
            if x is None:
                return None
            return model.tangent(tuple(_decode(a) for a in x['ket']), tuple(_decode(a) for a in x['bra']),
                                 _decode(x['u']), _decode(x['p']))
        sources = tuple(OrderedSource(*(unpack(s[k]) for k in PARTS), s['ordering'],
                                      s['right_identity'], s['preparation_authority']) for s in p['sources'])
        result = cls(model, mean, sources, identity=identity, limits=limits, reality=reality, initial_time=p['time'])
        count = p['completed_steps']
        active = sum(getattr(s, k) is not None for s in sources for k in PARTS)
        c = p['counts']
        if (type(count) is not int or not 0 <= count <= limits.max_steps
                or any(type(v) is not int or v < 0 for v in c.values())
                or not 4*count <= c['executed_mean'] <= c['reserved_mean'] <= limits.max_mean_calls
                or not 4*active*count <= c['executed_tangent'] <= c['reserved_tangent'] <= limits.max_tangent_calls
                or not np.isfinite(p['elapsed_seconds']) or p['elapsed_seconds'] < 0):
            raise ValueError('invalid cumulative budget state')
        result.boundary = replace(result.boundary, completed_steps=count)
        result.counts = dict(c)
        result._carried_seconds = p['elapsed_seconds']
        result.reality_diagnostics = dict(p['reality_diagnostics'])
        return result
