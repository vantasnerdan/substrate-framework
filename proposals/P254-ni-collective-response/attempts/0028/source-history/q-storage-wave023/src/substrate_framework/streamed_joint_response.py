"""Disk-backed all-edge RK4 tangent and its actual Hermitian dual.

The dual is the adjoint of the *discrete* RK4 map, not a separately discretized
backward ODE. Each reverse stage uses the same moving baseline stage and full
reduction barrier as its forward counterpart. No PW-square objects are formed.
Scratch files are reusable; completed checkpoints are immutable and hashed.
"""
from contextlib import contextmanager
from pathlib import Path
import json

import numpy as np

from .streamed_material_state import authenticate, file_binding


NAMES = ('u', 'p', 'coordinate', 'momentum')


class MappedResponse:
    def __init__(self, response, directory, budget, *, create=False):
        self.response, self.path = response, Path(directory).resolve()
        if create:
            self.path.mkdir(parents=True, exist_ok=False)
        self.ket, self.bra = [], []
        shapes = []
        for e in response.edges:
            ns, rs = response.sectors[e.source].shape
            nt, rt = response.sectors[e.target].shape
            shapes.extend(((nt, rs), (rt, ns)))
        for i, shape in enumerate(shapes):
            path = self.path/('leg-%04d.npy' % i)
            if create:
                budget.reserve(16*int(np.prod(shape))+512)
                value = np.lib.format.open_memmap(path, mode='w+', dtype='complex128', shape=shape)
                value[:] = 0
            else:
                value = np.load(path, mmap_mode='r', allow_pickle=False)
                if value.shape != shape or value.dtype != np.complex128:
                    raise ValueError('checkpoint original rectangular leg dimensions differ')
            (self.ket if i % 2 == 0 else self.bra).append(value)
        for name, size in zip(NAMES, (response.ndof, response.ndof, response.nfield, response.nfield)):
            setattr(self, name, np.zeros(size, complex))

    def zero(self):
        for value in self.values():
            value[:] = 0
        return self

    def values(self):
        yield from self.ket
        yield from self.bra
        for name in NAMES:
            yield getattr(self, name)

    def combine(self, *terms):
        for index, target in enumerate(self.values()):
            # One edge-sized temporary, including when target aliases an input.
            value = sum(a*tuple(x.values())[index] for a, x in terms)
            target[:] = value
        return self

    def dot(self, other):
        return sum(np.vdot(a, b) for a, b in zip(self.values(), other.values()))

    def publish(self, *, time, identity, budget, kind, provenance=None):
        if not np.isfinite(time):
            raise ValueError('finite checkpoint time required')
        for value in self.values():
            budget.check()
            if not np.isfinite(value).all():
                raise ValueError('nonfinite response is not an adoptable checkpoint')
        for value in self.ket+self.bra:
            value.flush()
        record = dict(schema='streamed-joint-response-v1', complete=True, time=float(time),
                      identity=identity, kind=kind, provenance=provenance, q=self.response.q_fractional.tolist(),
                      legs=[file_binding(self.path/('leg-%04d.npy' % i))
                            for i in range(2*len(self.ket))],
                      small={name: [[float(z.real), float(z.imag)] for z in getattr(self, name)] for name in NAMES})
        return budget.json(self.path/'checkpoint.json', record)

    @classmethod
    def restore(cls, response, binding, budget, identity):
        record = json.loads(authenticate(binding).read_text())
        if (record['schema'] != 'streamed-joint-response-v1' or not record['complete']
                or record['identity'] != identity or record['q'] != response.q_fractional.tolist()):
            raise ValueError('response checkpoint identity/q differs')
        for leg in record['legs']:
            authenticate(leg)
        result = cls(response, Path(binding['path']).parent, budget)
        for name in NAMES:
            setattr(result, name, np.asarray([complex(*z) for z in record['small'][name]]))
        return result, record


def resource_plan(response):
    sizes = [16*(response.sectors[e.target].shape[0]*response.sectors[e.source].shape[1]
                  + response.sectors[e.target].shape[1]*response.sectors[e.source].shape[0]) for e in response.edges]
    one = sum(sizes)
    return dict(one_q_leg_bytes=one, scratch_states=2, scratch_leg_bytes=2*one,
                checkpoint_leg_bytes=one, largest_edge_leg_bytes=max(sizes),
                policy='one q at a time; two reusable mapped scratch states; immutable selected checkpoints',
                rss='edge-local arrays/operators and grids; mmap page cache is reclaimable',
                no_dense_pw_square=True, no_all_q_allocation=True)

class StreamedJointResponse:
    """edge_stage(i,t) yields current ps,pt,Hs,Hs†,Ht,Ht†,g_ts,L_xyz.

    K_adjoint is the actual grid-Euclidean adjoint. The physical uniform grid
    measure is supplied explicitly, not inferred from scalar equality.
    """
    def __init__(self, response, edge_stage, *, kernel_adjoint, grid_measure, budget,
                 directory, identity, dynamics='screened-joint'):
        if dynamics not in ('screened-joint', 'bare-electronic') or grid_measure <= 0:
            raise ValueError('explicit physical measure and dynamics label required')
        self.r, self.edge_stage, self.ka = response, edge_stage, kernel_adjoint
        self.measure, self.budget, self.identity, self.dynamics = grid_measure, budget, identity, dynamics
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=False)
        self.work = [MappedResponse(response, self.directory/('scratch-%d' % i), budget, create=True) for i in range(2)]

    def _field_moments(self, axes, left, right):
        return np.asarray([np.vdot(a.apply(right), left) for a in axes])

    def rhs(self, time, x, out, *, dual=False):
        r = self.r
        n, g, ell = np.zeros(r.density_shape, complex), np.zeros(r.ndof, complex), np.zeros(3, complex)
        coupled = self.dynamics == 'screened-joint'
        # Complete reduction BEFORE K or any feedback. Dual Y is the adjoint
        # electronic drive, with D exactly once and NO original weight here.
        for i, e in enumerate(r.edges):
            self.budget.check()
            with self.edge_stage(i, time) as (ps, pt, hs, hsa, ht, hta, c, axes):
                if dual:
                    left, right = 1j*x.ket[i] @ r.sectors[e.source].d, ps
                    left2, right2 = -1j*pt @ r.sectors[e.target].d, x.bra[i].conj().T
                    weight = 1.
                else:
                    left, right, left2, right2 = x.ket[i], ps, pt, x.bra[i].conj().T
                    weight = r.sectors[e.source].weight
                if coupled:
                    n1, g1 = c.moments(left, right, weight=weight)
                    n2, g2 = c.moments(left2, right2, weight=weight)
                    n += n1+n2
                    g += g1+g2
                    if r.gamma:
                        ell += weight*(self._field_moments(axes, left, right)+self._field_moments(axes, left2, right2))
        potential = (self.ka(n, q_bohr_inv=r.q_bohr_inv) if dual else r.kernel(n, q_bohr_inv=r.q_bohr_inv)) if coupled else np.zeros_like(n)
        for i, e in enumerate(r.edges):
            self.budget.check()
            s, t = r.sectors[e.source], r.sectors[e.target]
            with self.edge_stage(i, time) as (ps, pt, hs, hsa, ht, hta, c, axes):
                if dual:
                    direction = -x.p if coupled else np.zeros(r.ndof)
                    field = -r.sqrt_n*r.polarizations.T @ x.momentum if r.gamma and coupled else np.zeros(3)
                else:
                    direction = x.u if coupled else np.zeros(r.ndof)
                    field = r.sqrt_n*r.polarizations.T @ x.coordinate if r.gamma and coupled else np.zeros(3)
                op = c.prepare(direction.reshape(-1, 3))
                v = potential[:, s.density_component]
                drive = c.bridge.local_apply(v, ps)+op.apply(ps)
                adj = c.bridge.local_adjoint_apply(v, pt)+op.adjoint_apply(pt)
                for coefficient, a in zip(field, axes):
                    drive += coefficient*a.apply(ps)
                    adj += coefficient.conjugate()*a.adjoint_apply(pt)
                if dual:
                    out.ket[i][:] = 1j*hta(x.ket[i])+s.weight*drive
                    out.bra[i][:] = -1j*hs(x.bra[i].conj().T).conj().T+s.weight*adj.conj().T
                else:
                    out.ket[i][:] = -1j*(ht(x.ket[i])+drive @ s.d)
                    out.bra[i][:] = 1j*(hsa(x.bra[i].conj().T).conj().T+t.d @ adj.conj().T)
        return self.small_rhs(x, out, g, ell, dual=dual)

    def small_rhs(self, x, out, g, ell, *, dual=False):
        r = self.r
        coupled = self.dynamics == 'screened-joint'
        for name in NAMES:
            getattr(out, name)[:] = 0
        if not coupled:
            return out
        if dual:
            out.u[:] = g-r.curvature.conj().T @ x.p
            out.p[:] = x.u/r.masses
            if r.gamma:
                out.u[:] -= r.sqrt_n*r.mixed.conj() @ x.momentum
                out.p[:] += r.sqrt_n*(r.charge @ x.momentum)/r.masses
                out.coordinate[:] = (r.sqrt_n*(r.polarizations @ ell-r.charge.T @ (x.u/r.masses)-r.mixed.conj().T @ x.p)
                    -r.ncells*(r.contact+r.charge.T @ (r.charge/r.masses[:, None])).conj().T @ x.momentum)
                out.momentum[:] = x.coordinate/r.field_mass
        else:
            out.u[:] = x.p/r.masses
            out.p[:] = -r.curvature @ x.u-g
            if r.gamma:
                out.u[:] -= r.sqrt_n*(r.charge @ x.coordinate)/r.masses
                out.p[:] -= r.sqrt_n*r.mixed @ x.coordinate
                out.coordinate[:] = x.momentum/r.field_mass
                out.momentum[:] = (-r.sqrt_n*(r.polarizations @ ell+r.mixed.T @ x.u-r.charge.T @ (x.p/r.masses))
                    -r.ncells*(r.contact+r.charge.T @ (r.charge/r.masses[:, None])) @ x.coordinate)
        return out

    def step(self, x, out, time, dt, *, dual=False):
        """Exact transpose of RK4 computational graph for the supplied history."""
        a, b = self.work
        if hasattr(self.edge_stage, 'set_interval'):
            self.edge_stage.set_interval(time, dt)
        if not dual:
            self.rhs(time, x, a)
            out.combine((1, x), (dt/6, a))
            b.combine((1, x), (dt/2, a))
            self.rhs(time+dt/2, b, a)
            out.combine((1, out), (dt/3, a))
            b.combine((1, x), (dt/2, a))
            self.rhs(time+dt/2, b, a)
            out.combine((1, out), (dt/3, a))
            b.combine((1, x), (dt, a))
            self.rhs(time+dt, b, a)
            out.combine((1, out), (dt/6, a))
        else:
            # k4bar=h*x/6; z4bar=A4†k4bar; k3bar=h*x/3+h*z4bar.
            a.combine((dt/6, x))
            self.rhs(time+dt, a, b, dual=True)
            out.combine((1, x), (1, b))
            a.combine((dt/3, x), (dt, b))
            self.rhs(time+dt/2, a, b, dual=True)
            out.combine((1, out), (1, b))
            a.combine((dt/3, x), (dt/2, b))
            self.rhs(time+dt/2, a, b, dual=True)
            out.combine((1, out), (1, b))
            a.combine((dt/6, x), (dt/2, b))
            self.rhs(time, a, b, dual=True)
            out.combine((1, out), (1, b))
        return out

    def terminal_probe(self, time, probe, out):
        probe.validate(self.r)
        for i, e in enumerate(self.r.edges):
            with self.edge_stage(i, time) as (ps, pt, *unused):
                w = self.r.sectors[e.source].weight
                out.ket[i][:] = w*probe.actions[i](ps)
                out.bra[i][:] = w*probe.adjoints[i](pt).conj().T
        for name in NAMES:
            getattr(out, name)[:] = getattr(probe, name)
        return out

    def propagate(self, initial, times, directory, *, dual=False):
        times = np.asarray(times, float)
        if len(times) < 2 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError('strictly increasing explicit history-covered integration grid required')
        current = initial
        indices = range(len(times)-2, -1, -1) if dual else range(len(times)-1)
        records = []
        for index in indices:
            out = MappedResponse(self.r, Path(directory)/('step-%06d' % index), self.budget, create=True)
            self.step(current, out, times[index], times[index+1]-times[index], dual=dual)
            records.append(out.publish(time=times[index if dual else index+1], identity=self.identity,
                                       budget=self.budget, kind='dual' if dual else 'forward',
                                       provenance=getattr(self, 'provenance', None)))
            current = out
        return current, records
