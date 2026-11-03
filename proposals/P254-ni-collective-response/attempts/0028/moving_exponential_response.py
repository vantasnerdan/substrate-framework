"""H0 interaction-picture midpoint for the actual moving JointBlochResponse.

Algebra (time is inverse Ry)
---------------------------
For each directed edge, X=L Psi_s^dagger+Psi_t R, with independent complex
L[nt,rs], R[rt,ns]. Split y'=A0 y+B(t)y, where
    A0 L=-i H0_t L,  A0 R=+i R H0_s,
and A0 is zero on (U,P,a,pi). B contains the FULL coupled density/ionic/field
feedback and CURRENT H(t)-H0, not a stationary or clamped physical model.
Writing z(s)=exp(-s A0)y(t+s), explicit midpoint on this interaction-picture
ODE gives the genuine exponential method
    v = E(h/2)(y + h/2 B(t)y),
    y_next = E(h/2)(E(h/2)y + h B(t+h/2)v).
E acts on full PW rows, including the empty complement; D remains the raw
original matrix and occurs in the response drive exactly once. Expanding gives
I+h(A0+B)+h^2/2*((A0+B)^2+Bdot)+O(h^3). No explicit polynomial is used for H0.
This is second order for smooth finite-dimensional moving stages. It removes
RK's free-H0 stability restriction, NOT all stiffness: residual/feedback scales
and commutators still constrain temporal accuracy and may cause order reduction.

Consumer contract
-----------------
JointResidual(response, split_stage, free_pair, prepare_interval=None) uses the
existing JointStage/JointTangent physics. split_stage(t) returns SplitJointStage
with the actual CURRENT JointStage and independently formed H(t)-H0 actions and
adjoints. Do not subtract large H0 actions numerically, freeze t0, phase-reset
frames, infer a real density, or replace complex perturbations by Hermitian ones.
The supplied mean trajectory must cover every requested time; piecewise saved
histories must supply prepare_interval(t,h) to reject steps crossing a saved mean
boundary. Refine the mean separately: response refinement on an unchanged saved
mean does not measure mean error, nor license extrapolation beyond its endpoint.

Use SourceLocalMidpoint(model).evolve(initial, times). Its result.state exposes
ket/bra indexed sequences and canonical arrays, compatible with moments/readout.
Initial can be JointTangent or SourceLocalTangent with an immutable edge getter.
There is no full-Q state allocation: each interval stores only two grid potentials
and small canonical arrays; electronic edges are replayed in sequence from the
original source. Replaying K completed intervals for each of K reductions costs
O(K^2) edge work; each replayed interval uses 3 free half-flow applications and
2 residual edge actions (each free flow acts independently on both legs). This
is a concrete storage/work trade, not a speedup claim. Native execution may need
a different checkpoint/replay allocation. No persistent files are written here.
The resident JointResidual keeps those grid contexts in memory. A Native generic
model must store/reload its immutable potentials within the actual grid/disk
allocation; source-local does NOT make O(steps*grid) payload disappear. The
context_bytes property counts numerical context/canonical payload only, not
Python metadata, provider caches, initial storage, or exponential workspace.

Free-flow contract: JointH0Flow takes endpoint callbacks implementing
exp(-i*tau*H0) and exp(-i*tau*H0^dagger) on arbitrary rectangular full-row blocks.
The latter propagates R^dagger, NOT the ket and NOT a conjugated tangent drive.
Each callback must budget its actual action error and H0 applications. Their
errors, stage interpolation error, replay roundoff, and finite-step error remain
separate; compare physical observables on h,h/2,h/4 and tighten free-action/mean
errors. refinement_difference is an EMPIRICAL second-order estimate, not a bound.

Native gap (not hidden by this finite-system API)
------------------------------------------------
MaterialJointHistory.__call__ yields (ps,pt,Hs,Hs_dag,Ht,Ht_dag,c,axes), not
split H-H0 actions or arbitrary-leg H0 exponentials. Its private _frame only
propagates prior mean columns under a saved FULL frozen H. Native needs budgeted
edge-local residual actions formed directly from saved potential/u/field plus
H0 exponential callbacks for BOTH legs, a lazy current JointStage provider,
and its existing complete-reduction barriers. SourceLocalQResponse's saved
contexts/replay are RK4-specific and cannot consume these contexts unchanged.
The existing exponential_action can inform an implementation, but this module
installs NO Native adapter and authorizes NO Native job. Native wall time,
workspace, disk, action tolerance and cost must be measured before allocation.
"""
from dataclasses import asdict, dataclass
from typing import Callable

import numpy as np

from substrate_framework.joint_bloch_response import JointStage


NAMES = ('u', 'p', 'coordinate', 'momentum')


@dataclass
class MethodCounts:
    context_reductions: int = 0
    residual_edge_actions: int = 0
    free_pair_actions: int = 0
    replayed_intervals: int = 0


class _CountedModel:
    """Operation counts, not Native H0 application counts or runtime forecasts."""
    def __init__(self, model):
        self.model, self.counts = model, MethodCounts()
        self.response = model.response
        self.prepare_interval = getattr(model, 'prepare_interval', None)

    def context(self, time, state):
        self.counts.context_reductions += 1
        return self.model.context(time, state)

    def edge_rhs(self, context, index, pair, small):
        self.counts.residual_edge_actions += 1
        return self.model.edge_rhs(context, index, pair, small)

    def free_pair(self, index, pair, tau):
        self.counts.free_pair_actions += 1
        return self.model.free_pair(index, pair, tau)


@dataclass(frozen=True)
class SmallState:
    u: np.ndarray
    p: np.ndarray
    coordinate: np.ndarray
    momentum: np.ndarray

    @classmethod
    def from_state(cls, state):
        return cls(*(np.asarray(getattr(state, name), dtype=complex) for name in NAMES))


def small_combination(*terms):
    return SmallState(*(sum(coefficient*getattr(state, name)
                            for coefficient, state in terms) for name in NAMES))


def pair_combination(*terms):
    return tuple(sum(coefficient*pair[leg] for coefficient, pair in terms)
                 for leg in (0, 1))


class _Legs:
    def __init__(self, state, leg):
        self.state, self.leg = state, leg

    def __len__(self):
        return self.state.count

    def __getitem__(self, index):
        if not 0 <= index < self.state.count:
            raise IndexError(index)
        return self.state.edge(index)[self.leg]


class SourceLocalTangent:
    """Immutable-value view with at most one cached electronic edge per view.

    edge(index) must return full-PW independent (ket,bra) arrays without mutating
    any input. Retain original input bytes for the entire replay lifetime.
    """
    def __init__(self, count, edge: Callable, small):
        self.count, self._edge, self.small = count, edge, SmallState.from_state(small)
        self.ket, self.bra = _Legs(self, 0), _Legs(self, 1)
        self._index, self._pair = None, None
        for name in NAMES:
            setattr(self, name, getattr(self.small, name))

    def edge(self, index):
        if not 0 <= index < self.count:
            raise IndexError(index)
        if self._index != index:
            self._index, self._pair = None, None
            self._pair = self._edge(index)
            self._index = index
        return self._pair

    def release(self):
        self._index, self._pair = None, None

    @classmethod
    def from_joint(cls, tangent):
        return cls(len(tangent.ket), lambda i: (tangent.ket[i], tangent.bra[i]), tangent)


@dataclass(frozen=True)
class SplitJointStage:
    current: JointStage
    residual_actions: tuple
    residual_adjoints: tuple


class JointH0Flow:
    """Exact/numerically controlled H0 flow; no saved-band projection permitted."""
    def __init__(self, response, propagate, propagate_adjoint):
        self.response = response
        self.propagate, self.propagate_adjoint = propagate, propagate_adjoint

    def __call__(self, index, pair, tau):
        edge = self.response.edges[index]
        left, right = pair
        return (self.propagate(edge.target, left, tau),
                self.propagate_adjoint(edge.source, right.conj().T, tau).conj().T)


@dataclass(frozen=True)
class FeedbackContext:
    time: float
    potential: np.ndarray
    small_rhs: SmallState


class JointResidual:
    """Two-pass B(t) using original JointBlochResponse moments and scalar.

    The first pass completes ALL edges before K; edge_rhs uses that immutable
    complete context. Only the H-H0 connection is split off; feedback is intact.
    """
    def __init__(self, response, split_stage, free_pair, *, prepare_interval=None):
        self.response, self.split_stage = response, split_stage
        self.free_pair, self.prepare_interval = free_pair, prepare_interval

    def stage(self, time):
        split = self.split_stage(float(time))
        stage, r = split.current, self.response
        if stage.time != float(time):
            raise ValueError('split provider returned the wrong current mean time')
        if any(len(actions) != len(r.sectors) or not all(callable(a) for a in actions)
               for actions in (split.residual_actions, split.residual_adjoints)):
            raise ValueError('complete independently formed H-H0 and adjoint actions required')
        return split

    def context(self, time, state):
        r, stage = self.response, self.stage(time).current
        density, gradient, ell = r.moments(stage, state)
        potential = np.asarray(r.kernel(density, q_bohr_inv=r.q_bohr_inv), dtype=complex)
        if potential.shape != r.density_shape or not np.isfinite(potential).all():
            raise ValueError('complete finite complex Kq density potential required')
        velocity = state.p/r.masses
        force = -r.curvature @ state.u-gradient
        fdot, pidot = np.zeros(0, complex), np.zeros(0, complex)
        if r.gamma:
            velocity = (state.p-r.sqrt_n*r.charge @ state.coordinate)/r.masses
            force -= r.sqrt_n*r.mixed @ state.coordinate
            fdot = state.momentum/r.field_mass
            pidot = (-r.sqrt_n*(r.polarizations @ ell+r.mixed.T @ state.u
                      -r.charge.T @ (state.p/r.masses))
                      -r.ncells*(r.contact+r.charge.T @ (r.charge/r.masses[:, None])) @ state.coordinate)
        return FeedbackContext(float(time), potential, SmallState(velocity, force, fdot, pidot))

    def edge_rhs(self, context, index, pair, small):
        r, split = self.response, self.stage(context.time)
        stage, edge = split.current, r.edges[index]
        source, target = r.sectors[edge.source], r.sectors[edge.target]
        ps, pt = stage.frames[edge.source], stage.frames[edge.target]
        left, right = pair
        # Enforce rectangular dimensions without materializing any other edge.
        if left.shape != (target.shape[0], source.shape[1]) or right.shape != (target.shape[1], source.shape[0]):
            raise ValueError('original independently rectangular full-PW legs required')
        ionic = edge.coupling.prepare(small.u.reshape(-1, 3))
        potential = context.potential[:, source.density_component]
        drive = edge.coupling.bridge.local_apply(potential, ps)+ionic.apply(ps)
        adjoint = edge.coupling.bridge.local_adjoint_apply(potential, pt)+ionic.adjoint_apply(pt)
        if r.gamma:
            field = r.sqrt_n*r.polarizations.T @ small.coordinate
            vertices = r.field_vertices[edge.source]
            drive += vertices.apply(field, ps)
            adjoint += vertices.apply(field, pt, adjoint=True)
        return (-1j*(split.residual_actions[edge.target](left)+drive @ source.d),
                1j*(split.residual_adjoints[edge.source](right.conj().T).conj().T
                    +target.d @ adjoint.conj().T))


@dataclass(frozen=True)
class MidpointInterval:
    time: float
    dt: float
    small: SmallState
    midpoint_small: SmallState
    start: FeedbackContext
    midpoint: FeedbackContext


class ReplaySolution:
    """Stores O(steps*grid) feedback, never O(steps*Q) electronic legs."""
    def __init__(self, model, initial, intervals, small, time):
        self.model, self.initial = model, initial
        self.intervals, self.time = tuple(intervals), float(time)
        self.state = SourceLocalTangent(initial.count, self.edge, small)

    def edge(self, index):
        pair = self.initial.edge(index)
        for interval in self.intervals:
            pair = replay_interval(self.model, index, pair, interval)
        return pair

    @property
    def operation_counts(self):
        """Live counts include later lazy readout/replay, not only evolve()."""
        return asdict(self.model.counts)

    @property
    def context_bytes(self):
        return sum(context.potential.nbytes
                   +sum(getattr(context.small_rhs, name).nbytes for name in NAMES)
                   for interval in self.intervals
                   for context in (interval.start, interval.midpoint)) + sum(
            getattr(small, name).nbytes for interval in self.intervals
            for small in (interval.small, interval.midpoint_small) for name in NAMES)


def replay_interval(model, index, pair, interval):
    if hasattr(model, 'counts'):
        model.counts.replayed_intervals += 1
    h = interval.dt
    k1 = model.edge_rhs(interval.start, index, pair, interval.small)
    midpoint = model.free_pair(index, pair_combination((1, pair), (h/2, k1)), h/2)
    k2 = model.edge_rhs(interval.midpoint, index, midpoint, interval.midpoint_small)
    free_half = model.free_pair(index, pair, h/2)
    return model.free_pair(index, pair_combination((1, free_half), (h, k2)), h/2)


class SourceLocalMidpoint:
    """Generic model protocol: response.edges, context(t,state),
    edge_rhs(context,index,pair,small), free_pair(index,pair,tau), and optional
    prepare_interval(t,h). A context exposes small_rhs and potential.nbytes for
    payload accounting; potential may be a budgeted disk-backed array. Native
    models own complete reductions, physical factors, and action error controls.
    Historical stage lookup must be deterministic during replay, independent of
    the most recently prepared interval.
    """
    def __init__(self, model):
        self.model = _CountedModel(model)

    def evolve(self, initial, times):
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError('finite strictly increasing history-covered time grid required')
        if not isinstance(initial, SourceLocalTangent):
            initial = SourceLocalTangent.from_joint(initial)
        if initial.count != len(self.model.response.edges):
            raise ValueError('complete original edge source required')
        intervals, small = [], initial.small
        current = ReplaySolution(self.model, initial, intervals, small, times[0])
        for time, end in zip(times[:-1], times[1:]):
            dt = float(end-time)
            if self.model.prepare_interval is not None:
                self.model.prepare_interval(float(time), dt)
            start = self.model.context(time, current.state)
            midpoint_small = small_combination((1, small), (dt/2, start.small_rhs))

            def midpoint_edge(index, current=current, start=start, small=small, dt=dt):
                pair = current.state.edge(index)
                derivative = self.model.edge_rhs(start, index, pair, small)
                return self.model.free_pair(index, pair_combination((1, pair), (dt/2, derivative)), dt/2)

            midpoint_state = SourceLocalTangent(initial.count, midpoint_edge, midpoint_small)
            midpoint = self.model.context(time+dt/2, midpoint_state)
            intervals.append(MidpointInterval(float(time), dt, small, midpoint_small, start, midpoint))
            small = small_combination((1, small), (dt, midpoint.small_rhs))
            midpoint_state.release()
            current.state.release()
            current = ReplaySolution(self.model, initial, intervals, small, end)
        return current


def refinement_difference(coarse, fine, observable, norm):
    """User-chosen physical observable/norm; /3 assumes observed order two.

    observable takes a ReplaySolution, including its time and lazy state. This
    does not establish the asymptotic regime, free-action accuracy or mean error.
    """
    difference = float(norm(observable(fine)-observable(coarse)))
    return dict(difference=difference, second_order_fine_error_estimate=difference/3)
