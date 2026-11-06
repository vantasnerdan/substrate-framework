"""Selected joint observables and TWO-sided original-preparation pullbacks.

A pullback(time,probe) returns the physical dual probe at t0 such that
readout_t(probe,T_q(t,0)x)=readout_0(pullback(time,probe),x) for EVERY initial
physical x. It must include all electron/ion/field mixing, the moving-frame
connection, and the full empty PW complement. No identity fallback exists.
The companion finite oracle constructs and exercises an independent DOP853
fundamental solution, its adjoint pullback, and direct two-sided covariance.
"""
import numpy as np

from .joint_ordered_sources import JointProbe
from .prepared_reference import authority, array


def readout(response, stage, tangent, probe):
    """Sesquilinear directed q readout; probe dagger is the opposite edge."""
    response.validate_stage(stage)
    probe.validate(response)
    electronic = 0j
    for e, l, b, a, ad in zip(response.edges, tangent.ket, tangent.bra, probe.actions, probe.adjoints):
        ps, pt = stage.frames[e.source], stage.frames[e.target]
        ap = array(a(ps), l.shape, 'A Psi_source')
        ah = array(ad(pt), b.conj().T.shape, 'A dagger Psi_target')
        electronic += response.sectors[e.source].weight*(np.vdot(ap, l)+np.einsum('ij,ji->', b, ah))
    ionic = np.vdot(probe.u, tangent.u)+np.vdot(probe.p, tangent.p)
    field = np.vdot(probe.coordinate, tangent.coordinate)+np.vdot(probe.momentum, tangent.momentum)
    return dict(electronic=complex(electronic), ionic=complex(ionic), field=complex(field),
                total=complex(electronic+ionic+field))


def source_readout(response, stage, source, probe):
    """Keys are INITIAL origins; each can contribute at every final location."""
    parts = {name: readout(response, stage, tangent, probe)['total']
             for name, tangent in source.parts.items()}
    return dict(parts, total=sum(parts.values()), ordering=source.ordering)


class TwoSidedJointCovariance:
    """C_AB(t,s) and C_BA(t,s), including arbitrary nonzero/equal second times.

    dynamics must distinguish 'screened-joint' (the coupled tangent) from
    'bare-electronic' (externally driven electrons, no feedback). A bare result
    is never labeled joint screened noise. The preparation is fixed at t0;
    pulling BOTH probes back retains generated cross correlations without
    storing the full covariance or recreating product states at s.
    """
    def __init__(self, preparation, pullback, *, dynamics, propagation_identity):
        if not callable(pullback) or dynamics not in ('screened-joint', 'bare-electronic'):
            raise ValueError('actual two-sided dual propagator and explicit dynamics required')
        authority(propagation_identity, 'propagator and moving baseline identity')
        self.preparation, self.pullback = preparation, pullback
        self.dynamics, self.propagation_identity = dynamics, propagation_identity

    def pair(self, left, right, time, second_time):
        if not np.isfinite([time, second_time]).all() or min(time, second_time) < 0:
            raise ValueError('explicit finite nonnegative times from original preparation required')
        a, b = self.pullback(float(time), left), self.pullback(float(second_time), right)
        r, stage = self.preparation.response, self.preparation.initial_stage
        a.validate(r)
        b.validate(r)
        ab = source_readout(r, stage, self.preparation.source(b, ordering='AB'), a)
        ba = source_readout(r, stage, self.preparation.source(b, ordering='BA'), a)
        return dict(AB=ab, BA=ba, commutator=ab['total']-ba['total'],
            symmetrized=(ab['total']+ba['total'])/2,
            time=float(time), second_time=float(second_time), dynamics=self.dynamics,
            preparation_identity=self.preparation.identity,
            propagation_identity=self.propagation_identity,
            convention='<O_A(q,t) O_B(q,s)^dagger> connected; opposite q explicit via dagger')


def total_force_probe(response, direction, *, identity):
    """direction†[-Gq-Cphysical(q)U-sqrtN*T*field] in original atom order."""
    v = array(direction, (response.ndof,), 'force direction')
    operators = tuple(e.coupling.prepare(v.reshape(-1, 3)) for e in response.edges)
    actions = tuple((lambda x, op=op: -op.apply(x)) for op in operators)
    adjoints = tuple((lambda x, op=op: -op.adjoint_apply(x)) for op in operators)
    field = -response.sqrt_n*response.mixed.conj().T @ v if response.gamma else np.zeros(0)
    return JointProbe(actions, adjoints, -response.curvature.conj().T @ v,
        np.zeros(response.ndof), field, np.zeros(response.nfield), identity, 'Ry/bohr')


def local_momentum_kick(response, cell, direction, amplitude):
    """One q component of a localized canonical coherent kick, not covariance."""
    position = array(cell, (3,), 'cell lattice index', real=True)
    if not np.array_equal(position, np.round(position)):
        raise ValueError('integer original cell lattice indices required')
    v = array(direction, (response.ndof,), 'canonical kick direction')
    if not np.isfinite(amplitude):
        raise ValueError('finite coherent kick amplitude required')
    zero = response.zero_tangent()
    p = amplitude*np.exp(-2j*np.pi*response.q_fractional @ position)*v/response.sqrt_n
    return response.tangent(zero.ket, zero.bra, zero.u, p, zero.coordinate, zero.momentum)


def spatial_pair(q_fractional, ordered_pairs, cell, other, *, ncells):
    """Full q/-q orthonormal reconstruction: sum_q exp(iq(R-R')) Cq/N.

    ordered_pairs already pairs q with -q; no conjugation/dagger is added.
    All N distinct q must be supplied; selected q is not a local reconstruction.
    """
    q = array(q_fractional, (ncells, 3), 'complete q grid', real=True)
    if len(np.unique(np.mod(q, 1), axis=0)) != ncells:
        raise ValueError('complete distinct periodic q grid required')
    values = np.asarray(ordered_pairs)
    if values.shape[0] != ncells:
        raise ValueError('one ordered pair per q required')
    delta = array(cell, (3,), 'cell', real=True)-array(other, (3,), 'other cell', real=True)
    return np.tensordot(np.exp(2j*np.pi*q @ delta)/ncells, values, axes=(0, 0))
