"""Initial tensor-product ordered sources, never a re-preparation at later times.

A directed probe P denotes O_P(q)=sum_s w_s Tr(P_ts† X_ts)+p† y_q.
AB is <O_A(q) O_B(q)†>_c, BA is <O_B(q)† O_A(q)>_c. Thus the
right coefficients enter sources WITHOUT conjugation. The dagger supplies
the explicitly opposite-q observable, not another photon or a second edge.
The source carries N*w_s and the readout w_s, exactly as the admitted raw
quasifree collective preparation. Local ionic W is already canonical Fourier
W, not W/N. The global field W has no matter-volume rescaling.
"""
from dataclasses import dataclass

import numpy as np

from .joint_bloch_response import JointTangent
from .prepared_reference import array, authority


@dataclass(frozen=True)
class JointProbe:
    actions: tuple
    adjoints: tuple
    u: np.ndarray
    p: np.ndarray
    coordinate: np.ndarray
    momentum: np.ndarray
    identity: str
    units: str

    def validate(self, response):
        authority(self.identity, 'directed observable identity')
        authority(self.units, 'observable units')
        if len(self.actions) != len(response.edges) or len(self.adjoints) != len(response.edges) or not all(callable(a) for a in self.actions+self.adjoints):
            raise ValueError('full-PW forward/independent adjoint per edge required')
        for name, size in (('u', response.ndof), ('p', response.ndof),
                           ('coordinate', response.nfield), ('momentum', response.nfield)):
            array(getattr(self, name), (size,), 'probe '+name)


@dataclass(frozen=True)
class JointOrderedSource:
    electronic: JointTangent
    ionic: JointTangent
    field: JointTangent
    ordering: str
    right_identity: str
    preparation_identity: str

    @property
    def parts(self):
        return dict(electronic=self.electronic, ionic=self.ionic, field=self.field)


def linear_combination(*terms):
    """Combine (complex scalar,tangent) pairs without tying the two legs."""
    if not terms:
        raise ValueError('at least one tangent required')
    first = terms[0][1]
    ket = tuple(sum(a*x.ket[i] for a, x in terms) for i in range(len(first.ket)))
    bra = tuple(sum(a*x.bra[i] for a, x in terms) for i in range(len(first.bra)))
    values = [sum(a*getattr(x, name) for a, x in terms)
              for name in ('u', 'p', 'coordinate', 'momentum')]
    return JointTangent(ket, bra, *values)


def _ordered(value, size, name):
    w = array(value, (2*size,)*2, name)
    if not np.allclose(w, w.conj().T, atol=1e-11, rtol=1e-11):
        raise ValueError(name+' must be Hermitian ordered covariance')
    if size and np.linalg.eigvalsh(w).min(initial=0) < -1e-10*max(1., np.linalg.norm(w, 2)):
        raise ValueError(name+' must be positive semidefinite; no clipping')
    return w


def _canonical(size):
    return np.block([[np.zeros((size, size)), np.eye(size)],
                     [-np.eye(size), np.zeros((size, size))]])


class JointProductPreparation:
    """Original quasifree D x saved ionic W x saved global field W.

    q and -q ionic matrices are both explicit: BA uses W_-q.T, never W_q†.
    Cross origins are exactly zero ONLY at this original t0 product boundary.
    A caller must supply authenticated original frames/D identities; no new
    quasifree state is fitted, D inverted, or empty PW complement discarded.
    """
    def __init__(self, response, initial_stage, *, ionic_ordered, ionic_minus_ordered,
                 field_ordered, electronic_identity, ionic_identity,
                 field_identity, product_identity):
        response.validate_stage(initial_stage)
        if initial_stage.time != 0:
            raise ValueError('ordered sources belong to the original t0 preparation only')
        for value in (electronic_identity, ionic_identity, field_identity, product_identity):
            authority(value, 'original tensor-product preparation identity')
        self.response, self.initial_stage = response, initial_stage
        self.electronic_identity, self.ionic_identity = electronic_identity, ionic_identity
        self.field_identity, self.identity = field_identity, product_identity
        self.ionic = _ordered(ionic_ordered, response.ndof, 'stored ionic Wq')
        self.ionic_minus = _ordered(ionic_minus_ordered, response.ndof, 'stored ionic Wminus')
        if not np.allclose(self.ionic-self.ionic_minus.T, 1j*_canonical(response.ndof), atol=1e-10, rtol=1e-10):
            raise ValueError('original q/-q ionic canonical commutator required')
        if response.gamma and not np.allclose(self.ionic, self.ionic_minus, atol=1e-12, rtol=1e-12):
            raise ValueError('Gamma is its own opposite q')
        self.field = _ordered(field_ordered, response.nfield, 'stored global field W')
        if not np.allclose(self.field-self.field.T, 1j*_canonical(response.nfield), atol=1e-10, rtol=1e-10):
            raise ValueError('original global Maxwell canonical commutator required')

    def source(self, right, *, ordering='AB'):
        if ordering not in ('AB', 'BA'):
            raise ValueError('AB or BA required')
        r, stage = self.response, self.initial_stage
        right.validate(r)
        ket, bra = [], []
        for e, action, adjoint in zip(r.edges, right.actions, right.adjoints):
            s, t = r.sectors[e.source], r.sectors[e.target]
            ps, pt = stage.frames[e.source], stage.frames[e.target]
            bp = array(action(ps), (pt.shape[0], ps.shape[1]), 'B Psi_s')
            ba = array(adjoint(pt), (ps.shape[0], pt.shape[1]), 'B adjoint Psi_t')
            overlap = pt.conj().T @ bp
            scale = r.ncells*s.weight
            if ordering == 'AB':
                ket.append(scale*bp @ s.d)
                bra.append(-scale*t.d @ overlap @ s.d @ ps.conj().T)
            else:
                ket.append(-scale*pt @ (t.d @ overlap @ s.d))
                bra.append(scale*t.d @ ba.conj().T)
        zero = r.zero_tangent()
        electronic = r.tangent(ket, bra, zero.u, zero.p, zero.coordinate, zero.momentum)
        wi = self.ionic if ordering == 'AB' else self.ionic_minus.T
        wf = self.field if ordering == 'AB' else self.field.T
        ions = wi @ np.concatenate((right.u, right.p))
        fields = wf @ np.concatenate((right.coordinate, right.momentum))
        ionic = r.tangent(zero.ket, zero.bra, ions[:r.ndof], ions[r.ndof:], zero.coordinate, zero.momentum)
        field = r.tangent(zero.ket, zero.bra, zero.u, zero.p, fields[:r.nfield], fields[r.nfield:])
        return JointOrderedSource(electronic, ionic, field, ordering, right.identity, self.identity)

    def coherent_impulse(self, right):
        """-i[O_B†,rho] and canonical kicks, NOT an ordered noise source."""
        ab, ba = self.source(right, ordering='AB'), self.source(right, ordering='BA')
        return linear_combination(*[term for name in ab.parts for term in
            ((-1j, ab.parts[name]), (1j, ba.parts[name]))])
