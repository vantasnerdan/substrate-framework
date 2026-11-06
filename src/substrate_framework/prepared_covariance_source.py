"""Selected ordered covariance columns and physical action-backed readout.

Uses the existing FactorPreparation quasifree marginals and IonicPreparation
canonical covariance. Explicit product or authorized raw-factor cross moments
are accepted; unknown cross never means zero. Joint states with non-quasifree
or correlated electronic-sector marginals are NOT representable by this API.
No finite-time propagation or arbitrary two-time pullback is supplied here.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .physical_bilinear_preparation import IonicPreparation
from .prepared_reference import array, authority, frozen
from .prepared_factor_dynamics import FactorTangent


@dataclass(frozen=True)
class PreparedProbe:
    """O=constant + sum_s w_s Tr(A_s rho_s) + u_coeff.u + p_coeff.p.

Electronic actions are UNWEIGHTED; collective fluctuations replace rho by
X=sqrt(N)*delta rho. Explicit adjoints act on arbitrary full-PW vectors, not
only saved columns. Constants affect means, never connected covariance.
"""
    actions: tuple
    adjoints: tuple
    u: np.ndarray
    p: np.ndarray
    constant: complex
    identity: str
    units: str

    def adjoint(self):
        return PreparedProbe(self.adjoints, self.actions, np.asarray(self.u).conj(),
                             np.asarray(self.p).conj(), complex(self.constant).conjugate(),
                             self.identity + '^dagger', self.units)

    def validate(self, model):
        authority(self.identity, 'probe identity')
        authority(self.units, 'probe units')
        if (len(self.actions) != len(model.sectors) or len(self.adjoints) != len(model.sectors)
                or not all(callable(x) for x in self.actions + self.adjoints)
                or not np.isfinite(self.constant)):
            raise ValueError('one genuine action/adjoint per explicit sector required')
        array(self.u, (model.reference.ndof,), 'probe u')
        array(self.p, (model.reference.ndof,), 'probe p')


@dataclass(frozen=True)
class FactorCrossPreparation:
    """<delta rho_X (u,p)_j>=Psi E_j Psi^dagger, per sector.

This is the raw nonorthogonal equivalent of CrossPreparation.matrices, with
its SAME collective scaling and explicit joint-state authority convention.
The supplied joint construction must preserve the admitted independent
quasifree electronic marginals and ionic covariance. A PSD block alone is not
such authority. Cross support outside this factor span is explicitly rejected
by this representation, not dropped. No additional electron-electron
cumulants or cross-sector correlations can be encoded here.
"""
    coefficients: tuple[np.ndarray, ...]
    state_authority: str
    representation_authority: str

    @classmethod
    def product(cls, model, state_authority):
        authority(state_authority, 'explicit electron/ion product preparation')
        return cls(tuple(np.zeros((2 * model.reference.ndof, len(s.d), len(s.d)), complex)
                         for s in model.sectors), state_authority,
                   'Explicit product; all electron-ion cross moments zero')

    def validate(self, model):
        authority(self.state_authority, 'actual joint-state construction')
        authority(self.representation_authority, 'joint representation and marginal authority')
        if len(self.coefficients) != len(model.sectors):
            raise ValueError('one raw-factor cross coefficient array per sector required')
        for s, e in zip(model.sectors, self.coefficients):
            e = array(e, (2 * model.reference.ndof, len(s.d), len(s.d)), 'cross coefficients')
            if not np.allclose(e, e.conj().swapaxes(-1, -2), atol=1e-12, rtol=1e-12):
                raise ValueError('commuting electron-ion Hermitian-coordinate cross moments required')


@dataclass(frozen=True)
class OrderedSource:
    """Parts denote INITIAL covariance origin, not final electronic/ionic location."""
    electronic: FactorTangent
    ionic: FactorTangent | None
    cross: FactorTangent | None
    ordering: str
    right_identity: str
    preparation_authority: str

    @property
    def complete(self):
        return self.ionic is not None and self.cross is not None


class PreparedCovariance:
    """Authorized initial covariance-source construction, never a new state fit."""
    def __init__(self, model, *, electronic_authority, ionic=None, cross=None):
        authority(electronic_authority, 'independent quasifree sector/number/spin authority')
        if ionic is not None and (not isinstance(ionic, IonicPreparation)
                                  or ionic.ordered.shape != (2 * model.reference.ndof,) * 2):
            raise ValueError('existing full-coordinate IonicPreparation required')
        if cross is not None:
            if ionic is None or not isinstance(cross, FactorCrossPreparation):
                raise ValueError('ionic authority and representable FactorCrossPreparation required')
            cross.validate(model)
        self.model, self.ionic, self.cross = model, ionic, cross
        self.electronic_authority = electronic_authority

    def source(self, right, *, ordering='AB'):
        """At t0: source=Cov(coordinate,B) for AB, Cov(B,coordinate) for BA.

AB electronic matrix: N*w_B*(I-rho) B rho.
BA electronic matrix: N*w_B*rho B (I-rho).
The full identity is applied implicitly to FULL-ROW B Psi. No projector onto
saved bands replaces it. Products/ionic transposes do not conjugate B.
"""
        if ordering not in ('AB', 'BA'):
            raise ValueError('ordering must be AB or BA')
        model = self.model
        right.validate(model)
        n = model.reference.ndof
        ket, bra = [], []
        for s, action, adjoint in zip(model.sectors, right.actions, right.adjoints):
            psi, d = s.initial_frame, s.d
            scale = model.reference.ncells * s.weight
            if ordering == 'AB':
                bpsi = array(action(psi), psi.shape, 'B Psi')
                ket.append(scale * (bpsi - psi @ (d @ (psi.conj().T @ bpsi))) @ d)
                bra.append(np.zeros_like(psi.T))
            else:
                bdagpsi = array(adjoint(psi), psi.shape, 'B adjoint Psi')
                ket.append(np.zeros_like(psi))
                bra.append(scale * d @ (bdagpsi - psi @ (d @ (psi.conj().T @ bdagpsi))).conj().T)
        electronic = model.tangent(tuple(ket), tuple(bra), np.zeros(n), np.zeros(n))
        ri = np.concatenate((right.u, right.p))
        ionic = cross_source = None
        if self.ionic is not None:
            sigma = self.ionic.ordered if ordering == 'AB' else self.ionic.ordered.T
            y = sigma @ ri
            zero = model.zero_tangent()
            ionic = FactorTangent(zero.ket, zero.bra, y[:n], y[n:])
        if self.cross is not None:
            ket, bra = [], []
            y = np.zeros(2 * n, complex)
            for s, e, action in zip(model.sectors, self.cross.coefficients, right.actions):
                psi = s.initial_frame
                ket.append(psi @ np.einsum('j,jab->ab', ri, e))
                bra.append(np.zeros_like(psi.T))
                overlap = psi.conj().T @ action(psi)
                y += s.weight * np.einsum('jab,ba->j', e, overlap)
            cross_source = model.tangent(tuple(ket), tuple(bra), y[:n], y[n:])
        return OrderedSource(electronic, ionic, cross_source, ordering, right.identity,
                             self.electronic_authority)


def readout(model, state, tangent, probe):
    """Physical bilinear trace dual, without a dense pulled probe."""
    probe.validate(model)
    value = np.asarray(probe.u) @ tangent.u + np.asarray(probe.p) @ tangent.p
    for s, psi, left, right, action in zip(model.sectors, state.frames, tangent.ket,
                                          tangent.bra, probe.actions):
        value += s.weight * (np.trace(psi.conj().T @ action(left))
                             + np.trace(right @ action(psi)))
    return complex(value)


def probe_mean(model, state, probe):
    """Uncentered mean, separately normalized from collective covariance."""
    probe.validate(model)
    value = probe.constant + np.asarray(probe.u) @ state.u + np.asarray(probe.p) @ state.p
    for s, psi, action in zip(model.sectors, state.frames, probe.actions):
        value += s.weight * np.trace(s.d @ (psi.conj().T @ action(psi)))
    return complex(value)


def ordered_readout(model, state, source, probe):
    parts = {name: None if getattr(source, name) is None else
             readout(model, state, getattr(source, name), probe)
             for name in ('electronic', 'ionic', 'cross')}
    parts.update(total=sum(parts.values()) if source.complete else None,
                 complete=source.complete, ordering=source.ordering,
                 scope='TOTAL supplied preparation' if source.complete else 'PARTIAL contributions')
    return parts


def total_force_probe(model, direction, *, identity):
    """All-coordinate total mean force and its fluctuation, not just -g Psi."""
    r = model.reference
    direction = frozen(array(direction, (r.ndof,), 'force direction'))
    actions = tuple((lambda vectors, c=s.coupling: -c.ionic_action(direction, vectors))
                    for s in model.sectors)
    adjoints = tuple((lambda vectors, c=s.coupling: -c.ionic_action(direction.conj(), vectors))
                     for s in model.sectors)
    constant = direction @ (r.ionic_force + r.support_force - r.held_gradient)
    return PreparedProbe(actions, adjoints, -r.curvature.T @ direction, np.zeros(r.ndof),
                         complex(constant), identity, 'Ry/bohr')
