"""Directed finite-q response of the declared prepared finite scalar.

X_ts=L_ts Psi_s†+Psi_t R_ts; L(nt,rs) and R(rt,ns) are independent.
With hbar=1 and H in Ry, time is inverse Ry. No equilibrium, gap division,
PW-square density, band projection, fitted damping or Hermitization is used.
The supplied Cq is unmassweighted Ry/bohr², Kq maps bohr^-3 to Ry, and the
canonical ionic equations are udot=p/M, pdot=-Cq u-Gq. The mean and both
endpoint Hamiltonians always come from ONE current MeanEvaluation.
"""
from dataclasses import dataclass
from itertools import product

import numpy as np

from .prepared_factor_dynamics import PreparedFactorDynamics
from .prepared_reference import array, authority, frozen
from .qe_bloch_edge import BlochEdge


@dataclass(frozen=True)
class ResponseEdge:
    """Indices into the mean sectors, and an existing BlochFiniteQIonic API.

    Finite mathematical adapters may implement that same action/moment API;
    their edge must still be a checked BlochEdge, not guessed endpoint labels.
    """
    source: int
    target: int
    coupling: object


@dataclass(frozen=True)
class BlochTangent:
    ket: tuple[np.ndarray, ...]
    bra: tuple[np.ndarray, ...]
    u: np.ndarray
    p: np.ndarray


class PreparedBlochResponse:
    """Complete same-spin q permutation, coupled to the existing moving mean.

    kernel is callable(density, q_bohr_inv=q), e.g. MaterialReference.apply_kernel.
    curvature is the supplied Fourier coefficient Cq of the SAME fixed scalar,
    not a new-state electronic Hessian. The caller declares its authority and
    must supply its exact opposite-q conjugate in a separate -q calculation.
    Physical real perturbations may relate opposite-q solutions afterwards;
    no such constraint is imposed on complex ordered sources here. At q=0 or
    self-inverse q each source sector is still visited exactly once.
    """
    def __init__(self, mean_model, edges, *, kernel, curvature, authority_text):
        if not isinstance(mean_model, PreparedFactorDynamics):
            raise TypeError('existing PreparedFactorDynamics required')
        authority(authority_text, 'common finite-q scalar authority')
        if not callable(kernel):
            raise TypeError('physical Kq action required')
        self.mean_model = mean_model
        self.edges = tuple(edges)
        self.kernel = kernel
        self.authority = authority_text
        r, sectors = mean_model.reference, mean_model.sectors
        if r.ndof % 3:
            raise ValueError('finite-q ionic API requires complete atom-order xyz coordinates')
        count = len(sectors)
        if (len(self.edges) != count
                or sorted(e.source for e in self.edges) != list(range(count))
                or sorted(e.target for e in self.edges) != list(range(count))):
            raise ValueError('one outgoing AND incoming edge per mean sector required')
        first = self.edges[0].coupling.edge
        if not isinstance(first, BlochEdge):
            raise TypeError('checked BlochEdge required')
        self.q_fractional = frozen(first.q_fractional)
        self.q_bohr_inv = frozen(first.q_cart_bohr_inv)
        mesh = first.mesh
        if int(np.prod(mesh)) != r.ncells:
            raise ValueError('uniform complete mesh and mean N differ')
        slots = set()
        for link in self.edges:
            c, edge = link.coupling, link.coupling.edge
            if not isinstance(edge, BlochEdge):
                raise TypeError('checked BlochEdge required')
            edge.validate(c.inputs.reference)
            s, t = sectors[link.source], sectors[link.target]
            if (edge.mesh != mesh or not np.array_equal(edge.q_fractional, self.q_fractional)
                    or not np.array_equal(edge.q_cart_bohr_inv, self.q_bohr_inv)
                    or s.basis_identity != edge.source_basis.identity
                    or t.basis_identity != edge.target_basis.identity
                    or s.initial_frame.shape[0] != edge.source_dimension
                    or t.initial_frame.shape[0] != edge.target_dimension
                    or s.density_component != edge.ispin-1
                    or t.density_component != edge.ispin-1
                    or s.weight != t.weight
                    or s.weight != edge.source_basis.slot.weight
                    or t.weight != edge.target_basis.slot.weight
                    or c.ndof != r.ndof or c.bridge.source_fft.nnr != r.nbar.shape[0]
                    or c.bridge.grid_measure_bohr3 != mean_model.measure):
                raise ValueError('q/endpoint/spin/weight/coordinate/grid authority mismatch')
            slots.add((edge.ispin, tuple(edge.source_basis.slot.mesh_index)))
        expected = {(spin+1, index) for spin in range(r.nbar.shape[1])
                    for index in product(*(range(n) for n in mesh))}
        if slots != expected or len(slots) != count:
            raise ValueError('every explicit spin/full-grid source slot required exactly once')
        c = array(curvature, (r.ndof, r.ndof), 'unmassweighted Cq')
        if not np.allclose(c, c.conj().T, rtol=1e-12, atol=1e-12):
            raise ValueError('Cq must be Hermitian; no symmetrization is performed')
        if np.all(self.q_fractional == 0) and not np.allclose(c, r.curvature, rtol=1e-12, atol=1e-12):
            raise ValueError('Gamma curvature differs from the mean finite scalar')
        self.curvature = frozen(c)

    def tangent(self, ket, bra, u, p):
        """Validate unequal endpoint PW/rank shapes; never impose R=L†."""
        if len(ket) != len(self.edges) or len(bra) != len(self.edges):
            raise ValueError('one independent ket/bra pair per directed edge required')
        sectors = self.mean_model.sectors
        left, right = [], []
        for e, l, b in zip(self.edges, ket, bra):
            ns, rs = sectors[e.source].initial_frame.shape
            nt, rt = sectors[e.target].initial_frame.shape
            left.append(array(l, (nt, rs), 'directed ket'))
            right.append(array(b, (rt, ns), 'directed bra'))
        ndof = self.mean_model.reference.ndof
        return BlochTangent(tuple(left), tuple(right), array(u, (ndof,), 'q displacement'),
                            array(p, (ndof,), 'q momentum'))

    def zero_tangent(self):
        sectors = self.mean_model.sectors
        return self.tangent(
            tuple(np.zeros((sectors[e.target].initial_frame.shape[0],
                            sectors[e.source].initial_frame.shape[1]), complex) for e in self.edges),
            tuple(np.zeros((sectors[e.target].initial_frame.shape[1],
                            sectors[e.source].initial_frame.shape[0]), complex) for e in self.edges),
            np.zeros(self.mean_model.reference.ndof), np.zeros(self.mean_model.reference.ndof))

    def _evaluation(self, evaluation):
        if evaluation.model is not self.mean_model:
            raise ValueError('evaluation must belong to the same moving mean model')

    def moments(self, evaluation, tangent):
        """Return n_q and G_q with source-sector weight exactly once.

        G_q[a]=sum_s w_s Tr(g_ts[e_a]† X_ts), so vdot(u,G_q) is the
        directed pairing. There is no D here: the independent legs already
        represent delta rho, NOT delta Psi. No conjugate/opposite edge added.
        """
        self._evaluation(evaluation)
        tangent = self.tangent(tangent.ket, tangent.bra, tangent.u, tangent.p)
        density = np.zeros(self.mean_model.reference.nbar.shape, complex)
        gradient = np.zeros(self.mean_model.reference.ndof, complex)
        for e, l, b in zip(self.edges, tangent.ket, tangent.bra):
            ps, pt = evaluation.state.frames[e.source], evaluation.state.frames[e.target]
            w = self.mean_model.sectors[e.source].weight
            n1, g1 = e.coupling.moments(l, ps, weight=w)
            n2, g2 = e.coupling.moments(pt, b.conj().T, weight=w)
            density += n1 + n2
            gradient += g1 + g2
        return density, gradient

    def rhs(self, evaluation, tangent):
        """Return leg derivatives including the moving-frame connection.

        Reconstruct Xdot with Ldot Psi_s†+L Psidot_s†+Psidot_t R+Psi_t Rdot.
        It equals -i(H_t X-X H_s+deltaH_ts rho_s-rho_t deltaH_ts).
        D_s occurs on the ket drive and D_t on the bra drive; they are never
        interchanged, inverted, replaced by MV labels, or multiplied into X twice.
        The independent adjoint action is used only to evaluate Psi_t† deltaH,
        not to assume the complex deltaH or source is Hermitian.
        """
        density, gradient = self.moments(evaluation, tangent)
        field = array(self.kernel(density, q_bohr_inv=self.q_bohr_inv),
                      density.shape, 'Kq density')
        ket, bra = [], []
        model = self.mean_model
        for e, l, b in zip(self.edges, tangent.ket, tangent.bra):
            s, t = model.sectors[e.source], model.sectors[e.target]
            ps, pt = evaluation.state.frames[e.source], evaluation.state.frames[e.target]
            bridge = e.coupling.bridge
            ionic = e.coupling.prepare(tangent.u.reshape(-1, 3))
            v = field[:, s.density_component]
            drive = bridge.local_apply(v, ps) + ionic.apply(ps)
            drive_adjoint = bridge.local_adjoint_apply(v, pt) + ionic.adjoint_apply(pt)
            ket.append(-1j*(model._h_action(e.target, l, evaluation) + drive @ s.d))
            bra.append(1j*(model._h_action(e.source, b.conj().T, evaluation).conj().T
                            + t.d @ drive_adjoint.conj().T))
        return BlochTangent(tuple(ket), tuple(bra), tangent.p/model.reference.masses,
                            -self.curvature @ tangent.u-gradient)

    def vertex_source(self, evaluation, actions, adjoint_actions, *, ordering):
        """Raw quasifree directed source for a supplied B_ts, all full PW rows.

        AB: (I-rho_t) B_ts rho_s. BA: rho_t B_ts (I-rho_s).
        AB-BA is [B,rho]_ts; -i times it is the impulse response source.
        These are UNWEIGHTED single-sector electronic sources. The readout
        inserts w_s; a chosen collective/probe normalization is caller-owned.
        I is the full PW identity, including the initially empty complement.
        Ionic/cross preparation is not inferred; returned ionic legs are zero.
        """
        self._evaluation(evaluation)
        if ordering not in ('AB', 'BA') or len(actions) != len(self.edges) or len(adjoint_actions) != len(self.edges):
            raise ValueError('AB/BA and one forward/adjoint callable per edge required')
        ket, bra = [], []
        for e, action, adjoint in zip(self.edges, actions, adjoint_actions):
            s, t = self.mean_model.sectors[e.source], self.mean_model.sectors[e.target]
            ps, pt = evaluation.state.frames[e.source], evaluation.state.frames[e.target]
            bp = array(action(ps), (pt.shape[0], ps.shape[1]), 'B Psi_s')
            badj = array(adjoint(pt), (ps.shape[0], pt.shape[1]), 'B† Psi_t')
            overlap = pt.conj().T @ bp
            if ordering == 'AB':
                ket.append(bp @ s.d)
                bra.append(-t.d @ overlap @ s.d @ ps.conj().T)
            else:
                ket.append(-pt @ (t.d @ overlap @ s.d))
                bra.append(t.d @ badj.conj().T)
        ndof = self.mean_model.reference.ndof
        return self.tangent(ket, bra, np.zeros(ndof), np.zeros(ndof))

    def readout(self, evaluation, tangent, actions, adjoint_actions, *, ionic_u=None, ionic_p=None):
        """sum_s w_s Tr(A_ts† X_ts) + vdot(a_u,u)+vdot(a_p,p).

        This is a declared directed SESQUILINEAR observable pairing, not an
        implicit Gamma bilinear convention. For ordered AB use the explicitly
        chosen opposite-edge observable as A†. Return raw complex terms.
        """
        self._evaluation(evaluation)
        tangent = self.tangent(tangent.ket, tangent.bra, tangent.u, tangent.p)
        if len(actions) != len(self.edges) or len(adjoint_actions) != len(self.edges):
            raise ValueError('one action and independent adjoint per edge required')
        electronic = 0j
        for e, l, b, action, adjoint in zip(self.edges, tangent.ket, tangent.bra, actions, adjoint_actions):
            ps, pt = evaluation.state.frames[e.source], evaluation.state.frames[e.target]
            ap = array(action(ps), l.shape, 'A Psi_s')
            ad = array(adjoint(pt), b.conj().T.shape, 'A† Psi_t')
            electronic += self.mean_model.sectors[e.source].weight*(np.vdot(ap, l)+np.einsum('ij,ji->', b, ad))
        ionic = 0j
        for coefficient, value in ((ionic_u, tangent.u), (ionic_p, tangent.p)):
            if coefficient is not None:
                ionic += np.vdot(array(coefficient, value.shape, 'ionic readout'), value)
        return dict(electronic=complex(electronic), ionic=complex(ionic), total=complex(electronic+ionic))
