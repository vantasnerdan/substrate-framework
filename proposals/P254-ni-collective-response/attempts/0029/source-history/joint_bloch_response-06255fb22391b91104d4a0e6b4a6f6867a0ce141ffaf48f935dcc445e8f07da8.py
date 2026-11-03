"""Joint finite-q tangent about explicitly supplied CURRENT physical endpoints.

No mean runtime, PreparedSector, PW-square density or copy of original frames is
created here. X_ts=L_ts Psi_s†+Psi_t R_ts has independent full-PW legs. Matter
coordinates are orthonormal Fourier U,P (sqrt(N) times homogeneous variations
at Gamma); the Maxwell coordinate is global. C_common(q)+DeltaC is used once
at every q. The uniform mean field MUST already be in both endpoint H actions.
See A28/joint-response-registration.json for scalar/data identities and oracle.
"""
from dataclasses import dataclass

import numpy as np

from .prepared_reference import array, authority


@dataclass(frozen=True)
class JointSector:
    """Copy-free endpoint specification; D is the original raw-frame matrix."""
    shape: tuple
    d: np.ndarray
    weight: float
    density_component: int
    basis_identity: str

    def __post_init__(self):
        if len(self.shape) != 2 or min(self.shape) < 1 or not np.isfinite(self.weight) or self.weight <= 0:
            raise ValueError('positive full-PW/rank dimensions and original weight required')
        d = array(self.d, (self.shape[1],)*2, 'original complex D')
        if not np.allclose(d, d.conj().T, atol=1e-12, rtol=1e-12):
            raise ValueError('Hermitian original D required; no repair')
        authority(self.basis_identity, 'endpoint basis identity')


@dataclass(frozen=True)
class JointStage:
    """H callbacks act on arbitrary full-row blocks, not just saved columns.

    adjoints must independently implement H†. Physical H is Hermitian, but
    the explicit callback keeps right multiplication correct without assuming
    a complex tangent drive is Hermitian. Frames are references, not copies.
    """
    frames: tuple
    h_actions: tuple
    h_adjoints: tuple
    time: float
    identity: str


@dataclass(frozen=True)
class JointTangent:
    ket: tuple
    bra: tuple
    u: np.ndarray
    p: np.ndarray
    coordinate: np.ndarray
    momentum: np.ndarray


@dataclass(frozen=True)
class CartesianFieldActions:
    """Complex-linear extension of three real Cartesian physical L vertices."""
    axes: tuple

    @classmethod
    def from_vertex(cls, vertex):
        return cls(tuple(vertex.prepare(axis) for axis in np.eye(3)))

    def apply(self, direction, vectors, *, adjoint=False):
        direction = array(direction, (3,), 'complex tangent field direction')
        result = np.zeros_like(vectors, dtype=complex)
        for coefficient, action in zip(direction, self.axes):
            result += (coefficient.conjugate()*action.adjoint_apply(vectors) if adjoint
                       else coefficient*action.apply(vectors))
        return result


class JointBlochResponse:
    """Mathematical q block, independent of any resident or streamed mean.

    edges use ResponseEdge(source,target,coupling). Physical rectangular
    coupling supplies prepare(u).apply/adjoint_apply, moments(l,r,weight=),
    and bridge.local_apply/local_adjoint_apply. Checked BlochEdge metadata,
    when supplied by physical couplings, is matched to the explicit endpoints.
    Finite-system adapters implement the same algebra and declare their own
    source authority. Complete outgoing/incoming permutations are mandatory.
    """
    def __init__(self, sectors, edges, *, ncells, q_fractional, q_bohr_inv,
                 density_shape, kernel, common_curvature, coefficients, masses,
                 charges, polarizations, field_mass, field_vertices=(),
                 scalar_identity, coverage_identity):
        authority(scalar_identity, 'joint scalar identity')
        authority(coverage_identity, 'complete explicit slot/q coverage')
        self.scalar_identity, self.coverage_identity = scalar_identity, coverage_identity
        self.sectors, self.edges = tuple(sectors), tuple(edges)
        self.ncells, self.kernel = ncells, kernel
        self.q_fractional = array(q_fractional, (3,), 'q fractional', real=True)
        self.q_bohr_inv = array(q_bohr_inv, (3,), 'q inverse bohr', real=True)
        self.gamma = bool(np.all(self.q_fractional == 0))
        self.density_shape = tuple(density_shape)
        self.masses = np.asarray(masses, dtype=float)
        self.ndof = len(self.masses)
        self.coefficients = coefficients
        self.curvature = array(common_curvature, (self.ndof,)*2, 'common Cq') + array(
            coefficients.delta_curvature, (self.ndof,)*2, 'fixed DeltaC')
        if not np.allclose(self.curvature, self.curvature.conj().T, atol=1e-11, rtol=1e-11):
            raise ValueError('physical Cq must be Hermitian, not repaired or made real')
        if type(ncells) is not int or ncells < 1 or self.ndof % 3 or np.any(self.masses <= 0) or not np.isfinite(self.masses).all():
            raise ValueError('positive N and full Cartesian Ry masses required')
        if not callable(kernel) or len(self.density_shape) != 2:
            raise ValueError('physical Kq and (grid,component) density shape required')
        self.sqrt_n = np.sqrt(ncells)
        self._validate_edges()
        e = np.asarray(polarizations, dtype=float)
        if e.ndim != 2 or e.shape[1] != 3 or not np.allclose(e @ e.T, np.eye(len(e)), atol=1e-12):
            raise ValueError('real orthonormal Maxwell polarizations required')
        z = array(charges, (self.ndof//3,), 'original valence charges', real=True)
        if np.any(z <= 0) or not np.isfinite(field_mass) or field_mass <= 0:
            raise ValueError('positive charges and Maxwell mass required')
        self.polarizations, self.field_mass = e, field_mass
        self.nfield = len(e) if self.gamma else 0
        self.charge = (z[:, None, None]*np.eye(3)).reshape(self.ndof, 3) @ e.T
        self.mixed = coefficients.mixed @ e.T
        self.contact = e @ coefficients.contact @ e.T
        self.field_vertices = tuple(field_vertices)
        if self.gamma and (len(self.field_vertices) != len(self.sectors)
                           or any(len(v.axes) != 3 for v in self.field_vertices)):
            raise ValueError('three physical real Cartesian L actions per Gamma endpoint required')

    def _validate_edges(self):
        count = len(self.sectors)
        if not count or sorted(e.source for e in self.edges) != list(range(count)) or sorted(e.target for e in self.edges) != list(range(count)):
            raise ValueError('one outgoing and incoming edge per explicit sector required')
        slots, mesh, physical_count = set(), None, 0
        for link in self.edges:
            s, t = self.sectors[link.source], self.sectors[link.target]
            if s.weight != t.weight or s.density_component != t.density_component or not 0 <= s.density_component < self.density_shape[1]:
                raise ValueError('same-spin same-weight directed endpoints required')
            if self.gamma and link.source != link.target:
                raise ValueError('Gamma must preserve each explicit endpoint')
            edge = getattr(link.coupling, 'edge', None)
            if edge is not None:
                physical_count += 1
                if mesh is None:
                    mesh = tuple(edge.mesh)
                if tuple(edge.mesh) != mesh:
                    raise ValueError('all physical edges must use one complete mesh')
                edge.validate(link.coupling.inputs.reference)
                slots.add((edge.ispin-1, tuple(edge.source_basis.slot.mesh_index)))
                if (edge.source_dimension != s.shape[0] or edge.target_dimension != t.shape[0]
                    or edge.source_basis.identity != s.basis_identity or edge.target_basis.identity != t.basis_identity
                    or not np.array_equal(edge.q_fractional, self.q_fractional)
                    or not np.array_equal(edge.q_cart_bohr_inv, self.q_bohr_inv)
                    or int(np.prod(edge.mesh)) != self.ncells
                    or edge.ispin-1 != s.density_component
                    or edge.source_basis.slot.weight != s.weight
                    or edge.target_basis.slot.weight != t.weight):
                    raise ValueError('physical edge and explicit endpoint/q identities differ')
        if physical_count:
            expected = {(spin, index) for spin in range(self.density_shape[1])
                        for index in np.ndindex(mesh)}
            if physical_count != count or count != self.ncells*self.density_shape[1] or slots != expected:
                raise ValueError('every explicit physical spin/full-grid slot required once')

    def validate_stage(self, stage):
        authority(stage.identity, 'current moving stage identity')
        if not np.isfinite(stage.time) or any(len(x) != len(self.sectors) for x in (stage.frames, stage.h_actions, stage.h_adjoints)):
            raise ValueError('complete current endpoint stage required')
        for sector, frame, h, ha in zip(self.sectors, stage.frames, stage.h_actions, stage.h_adjoints):
            array(frame, sector.shape, 'current full-PW frame')
            if not callable(h) or not callable(ha):
                raise TypeError('full-PW H and adjoint callbacks required')

    def tangent(self, ket, bra, u, p, coordinate, momentum):
        if len(ket) != len(self.edges) or len(bra) != len(self.edges):
            raise ValueError('complete independent ket/bra legs required')
        left, right = [], []
        for e, l, b in zip(self.edges, ket, bra):
            ns, rs = self.sectors[e.source].shape
            nt, rt = self.sectors[e.target].shape
            left.append(array(l, (nt, rs), 'directed ket'))
            right.append(array(b, (rt, ns), 'independent directed bra'))
        return JointTangent(tuple(left), tuple(right), array(u, (self.ndof,), 'Uq'),
            array(p, (self.ndof,), 'Pq'), array(coordinate, (self.nfield,), 'global field coordinate'),
            array(momentum, (self.nfield,), 'global field momentum'))

    def zero_tangent(self):
        return self.tangent([np.zeros((self.sectors[e.target].shape[0], self.sectors[e.source].shape[1]), complex) for e in self.edges],
            [np.zeros((self.sectors[e.target].shape[1], self.sectors[e.source].shape[0]), complex) for e in self.edges],
            np.zeros(self.ndof), np.zeros(self.ndof), np.zeros(self.nfield), np.zeros(self.nfield))

    def moments(self, stage, tangent):
        """Weighted n,G,ell of X; D is already in the legs, never reapplied."""
        self.validate_stage(stage)
        density = np.zeros(self.density_shape, complex)
        gradient, ell = np.zeros(self.ndof, complex), np.zeros(3, complex)
        for e, l, b in zip(self.edges, tangent.ket, tangent.bra):
            ps, pt = stage.frames[e.source], stage.frames[e.target]
            w = self.sectors[e.source].weight
            n1, g1 = e.coupling.moments(l, ps, weight=w)
            n2, g2 = e.coupling.moments(pt, b.conj().T, weight=w)
            density += n1+n2
            gradient += g1+g2
            if self.gamma:
                for j, action in enumerate(self.field_vertices[e.source].axes):
                    ell[j] += w*(np.vdot(action.apply(ps), l)
                        + np.einsum('ij,ji->', b, action.adjoint_apply(pt)))
        return density, gradient, ell

    def rhs(self, stage, tangent):
        """Differentiate the settled joint scalar, retaining moving connections."""
        tangent = self.tangent(tangent.ket, tangent.bra, tangent.u, tangent.p,
                               tangent.coordinate, tangent.momentum)
        density, gradient, ell = self.moments(stage, tangent)
        potential = array(self.kernel(density, q_bohr_inv=self.q_bohr_inv), self.density_shape, 'Kq n')
        ket, bra = [], []
        for e, l, b in zip(self.edges, tangent.ket, tangent.bra):
            s, t = self.sectors[e.source], self.sectors[e.target]
            ps, pt = stage.frames[e.source], stage.frames[e.target]
            ionic = e.coupling.prepare(tangent.u.reshape(-1, 3))
            v = potential[:, s.density_component]
            drive = e.coupling.bridge.local_apply(v, ps)+ionic.apply(ps)
            adjoint = e.coupling.bridge.local_adjoint_apply(v, pt)+ionic.adjoint_apply(pt)
            if self.gamma:
                a = self.sqrt_n*self.polarizations.T @ tangent.coordinate
                vertices = self.field_vertices[e.source]
                drive += vertices.apply(a, ps)
                adjoint += vertices.apply(a, pt, adjoint=True)
            ket.append(-1j*(stage.h_actions[e.target](l)+drive @ s.d))
            bra.append(1j*(stage.h_adjoints[e.source](b.conj().T).conj().T+t.d @ adjoint.conj().T))
        velocity = tangent.p/self.masses
        force = -self.curvature @ tangent.u-gradient
        fdot, pidot = np.zeros(0, complex), np.zeros(0, complex)
        if self.gamma:
            velocity = (tangent.p-self.sqrt_n*self.charge @ tangent.coordinate)/self.masses
            force -= self.sqrt_n*self.mixed @ tangent.coordinate
            fdot = tangent.momentum/self.field_mass
            pidot = (-self.sqrt_n*(self.polarizations @ ell+self.mixed.T @ tangent.u
                         - self.charge.T @ (tangent.p/self.masses))
                     - self.ncells*(self.contact+self.charge.T @ (self.charge/self.masses[:, None])) @ tangent.coordinate)
        return JointTangent(tuple(ket), tuple(bra), velocity, force, fdot, pidot)

    def current_variation(self, stage, tangent):
        """Global Hartree current fluctuation at Gamma, separate from noise."""
        if not self.gamma:
            raise ValueError('no independent nonGamma photon/global current variable')
        _, _, ell = self.moments(stage, tangent)
        electronic = -0.5*(self.sqrt_n*(self.polarizations @ ell+self.mixed.T @ tangent.u)
                           + self.ncells*self.contact @ tangent.coordinate)
        ionic = 0.5*(self.sqrt_n*self.charge.T @ (tangent.p/self.masses)
                     - self.ncells*self.charge.T @ ((self.charge @ tangent.coordinate)/self.masses))
        return dict(electronic=electronic, ionic=ionic, total=electronic+ionic)
