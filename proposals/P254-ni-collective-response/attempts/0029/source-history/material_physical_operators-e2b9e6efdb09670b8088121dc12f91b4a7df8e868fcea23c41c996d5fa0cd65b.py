"""Unanchored physical-UPF factors shared by H0, spatial g and q0 L.

The saved local potential, K, original atoms, basis and complex state D are not
changed. No legacy interpolated projector bank is allocated in construction.
See P254 A28/physical-spatial-registration.json for the source-bound contract.

Scalar cutover: use these H0/g operators instead of A25+DeltaV0/Delta_g.
Keep PreparedReference.C and the existing fixed DeltaC coefficient OR install
Cphysical(q)=Ccommon(q)+DeltaC at every matter q and pass zero DeltaC to the
joint coordinate correction, never both. Contact and mixed remain fixed scalar
coefficients. PhysicalFieldVertices is the field-only vertex for either route;
UPFPreparedVertices on A25 is an independent comparison, not a runtime addition.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_bloch_frame import BlochBasis
from .qe_real_bloch_grid import RealBlochBasis
from .qe_bloch_dyads import BlochDyadCoupling
from .qe_bloch_finiteq_ionic import BlochFiniteQIonic
from .qe_material_grid import PhysicalRhoFFT
from .qe_material_hamiltonian import PWFactors, SpeciesProjectors, MaterialHamiltonian
from .qe_material_field_derivatives import UPFRadial, centered_derivatives
from .qe_nc_projectors import place_projectors
from .material_mixed_field_ionic import first_cartesian_jets
from .material_frozen_action import FrozenMaterialAction
from .material_joint_field_scalar import JointFieldCoefficients


REPRESENTATION = 'direct-physical-UPF-zero-field-v1'


def build_physical_factors(inputs, basis):
    """Physical factors for checked unitary or real-collinear original bases."""
    if not isinstance(basis, (BlochBasis, RealBlochBasis)):
        raise ValueError('checked BlochBasis or RealBlochBasis required')
    ref = inputs.reference
    basis.validate(ref)
    real_basis = isinstance(basis, RealBlochBasis)
    if real_basis:
        actual = {(name, table.upf.sha256) for name, table in inputs.tables.items()}
        expected = {(s.name, next(x['sha256'] for x in basis.scope.upf_sources
                     if x['path'].endswith('/'+s.pseudo_file))) for s in ref.species}
        if actual != expected:
            raise ValueError('physical factor UPFs differ from admitted real-Hamiltonian scope')
    volume = abs(float(np.linalg.det(ref.cell_bohr)))
    fft = PhysicalRhoFFT.create(ref.fft_smooth, basis.miller,
                               basis.reciprocal_bohr_inv, volume)
    p = (basis.reciprocal_bohr_inv @ basis.miller+basis.k_cart_bohr_inv[:, None]).T
    groups = []
    for species in ref.species:
        table = inputs.tables[species.name]
        if table.omega_bohr3 != volume:
            raise ValueError('original radial volume must equal physical cell volume')
        center = centered_derivatives(UPFRadial(table), p, np.zeros(3))[0]
        atoms = tuple(i for i, name in enumerate(ref.atom_names) if name == species.name)
        nh = table.upf.nh
        beta = np.empty((len(p), len(atoms)*nh), complex, order='F')
        for slot, atom in enumerate(atoms):
            beta[:, slot*nh:(slot+1)*nh] = place_projectors(center, p, ref.positions_bohr[:, atom])
        groups.append(SpeciesProjectors(species.name, tuple(a+1 for a in atoms),
                                        beta, table.upf.expanded_d))
    kind = 'real-collinear-orbital-conjugation-basis' if real_basis else 'unitary-reconstructed-basis'
    provenance = dict(kind=kind, identity=basis.identity,
                      source=dict(basis.slot.source_binding), xml_source=dict(basis.reference_source),
                      ispin=basis.slot.ispin, mesh_index=basis.slot.mesh_index,
                      source_ik=basis.slot.source_ik, operation=basis.slot.operation.index,
                      gauge_policy=basis.slot.gauge_policy, representation=REPRESENTATION)
    if real_basis:
        provenance.update(scope=basis.scope.identity, orbital_conjugation=basis.slot.orbital_conjugation)
    return PWFactors(fft, np.einsum('pi,pi->p', p, p), tuple(groups),
                     basis.k_cart_bohr_inv, None, inputs.sources+(provenance,))


def physical_hamiltonian(background, basis, ispin, block_size=8):
    if ispin != basis.slot.ispin:
        raise ValueError('explicit spin must match the checked basis slot')
    return MaterialHamiltonian(background, build_physical_factors(background.inputs, basis),
                               ispin, block_size)


def physical_coupling(inputs, basis, factors, block_size=8):
    """Gamma/q=0 coupling, sharing factors by object identity with physical H0."""
    _require_physical(factors)
    return BlochDyadCoupling(inputs, basis, factors=factors, block_size=block_size)


def physical_finiteq(inputs, edge, source_factors, target_factors, block_size=8):
    """Rectangular coupling with both actual endpoint physical banks injected."""
    _require_physical(source_factors)
    _require_physical(target_factors)
    return BlochFiniteQIonic(inputs, edge, source_factors=source_factors,
                            target_factors=target_factors, block_size=block_size)


def _require_physical(factors):
    if not factors.sources or factors.sources[-1].get('representation') != REPRESENTATION:
        raise ValueError('direct physical-UPF factors required; A25 correction must not be reapplied')


def _direction(value):
    a = np.asarray(value)
    if a.shape != (3,) or np.iscomplexobj(a) or not np.isfinite(a).all():
        raise ValueError('finite real Cartesian field shift/direction (3,) required')
    return np.asarray(a, dtype=float)


def _vectors(value, rows):
    v = np.asarray(value)
    single = v.ndim == 1
    x = v[:, None] if single else v
    if x.ndim != 2 or len(x) != rows or not np.isfinite(x).all():
        raise ValueError('finite full-PW vectors required')
    return x, single


def _weighted_overlap(species, overlap, width):
    shape = (len(species.atom_indices), species.active_d.shape[0], width)
    return np.einsum('ij,ajb->aib', species.active_d,
                     overlap.reshape(shape), optimize=True).reshape(-1, width)


@dataclass(frozen=True)
class PhysicalFieldAction:
    """L(direction), not H(a)-H_A25 and not a rho-dependent contact."""
    momenta: np.ndarray
    direction: np.ndarray
    terms: tuple
    block_size: int

    def apply(self, vectors):
        x, single = _vectors(vectors, len(self.momenta))
        out = np.empty(x.shape, complex)
        kinetic = 2*(self.momenta @ self.direction)
        for start in range(0, x.shape[1], self.block_size):
            sl = slice(start, start+self.block_size)
            v = x[:, sl]
            value = np.asarray(kinetic[:, None]*v, dtype=complex)
            for species, derivative in self.terms:
                base = _weighted_overlap(species, zgemm(1., species.beta, v, trans_a=2), v.shape[1])
                first = _weighted_overlap(species, zgemm(1., derivative, v, trans_a=2), v.shape[1])
                value += zgemm(1., derivative, base)+zgemm(1., species.beta, first)
            out[:, sl] = value
        return out[:, 0] if single else out

    def adjoint_apply(self, vectors):
        # Real field directions and original Hermitian UPF D make L Hermitian.
        return self.apply(vectors)


class PhysicalFieldVertices:
    """Physical first-field actions and unweighted original-complex-D reductions.

    evaluate has the existing JointPreparedFieldDynamics vertex interface but
    returns exactly zero zero-field corrections: H0 and g are already physical.
    prepare caches one directional derivative bank; no bank for exact zero a.
    No magnitude cutoff and no second-coordinate operator is used.
    """
    def __init__(self, hamiltonian, *, projector_block_size=64):
        _require_physical(hamiltonian.factors)
        if type(projector_block_size) is not int or projector_block_size < 1:
            raise ValueError('positive projector block size required')
        self.hamiltonian = hamiltonian
        f = hamiltonian.factors
        self.momenta = (f.fft.reciprocal_bohr_inv @ f.fft.miller+f.k_bohr_inv[:, None]).T
        self.block_size = projector_block_size
        for species in f.species:
            if not np.array_equal(species.active_d, species.active_d.conj().T):
                raise ValueError('original Hermitian UPF channel D required')

    def prepare(self, direction):
        a = _direction(direction).copy()
        a.setflags(write=False)
        terms = []
        inputs = self.hamiltonian.background.inputs
        if np.any(a != 0):
            for species in self.hamiltonian.factors.species:
                center = centered_derivatives(UPFRadial(inputs.tables[species.name]), self.momenta, a)[1]
                nh = species.active_d.shape[0]
                derivative = np.empty_like(species.beta, order='F')
                for slot, atom in enumerate(species.atom_indices):
                    derivative[:, slot*nh:(slot+1)*nh] = place_projectors(
                        center, self.momenta, inputs.reference.positions_bohr[:, atom-1])
                derivative.setflags(write=False)
                terms.append((species, derivative))
        return PhysicalFieldAction(self.momenta, a, tuple(terms), self.block_size)

    def moments(self, frame, density_coefficients):
        """Return <L_x,L_y,L_z> and current=-<L>/2, without 3 action frames."""
        psi, _ = _vectors(frame, len(self.momenta))
        d = np.asarray(density_coefficients)
        if d.shape != (psi.shape[1],)*2 or not np.isfinite(d).all():
            raise ValueError('original finite complex raw-frame D required')
        diagonal = np.zeros(len(psi), complex)
        for start in range(0, psi.shape[1], self.block_size):
            sl = slice(start, start+self.block_size)
            weighted = zgemm(1., psi, d[:, sl])
            diagonal += np.einsum('pi,pi->p', psi[:, sl].conj(), weighted)
        gradient = 2*self.momenta.T @ diagonal
        inputs = self.hamiltonian.background.inputs
        for species in self.hamiltonian.factors.species:
            jets = first_cartesian_jets(UPFRadial(inputs.tables[species.name]), self.momenta)
            x = zgemm(1., species.beta, psi, trans_a=2)
            dx = _weighted_overlap(species, x, psi.shape[1])
            nh = species.active_d.shape[0]
            derivative = np.empty_like(species.beta, order='F')
            weighted_base = dx @ d
            for axis in range(3):
                for slot, atom in enumerate(species.atom_indices):
                    derivative[:, slot*nh:(slot+1)*nh] = place_projectors(
                        jets[axis+1], self.momenta, inputs.reference.positions_bohr[:, atom-1])
                xa = zgemm(1., derivative, psi, trans_a=2)
                dxa = _weighted_overlap(species, xa, psi.shape[1])
                gradient[axis] += np.vdot(xa, weighted_base)+np.vdot(x, dxa @ d)
        return dict(field_gradient=gradient, electronic_current_hartree=-gradient/2)

    def evaluate(self, frame, density_coefficients, displacement, shift):
        ndof = 3*self.hamiltonian.background.inputs.reference.nat
        u = np.asarray(displacement)
        if u.shape not in ((ndof,), (ndof//3, 3)) or not np.isfinite(u).all():
            raise ValueError('finite original atom-order displacement required')
        values = self.moments(frame, density_coefficients)
        return dict(action=self.prepare(shift).apply(frame), delta_energy=0j,
                    delta_gradient=np.zeros(ndof, complex), **values)


class PhysicalFrozenAction:
    """Existing fused Hphysical+V+gphysical(u), plus ONLY L(a)."""
    def __init__(self, h0, coupling=None, *, potential=None, direction=None,
                 shift=(0., 0., 0.), fft_block_size=8, projector_block_size=64):
        self.frozen = FrozenMaterialAction(h0, coupling, potential=potential,
            direction=direction, fft_block_size=fft_block_size,
            projector_block_size=projector_block_size)
        a = _direction(shift)
        vertices = PhysicalFieldVertices(h0, projector_block_size=projector_block_size)
        self.field = vertices.prepare(a) if np.any(a != 0) else None

    def apply(self, vectors):
        out = self.frozen.apply(vectors)
        if self.field is not None:
            out += self.field.apply(vectors)
        return out

    def adjoint_apply(self, vectors):
        out = self.frozen.adjoint_apply(vectors)
        if self.field is not None:
            out += self.field.adjoint_apply(vectors)
        return out


def spatial_coefficient_handoff(common_q, coefficients):
    """Raw complex C(q) corrected ONCE and field-only residual coefficients.

    common_q has any mesh prefix followed by (ndof,ndof). Returns new arrays;
    never mutates existing references, projects small imaginary parts, applies
    ASR, or changes source state. For PreparedReference's strictly real C API,
    this raw array is not silently coerced: keep its existing C and pass the
    original coefficients to JointPreparedFieldDynamics instead. That route
    carries the same fixed DeltaC exactly once without altering the real API.
    """
    c = np.asarray(common_q)
    delta = coefficients.delta_curvature
    if c.ndim < 2 or c.shape[-2:] != delta.shape or not np.isfinite(c).all():
        raise ValueError('complete original C(q) with final full coordinate axes required')
    corrected = c+delta
    residual = JointFieldCoefficients(np.zeros_like(delta), coefficients.contact,
                                      coefficients.mixed, coefficients.coefficient_authority)
    return corrected, residual
