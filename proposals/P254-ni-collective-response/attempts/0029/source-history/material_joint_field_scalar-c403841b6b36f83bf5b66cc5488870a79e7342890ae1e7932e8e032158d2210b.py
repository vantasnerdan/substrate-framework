"""ONE conditional joint-second-order scalar; see A28 preregistration.

E_total=N*(E_PreparedReference + <DeltaV0> + u.<Delta_g> + a.<L>
 + u.DeltaC.u/2 + a.Cfield.a/2 + u.T.a
 + sum[(p-Za)^2-p^2]/(2M)) + Pi^2/(2Mfield).

H0_A25 remains named and unchanged; DeltaV0 is the physical UNANCHORED
UPF zero-field correction. Cfield,DeltaC,T are fixed prepared coefficients,
not rho-dependent operators. No delta_rho*a*a, delta_rho*u*a or
 delta_rho*u*u term is silently promoted to an all-field model.
"""
from dataclasses import dataclass

import numpy as np

from .prepared_reference import array, frozen, authority
from .prepared_factor_dynamics import FactorMean
from .physical_field_environment import UniformMaxwellField


@dataclass(frozen=True)
class JointFieldCoefficients:
    delta_curvature: np.ndarray
    contact: np.ndarray
    mixed: np.ndarray
    coefficient_authority: str

    def __post_init__(self):
        authority(self.coefficient_authority, 'fixed prepared coefficient authority')
        t = np.asarray(self.mixed)
        if t.ndim != 2 or t.shape[1] != 3 or t.shape[0] % 3:
            raise ValueError('mixed coefficient must have (3*nat,3) shape')
        for name, shape in [('mixed', t.shape), ('contact', (3,3)),
                            ('delta_curvature', (len(t),len(t)))]:
            object.__setattr__(self, name, frozen(array(getattr(self,name), shape, name)))
        for value in (self.contact, self.delta_curvature):
            if not np.allclose(value, value.T, rtol=1e-12, atol=1e-12):
                raise ValueError('symmetric fixed Taylor Hessian required; no repair')


def onsite_matrix(onsite):
    blocks = np.asarray(onsite)
    if blocks.ndim != 3 or blocks.shape[1:] != (3,3):
        raise ValueError('one complete onsite block per original atom required')
    result = np.zeros((3*len(blocks),3*len(blocks)), dtype=blocks.dtype)
    for index, block in enumerate(blocks):
        result[3*index:3*index+3,3*index:3*index+3] = block
    return result


def coordinate_correction(coefficients, displacement, momentum, shift, masses,
                          charges, delta_gradient, field_gradient):
    """Per-cell scalar addition and its exact derivatives; p is canonical.

    The return kinetic correction subtracts the neutral kinetic already in
    PreparedReference, never adds a second ionic kinetic term.
    """
    ndof = len(coefficients.mixed)
    u = array(displacement, (ndof,), 'u', real=True)
    p = array(momentum, (ndof,), 'canonical p', real=True)
    a = array(shift, (3,), 'a', real=True)
    m = array(masses, (ndof,), 'Ry masses', real=True)
    z = array(charges, (ndof//3,), 'original UPF valence charges', real=True)
    if np.any(m <= 0) or np.any(z <= 0):
        raise ValueError('positive original masses and valence charges required')
    dg = array(delta_gradient, (ndof,), 'delta positive gradient')
    ell = array(field_gradient, (3,), 'physical electronic field gradient')
    za = (z[:, None]*a).reshape(-1)
    velocity = (p-za)/m
    ion_field_gradient = -np.sum(z[:, None]*velocity.reshape(-1,3), axis=0)
    electron_field_gradient = ell+coefficients.contact @ a+coefficients.mixed.T @ u
    energy = (u @ dg + a @ ell + u @ coefficients.delta_curvature @ u/2
              + a @ coefficients.contact @ a/2 + u @ coefficients.mixed @ a
              + np.sum(((p-za)**2-p**2)/m)/2)
    return dict(energy=energy, velocity=velocity,
                force_correction=-dg-coefficients.delta_curvature @ u-coefficients.mixed @ a,
                electronic_field_gradient=electron_field_gradient,
                ionic_field_gradient=ion_field_gradient,
                field_gradient=electron_field_gradient+ion_field_gradient,
                electronic_current_hartree=-electron_field_gradient/2,
                ionic_current_hartree=-ion_field_gradient/2)


@dataclass(frozen=True)
class JointFieldMean:
    matter: FactorMean
    coordinate: np.ndarray
    momentum: np.ndarray


class JointPreparedFieldDynamics:
    """Compose existing PreparedFactorDynamics without changing H0/K modules.

    Homogeneous mean u,p: Poisson bracket 1/N. Independent electronic slot
    bracket 1/(N*w); field bracket 1. Thus E_total below generates all RHS.
    Spatially inhomogeneous ions/finite-q photons are outside this API.
    Held sectors require separate physical vertices; this class rejects them.
    """
    def __init__(self, matter, vertices, coefficients, field, charges,
                 *, scalar_authority):
        authority(scalar_authority, 'joint scalar authority')
        self.matter, self.vertices = matter, tuple(vertices)
        self.coefficients, self.field = coefficients, field
        self.scalar_authority = scalar_authority
        r = matter.reference
        if (not isinstance(field, UniformMaxwellField) or len(self.vertices) != len(matter.sectors)
                or coefficients.mixed.shape != (r.ndof,3)):
            raise ValueError('matching actual sectors, full coefficients and q0 field required')
        if np.any(r.held_density) or np.any(r.held_gradient):
            raise ValueError('held physical electronic vertices unavailable; empty held sector required')
        cell_volume = matter.sectors[0].coupling.cell_volume_bohr3
        if not np.isclose(field.volume_bohr3, r.ncells*cell_volume, rtol=1e-12, atol=0):
            raise ValueError('Maxwell volume must match N*admitted matter cell volume')
        self.charges = frozen(array(charges, (r.ndof//3,), 'UPF valence charges', real=True))
        if np.any(self.charges <= 0):
            raise ValueError('positive original valence charges required')

    def mean(self, frames, u, p, coordinate, momentum):
        return JointFieldMean(self.matter.mean(frames,u,p),
            array(coordinate, (self.field.size,), 'q0 coordinate', real=True),
            array(momentum, (self.field.size,), 'q0 momentum', real=True))

    def evaluate(self, state):
        r = self.matter.reference
        q = array(state.coordinate, (self.field.size,), 'field coordinate', real=True)
        pi = array(state.momentum, (self.field.size,), 'field momentum', real=True)
        a = self.field.cartesian_shift(q)
        base = self.matter.evaluate_mean(state.matter)
        dg, ell, de = np.zeros(r.ndof,complex), np.zeros(3,complex), 0j
        hframes = []
        for s, vertex, psi, hpsi in zip(self.matter.sectors,self.vertices,
                                       state.matter.frames,base.h_frames):
            values = vertex.evaluate(psi,s.d,state.matter.u,a)
            dg += s.weight*values['delta_gradient']
            ell += s.weight*values['field_gradient']
            de += s.weight*values['delta_energy']
            hframes.append(hpsi+values['action'])
        correction = coordinate_correction(self.coefficients,state.matter.u,
            state.matter.p,a,r.masses,self.charges,dg,ell)
        force = base.rhs.p+correction['force_correction']
        matter_rhs = FactorMean(tuple(-1j*h for h in hframes),correction['velocity'],force)
        field_rhs = (pi/self.field.mass,
                     -r.ncells*self.field.polarizations @ correction['field_gradient'])
        return dict(energy=r.ncells*(base.energy+de+correction['energy'])+self.field.energy(pi),
                    energy_per_cell_without_maxwell=base.energy+de+correction['energy'],
                    matter_rhs=matter_rhs, field_rhs=field_rhs, h_frames=tuple(hframes),
                    electronic_current_hartree=r.ncells*correction['electronic_current_hartree'],
                    ionic_current_hartree=r.ncells*correction['ionic_current_hartree'],
                    electric_field_hartree=self.field.electric_field_hartree(pi),
                    scalar_authority=self.scalar_authority,
                    scope='joint-second-order q0; no propagated noise or all-field gauge claim')


class MissingFieldPreparation(ValueError):
    """A normalizable Gaussian requires explicitly supplied positive stiffness."""


def conditional_field_gaussian(field, contact, masses, charges, ncells, *,
                               state_authority, reality_relative_tolerance=1e-10):
    """Centered trap vacuum released under the AFFINE joint scalar, not equilibrium.

    Ionic canonical momenta, not velocities, are clamped; their diamagnetic
    stiffness must therefore be retained. No eigenvalue clipping or damping.
    """
    authority(state_authority, 'conditional product preparation hypothesis')
    c = array(contact,(3,3),'fixed prepared contact')
    z = np.asarray(charges)
    m = array(masses,(3*len(z),),'Cartesian Ry masses',real=True).reshape(-1,3)
    z = array(z,(len(z),),'original valence charges',real=True)
    if (type(ncells) is not int or ncells < 1 or np.any(m <= 0) or np.any(z <= 0)
            or not 0 <= reality_relative_tolerance < 1):
        raise ValueError('positive cell count, masses/charges and explicit reality tolerance required')
    raw = ncells*field.polarizations @ (c+np.diag(np.sum(z[:,None]**2/m,axis=0))) @ field.polarizations.T
    scale = max(1.,float(np.max(np.abs(raw))))
    imaginary = float(np.max(np.abs(raw.imag)))
    asymmetry = float(np.max(np.abs(raw-raw.T)))
    if imaginary > reality_relative_tolerance*scale or asymmetry > reality_relative_tolerance*scale:
        raise MissingFieldPreparation('clamped stiffness is not real symmetric within declared roundoff tolerance; supply physical field preparation, no Hessian repair')
    # Real extraction only after reported roundoff test; no symmetry averaging.
    values, vectors = np.linalg.eigh(raw.real)
    if np.any(values <= 0):
        raise MissingFieldPreparation('clamped q0 stiffness is not strictly positive: eigenvalues='+repr(values.tolist())+'; supply an explicit external preparation/trap or normalizable field state, not clipping')
    qq = (vectors*(1/(2*np.sqrt(field.mass*values)))) @ vectors.T
    pp = (vectors*(np.sqrt(field.mass*values)/2)) @ vectors.T
    cross = .5j*np.eye(field.size)
    ordered = np.block([[qq,cross],[-cross,pp]])
    return dict(coordinate_mean=np.zeros(field.size), momentum_mean=np.zeros(field.size),
                clamped_stiffness_raw=raw, stiffness_eigenvalues=values,
                stiffness_imaginary_max=imaginary, stiffness_asymmetry_max=asymmetry,
                reality_relative_tolerance=reality_relative_tolerance,
                covariance_qq=qq,covariance_pp=pp,ordered_covariance=ordered,
                symmetrized_covariance=np.block([[qq,np.zeros_like(qq)],[np.zeros_like(pp),pp]]),
                frequencies_ry=np.sqrt(values/field.mass), state_authority=state_authority,
                quadratic_zero_point_energy_ry=float(np.sum(np.sqrt(values/field.mass))/2),
                hypothesis='centered sudden-release electronic/ionic/field tensor product; initial linear drift retained',
                coupled_equilibrium=False, propagated_noise=False, finite_q_photons=False)
