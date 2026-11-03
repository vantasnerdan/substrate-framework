"""Finite periodic q=0 Maxwell coordinate and the physical UPF optical family.

Coordinates a=e A_SI/hbar (bohr^-1), tau=t Ry/hbar, hbar=1. The
SI Maxwell energy, Ry=e^2/(8*pi*epsilon0*a0), gives M_a=Omega/(8*pi):
H_EM=Pi.Pi/(2*M_a), adot=Pi/M_a, E_Ha=-adot/2. There is no magnetic
restoring term for a uniform periodic field and no inferred cavity frequency.
The three harmonic torus coordinates are divergence-free; no propagation
wavevector or finite-q matter coupling is implied. Polarizations select an
explicit orthonormal subspace. Quantum [a_i,Pi_j]=i delta_ij is compatible
with this Hamiltonian, but no quantum field state/vacuum is inferred (the
free q=0 coordinate has no normalizable oscillator vacuum).

For electronic H in Ry, J_Ha=-dH/da/2. Ionic canonical P has kinetic
(P-Z*a)^2/(2*M_Ry), M_Ry=m/(2*m_e). All gradients, forces, and currents
must be derivatives of ONE scalar. A finite electronic/ionic environment
can exchange energy and recur; this module adds no irreversible loss.
Sources and symbolic unit derivation: A28 physical-environment-registration.
"""
from dataclasses import dataclass
import numpy as np

from .prepared_reference import array, authority, frozen
from .qe_material_field_derivatives import (
    MaterialFieldDerivatives, centered_derivatives, product_actions,
)
from .qe_nc_projectors import place_projectors


@dataclass(frozen=True)
class UniformMaxwellField:
    volume_bohr3: float
    polarizations: np.ndarray
    cell_authority: str

    def __post_init__(self):
        authority(self.cell_authority, 'finite periodic cell authority')
        if not np.isfinite(self.volume_bohr3) or self.volume_bohr3 <= 0:
            raise ValueError('positive finite physical volume required')
        e = np.asarray(self.polarizations)
        if e.ndim != 2 or e.shape[1] != 3 or not 1 <= len(e) <= 3:
            raise ValueError('one to three explicitly chosen Cartesian polarizations required')
        e = array(e, e.shape, 'polarizations', real=True)
        if not np.allclose(e @ e.T, np.eye(len(e)), rtol=1e-13, atol=1e-13):
            raise ValueError('orthonormal polarizations required; no normalization repair')
        object.__setattr__(self, 'polarizations', frozen(e))

    @classmethod
    def from_cell(cls, cell_bohr, *, polarizations, cell_authority):
        cell = array(cell_bohr, (3, 3), 'admitted cell', real=True)
        return cls(float(abs(np.linalg.det(cell))), polarizations, cell_authority)

    @property
    def mass(self):
        return self.volume_bohr3 / (8*np.pi)

    @property
    def size(self):
        return len(self.polarizations)

    def cartesian_shift(self, coordinate):
        return array(coordinate, (self.size,), 'field coordinate', real=True) @ self.polarizations

    def energy(self, momentum):
        p = array(momentum, (self.size,), 'field canonical momentum', real=True)
        return float(p @ p / (2*self.mass))

    def rhs(self, momentum, matter_gradient_cartesian):
        p = array(momentum, (self.size,), 'field canonical momentum', real=True)
        g = array(matter_gradient_cartesian, (3,), 'dE_matter/da', real=True)
        return p/self.mass, -self.polarizations @ g

    def electric_field_hartree(self, momentum):
        p = array(momentum, (self.size,), 'field canonical momentum', real=True)
        return -(p @ self.polarizations)/(2*self.mass)

    def normalization(self):
        return dict(volume_bohr3=self.volume_bohr3, mass=self.mass,
                    polarizations=self.polarizations.tolist(), cell_authority=self.cell_authority,
                    coordinate='a=e A_SI/hbar, numerical bohr^-1',
                    time='tau=t Ry/hbar', energy='Ry',
                    momentum='Pi=M_a da/dtau, canonical conjugate of a',
                    electric_field='Hartree atomic units: E=-da/dtau/2',
                    scope='uniform periodic q0 Maxwell coordinate; zero magnetic stiffness',
                    field_state=None, finite_q_matter_bridge=False)


def charged_ion_kinematics(momentum, masses_ry, charges_e, shift_bohr_inv):
    """Exact uniform minimal coupling, no charge/mass/ionic state defaults."""
    p = np.asarray(momentum)
    if p.ndim != 2 or p.shape[1] != 3 or len(p) < 1:
        raise ValueError('complete ionic canonical momenta (nat,3) required')
    p = array(p, p.shape, 'ion momentum', real=True)
    m = array(masses_ry, (len(p),), 'Ry ionic masses', real=True)
    z = array(charges_e, (len(p),), 'ionic charges in positive e', real=True)
    a = array(shift_bohr_inv, (3,), 'uniform shift', real=True)
    if np.any(m <= 0):
        raise ValueError('positive Ry ionic masses required')
    mechanical = p-z[:, None]*a
    velocity = mechanical/m[:, None]
    return dict(velocity=velocity, kinetic_energy=float(np.sum(mechanical*velocity)/2),
                field_gradient=-np.sum(z[:, None]*velocity, axis=0),
                field_contact=np.eye(3)*np.sum(z*z/m),
                current_hartree=np.sum(z[:, None]*velocity, axis=0)/2)


def joint_canonical_rhs(field, coordinate, field_momentum, *, electronic_gradient,
                        ion_momentum, masses_ry, charges_e, conservative_ionic_force,
                        scalar_authority):
    """Compose conservative field/ionic Hamilton equations at supplied matter state.

    electronic_gradient is the EXTENSIVE dE_e/da, with sector measures included
    exactly once by the caller. conservative_ionic_force=-dE_potential/dR
    includes electronic and ionic potential forces at THIS field and position;
    it is not licensed by a zero-field gradient checkpoint alone. Electronic
    frames obey i Psidot=H(a,R,rho) Psi for the SAME variational scalar, D fixed.
    Then dE_e/dtau=g_e.adot-F_e.Rdot, ionic potential/kinetic terms cancel
    the Rdot exchange and dH_EM/dtau=-(g_e+g_ion).adot. No numerical
    integrator or arbitrary supplied-force consistency is asserted here.
    """
    authority(scalar_authority, 'common electronic/ionic/field scalar authority')
    ions = charged_ion_kinematics(ion_momentum, masses_ry, charges_e,
                                 field.cartesian_shift(coordinate))
    ge = array(electronic_gradient, (3,), 'extensive electronic gradient', real=True)
    force = array(conservative_ionic_force, ions['velocity'].shape,
                  'same-scalar ionic force', real=True)
    adot, pdot = field.rhs(field_momentum, ge+ions['field_gradient'])
    cartesian_adot = adot @ field.polarizations
    matter_power = float((ge+ions['field_gradient']) @ cartesian_adot)
    return dict(field_coordinate=adot, field_momentum=pdot,
                ion_position=ions['velocity'], ion_momentum=force,
                field_power=-matter_power, matter_field_power=matter_power,
                conservative_potential_power=-float(np.sum(force*ions['velocity'])),
                scalar_authority=scalar_authority)


class UPFUniformCoupling:
    """Analytic action/current/contact at any uniform a of ONE physical family.

    Unlike the original optional anchored family this ALWAYS uses unanchored
    physical UPF H(0). A constant undressed A25 offset is never added. Exact
    radial quadrature or a caller-declared Hermite radial approximation is
    inherited from MaterialFieldDerivatives. Finite-PW truncation itself does
    not confer an exact continuum Ward identity or finite-q transversality.
    """
    def __init__(self, derivatives):
        if not isinstance(derivatives, MaterialFieldDerivatives):
            raise TypeError('physical MaterialFieldDerivatives required')
        self.derivatives = derivatives

    def action(self, vectors, shift_bohr_inv):
        return self.derivatives.field_action(vectors, shift_bohr_inv, anchored=False)

    def derivative_actions(self, vectors, shift_bohr_inv, left, right=None):
        f = self.derivatives
        vectors, single = f.original._vectors(vectors)
        shift = array(shift_bohr_inv, (3,), 'uniform shift', real=True)
        e = array(left, (3,), 'left Cartesian derivative', real=True)
        r = e if right is None else array(right, (3,), 'right Cartesian derivative', real=True)
        shifted = f.momenta+shift
        linear = np.asarray(2*(shifted @ e)[:, None]*vectors, dtype=complex)
        contact = np.asarray(2*(e @ r)*vectors, dtype=complex)
        inputs = f.hamiltonian.background.inputs
        for species in f.hamiltonian.factors.species:
            centered = centered_derivatives(f.radial[species.name], shifted, e, r)
            for atom in species.atom_indices:
                # Uniform Peierls phase cancels in each projector dyad. Keep
                # placement at unshifted p, exactly as physical field_action.
                placed = tuple(place_projectors(b, f.momenta,
                    inputs.reference.positions_bohr[:, atom-1]) for b in centered)
                for start in range(0, vectors.shape[1], f.hamiltonian.block_size):
                    stop = min(start+f.hamiltonian.block_size, vectors.shape[1])
                    l, c = product_actions(*placed, species.active_d, vectors[:, start:stop])
                    linear[:, start:stop] += l
                    contact[:, start:stop] += c
        return (linear[:, 0], contact[:, 0]) if single else (linear, contact)

    def current_action(self, vectors, shift_bohr_inv, direction):
        return -self.derivative_actions(vectors, shift_bohr_inv, direction)[0]/2

    def moments(self, frame, density_coefficients, shift_bohr_inv, *, block_size=8,
                row_block_size=2048):
        """Full raw-column unweighted energy/current/contact, no band truncation.

        D must be the admitted preparation's coefficient matrix, not population
        labels. Imaginary diagnostics survive; callers may reject, not silently
        discard, non-real conservative gradients. Peak action storage is O(P*b).
        """
        psi = np.asarray(frame)
        if (psi.ndim != 2 or psi.shape[0] != len(self.derivatives.momenta)
                or not 0 < psi.shape[1] <= psi.shape[0]):
            raise ValueError('complete physical PW frame required')
        d = array(density_coefficients, (psi.shape[1], psi.shape[1]), 'prepared D')
        if type(block_size) is not int or not 1 <= block_size <= 16:
            raise ValueError('column block in 1..16 required')
        if type(row_block_size) is not int or row_block_size < 1:
            raise ValueError('positive row block required')
        def pair(action, start, stop):
            value = 0j
            # Tr(D Psi† A Psi); restrict only the output column leg, never D.
            for row in range(0, len(psi), row_block_size):
                end = min(row+row_block_size, len(psi))
                weighted = psi[row:end] @ d[start:stop, :].conj().T
                value += np.vdot(weighted, action[row:end])
            return value
        energy, gradient, contact = 0j, np.zeros(3, complex), np.zeros((3, 3), complex)
        axes = np.eye(3)
        for start in range(0, psi.shape[1], block_size):
            stop = min(start+block_size, psi.shape[1])
            vectors = psi[:, start:stop]
            energy += pair(self.action(vectors, shift_bohr_inv), start, stop)
            for i in range(3):
                for j in range(i, 3):
                    l, c = self.derivative_actions(vectors, shift_bohr_inv, axes[i], axes[j])
                    if i == j:
                        gradient[i] += pair(l, start, stop)
                    contact[i, j] += pair(c, start, stop)
                    if i != j:
                        contact[j, i] = contact[i, j]  # exact mixed-partial identity
        return dict(energy_ry=energy, gradient=gradient, current_hartree=-gradient/2,
                    contact=contact, normalization='one unweighted spin/k sector',
                    family='unanchored physical UPF; no constant A25 offset')
