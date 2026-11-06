"""Uniform-field ordinary-NC current and contact, in the accepted full PW basis.

Declared family (fixed PW rows, fixed saved local potential and ionic positions):
    H_Ry(a) = |p+a|^2 + V_local + sum_I B_I(a) D_I B_I(a)^dagger,
    B_I(a)[G,c] = exp(-i p_G.tau_I) beta_c(p_G+a), p_G=k+G.
The omitted exp(-i a.tau_I) is a common phase per atom and cancels exactly
inside B D B^dagger, including off-diagonal radial channels. This is uniform
Peierls dressing exp[-i a.(r-r')] of the NC nonlocal kernel. No finite-q
transverse extension, ionic charge, photon mode normalization or bath is implied.

Here a=e A_phys/(hbar c), in bohr^-1 (Gaussian electromagnetic convention),
with e the positive elementary charge: electrons have kinetic p+a. The returned
L=dH_Ry/da has units Ry bohr and C=d2H_Ry/da2 has units Ry bohr^2. The field
energy expansion is H(0)+a_i L_i+1/2 a_i C_ij a_j. The generalized electronic
field force is -L, NOT +L; Hartree atomic velocity is L/2 and charge current
is -L/2 in e*atomic-velocity units. QE commutator_Hx_psi is -i L and
compute_ppsi is L/2. A density/current density additionally needs the state,
k weights and volume; these are deliberately not supplied here.

Both radial and angular projector dependence are differentiated together by
symmetric differences of the EXISTING interpolated factors, not by changing
that Hamiltonian or differentiating the underlying radial quadrature instead.
The exact product rule retains B''DB† + 2B'DB'† + BDB''† in the contact.
Nominal error is O(h^2) only in smooth interpolation intervals. The accepted
forward cubic table is continuous but not generally C1 at its knots; origin
and harmonic branch tolerances can also spoil a differentiable limit. Always
inspect h,h/2,h/4 action refinement. No extrapolation or convergence claim is
made automatically, and a failed refinement is not a physical contact value.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_material_hamiltonian import MaterialHamiltonian
from .qe_nc_projectors import centered_projectors, place_projectors


def _direction(value):
    value = np.asarray(value, dtype=float)
    if value.shape != (3,) or not np.isfinite(value).all():
        raise ValueError("finite Cartesian vector (3,) required")
    return value


def projector_derivatives(evaluate, momenta, direction, step_bohr_inv):
    """Return B,B',B'' for evaluate(p+t*direction), with a declared finite h.

    No normalization of direction: derivatives are with respect to t in
    bohr^-1. Keeping its magnitude makes mixed-contact polarization valid.
    evaluate returns centered (npw,nh) projectors; it must be the same function
    as the declared H(a). Use several steps, not this one stencil in isolation.
    """
    direction = _direction(direction)
    h = float(step_bohr_inv)
    if not np.isfinite(h) or h <= 0:
        raise ValueError("positive finite derivative step in bohr^-1 required")
    center = evaluate(momenta)
    plus = evaluate(momenta + h*direction)
    minus = evaluate(momenta - h*direction)
    return center, (plus-minus)/(2*h), ((plus-center)+(minus-center))/(h*h)


@dataclass(frozen=True)
class ProjectorDerivative:
    """One atom's complete NC product derivative; D may be off diagonal."""
    beta: np.ndarray
    first: np.ndarray
    second: np.ndarray
    active_d: np.ndarray

    def actions(self, vectors):
        """Two full-row actions, without a PW-square matrix or band projection."""
        overlap = zgemm(1., self.beta, vectors, trans_a=2)
        overlap_first = zgemm(1., self.first, vectors, trans_a=2)
        overlap_second = zgemm(1., self.second, vectors, trans_a=2)
        weighted = self.active_d @ overlap
        weighted_first = self.active_d @ overlap_first
        first = zgemm(1., self.first, weighted)
        first += zgemm(1., self.beta, weighted_first)
        second = zgemm(1., self.second, weighted)
        second += zgemm(2., self.first, weighted_first)
        second += zgemm(1., self.beta, self.active_d @ overlap_second)
        return first, second


@dataclass(frozen=True)
class MaterialOpticalCurrent:
    """Reuse an accepted MaterialHamiltonian, its background and its PW rows.

    Additional storage is species-centered derivatives plus one atom's placed
    derivative factors, never three copies of the full species/atom bank.
    Only caller-supplied columns are acted upon, in H0's bounded column blocks.
    """
    hamiltonian: MaterialHamiltonian

    def __post_init__(self):
        if not isinstance(self.hamiltonian, MaterialHamiltonian):
            raise ValueError("an accepted MaterialHamiltonian instance is required")

    @property
    def momenta(self):
        factors = self.hamiltonian.factors
        return (factors.fft.reciprocal_bohr_inv @ factors.fft.miller
                + factors.k_bohr_inv[:, None]).T

    def _vectors(self, vectors):
        vectors = np.asarray(vectors)
        single = vectors.ndim == 1
        block = vectors[:, None] if single else vectors
        if (block.ndim != 2 or block.shape[0] != self.hamiltonian.factors.kinetic_ry.size
                or not np.isfinite(block).all()):
            raise ValueError("finite full-PW vectors (npw,) or (npw,ncolumn) required")
        return block, single

    def field_action(self, vectors, shift_bohr_inv):
        """Exact declared finite H(a) action; H(0) is the accepted H0 verbatim.

        No new DFT or density update, shifted cutoff, occupied-band truncation,
        ultrasoft overlap term, or redefinition of the saved local field.
        """
        shift = _direction(shift_bohr_inv)
        vectors, single = self._vectors(vectors)
        output = self.hamiltonian.apply(vectors)
        if not np.any(shift):
            return output[:, 0] if single else output
        p = self.momenta
        delta_kinetic = 2*(p @ shift) + shift @ shift
        inputs = self.hamiltonian.background.inputs
        for start in range(0, vectors.shape[1], self.hamiltonian.block_size):
            stop = min(start+self.hamiltonian.block_size, vectors.shape[1])
            output[:, start:stop] += delta_kinetic[:, None]*vectors[:, start:stop]
        for species in self.hamiltonian.factors.species:
            tables = inputs.tables[species.name]
            shifted = centered_projectors(tables, p+shift)
            nh = species.active_d.shape[0]
            for slot, atom in enumerate(species.atom_indices):
                beta = place_projectors(shifted, p, inputs.reference.positions_bohr[:, atom-1])
                original = species.beta[:, slot*nh:(slot+1)*nh]
                for start in range(0, vectors.shape[1], self.hamiltonian.block_size):
                    stop = min(start+self.hamiltonian.block_size, vectors.shape[1])
                    block = vectors[:, start:stop]
                    output[:, start:stop] += zgemm(1., beta, species.active_d @ zgemm(1., beta, block, trans_a=2))
                    output[:, start:stop] -= zgemm(1., original, species.active_d @ zgemm(1., original, block, trans_a=2))
        return output[:, 0] if single else output

    def derivative_actions(self, vectors, direction, *, step_bohr_inv):
        """Return (L_e vectors,C_ee vectors) at a=0 for one explicit stencil h.

        Kinetic L=2 p.e and C=2 e.e are exact. Both nonlocal legs, the cross
        contact term and the complete off-diagonal active D remain present.
        Use refinement_actions for the mandatory controlled derivative study.
        """
        e = _direction(direction)
        vectors, single = self._vectors(vectors)
        p = self.momenta
        first = (2*(p @ e))[:, None]*vectors
        first = np.asarray(first, dtype=np.complex128)
        second = np.asarray(2*(e @ e)*vectors, dtype=np.complex128)
        inputs = self.hamiltonian.background.inputs
        for species in self.hamiltonian.factors.species:
            tables = inputs.tables[species.name]
            _, centered_first, centered_second = projector_derivatives(
                lambda momenta: centered_projectors(tables, momenta), p, e, step_bohr_inv)
            nh = species.active_d.shape[0]
            for slot, atom in enumerate(species.atom_indices):
                position = inputs.reference.positions_bohr[:, atom-1]
                term = ProjectorDerivative(
                    species.beta[:, slot*nh:(slot+1)*nh],
                    place_projectors(centered_first, p, position),
                    place_projectors(centered_second, p, position), species.active_d)
                for start in range(0, vectors.shape[1], self.hamiltonian.block_size):
                    stop = min(start+self.hamiltonian.block_size, vectors.shape[1])
                    linear, contact = term.actions(vectors[:, start:stop])
                    first[:, start:stop] += linear
                    second[:, start:stop] += contact
        if single:
            return first[:, 0], second[:, 0]
        return first, second

    def mixed_contact_action(self, vectors, left, right, *, step_bohr_inv):
        """C_ef via polarization, with independent finite-h errors exposed.

        C_ef = [C_(e+f)(e+f)-C_(e-f)(e-f)]/4. Directions are not normalized.
        Repeat at refined h before interpreting a mixed contact numerically.
        """
        left, right = _direction(left), _direction(right)
        _, plus = self.derivative_actions(vectors, left+right, step_bohr_inv=step_bohr_inv)
        _, minus = self.derivative_actions(vectors, left-right, step_bohr_inv=step_bohr_inv)
        return (plus-minus)/4

    def refinement_actions(self, vectors, direction, *, step_bohr_inv):
        """Yield h,h/2,h/4 actions; no silently accepted or extrapolated limit."""
        for divisor in (1, 2, 4):
            h = float(step_bohr_inv)/divisor
            first, second = self.derivative_actions(vectors, direction, step_bohr_inv=h)
            yield h, first, second
