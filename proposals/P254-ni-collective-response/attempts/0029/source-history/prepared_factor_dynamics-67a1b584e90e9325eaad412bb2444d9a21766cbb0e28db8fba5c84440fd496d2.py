"""Full-row moving factor mean and independent complex ket/bra tangent.

rho=Psi D Psi^dagger, Psi_dot=-i H Psi, D constant (no gap divisions).
A tangent is delta_rho=L Psi^dagger + Psi R with INDEPENDENT L and R.
R has shape (rank,npw), not the imposed adjoint of L. Therefore population
and arbitrary ordered non-Hermitian sources remain reconstructible even for
singular D. The redundant parametrization is deterministic, not amplitude noise.
No full PW-square matrix, saved-band closure or rank compression is used.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .prepared_reference import PreparedReference, array


@dataclass(frozen=True)
class FactorMean:
    frames: tuple[np.ndarray, ...]
    u: np.ndarray
    p: np.ndarray


@dataclass(frozen=True)
class FactorTangent:
    ket: tuple[np.ndarray, ...]
    bra: tuple[np.ndarray, ...]
    u: np.ndarray
    p: np.ndarray


@dataclass(frozen=True)
class MeanEvaluation:
    """Stage-local shared work. Do not mutate its state or arrays while reused."""
    model: object
    state: FactorMean
    active_density: np.ndarray
    active_gradient: np.ndarray
    delta_density: np.ndarray
    potential: np.ndarray
    h_frames: tuple[np.ndarray, ...]
    rhs: FactorMean
    acceleration: np.ndarray
    energy: complex


class PreparedFactorDynamics:
    """Actual RHS/energy/tangent consumer of the declared fixed scalar.

The mean uses per-cell densities and canonical mean u,p. Tangents can be
ordinary variations, or collective variables X=sqrt(N)*delta rho and the
matching canonical collective ionic u,p. Their linear equations coincide.
Mean and collective covariance normalizations must not be interchanged.
"""
    def __init__(self, reference: PreparedReference, sectors):
        self.reference = reference
        self.sectors = tuple(sectors)
        self.measure = reference.validate_sectors(self.sectors)

    def initial_mean(self, u, p):
        return self.mean(tuple(s.initial_frame for s in self.sectors), u, p)

    def mean(self, frames, u, p):
        if len(frames) != len(self.sectors):
            raise ValueError('one full-row frame per sector required')
        frames = tuple(array(x, s.initial_frame.shape, 'moving frame')
                       for s, x in zip(self.sectors, frames))
        return FactorMean(frames, array(u, (self.reference.ndof,), 'mean u', real=True),
                          array(p, (self.reference.ndof,), 'mean p', real=True))

    def tangent(self, ket, bra, u, p):
        if len(ket) != len(self.sectors) or len(bra) != len(self.sectors):
            raise ValueError('one ket/bra pair per sector required')
        return FactorTangent(
            tuple(array(x, s.initial_frame.shape, 'tangent ket') for s, x in zip(self.sectors, ket)),
            tuple(array(x, s.initial_frame.T.shape, 'independent tangent bra') for s, x in zip(self.sectors, bra)),
            array(u, (self.reference.ndof,), 'tangent u'),
            array(p, (self.reference.ndof,), 'tangent p'))

    def zero_tangent(self):
        return self.tangent(tuple(np.zeros_like(s.initial_frame) for s in self.sectors),
                            tuple(np.zeros_like(s.initial_frame.T) for s in self.sectors),
                            np.zeros(self.reference.ndof), np.zeros(self.reference.ndof))

    def _kernel(self, density):
        return array(self.reference.kernel(density), self.reference.nbar.shape, 'K density')

    def _h_action(self, index, vectors, evaluation):
        s = self.sectors[index]
        return (s.h0(vectors)
                + s.coupling.potential_action(evaluation.potential[:, s.density_component], vectors)
                + s.coupling.ionic_action(evaluation.state.u, vectors))

    def evaluate_mean(self, state):
        """Compute moments, K and H Psi ONCE for a stage and all its sources."""
        self.mean(state.frames, state.u, state.p)
        r = self.reference
        density = np.zeros(r.nbar.shape, complex)
        gradient = np.zeros(r.ndof, complex)
        h0_frames = []
        one_body_energy = 0j
        for s, psi in zip(self.sectors, state.frames):
            n, g = s.coupling.moments(psi, psi, s.d)
            density[:, s.density_component] += s.weight * n
            gradient += s.weight * g
            h0psi = s.h0(psi)
            h0_frames.append(h0psi)
            one_body_energy += s.weight * np.trace(s.d @ (psi.conj().T @ h0psi))
        delta = density + r.held_density - r.nbar
        potential = self._kernel(delta)
        hframes = tuple(h0psi + s.coupling.potential_action(potential[:, s.density_component], psi)
                        + s.coupling.ionic_action(state.u, psi)
                        for s, psi, h0psi in zip(self.sectors, state.frames, h0_frames))
        force = r.ionic_force + r.support_force - gradient - r.held_gradient - r.curvature @ state.u
        # Imaginary roundoff is retained in evaluated forces/energy; no physical
        # coefficient is Hermitized. A later real mean stepper owns its checks.
        rhs = FactorMean(tuple(-1j * h for h in hframes), state.p / r.masses, force)
        energy = (one_body_energy + self.measure * np.sum(delta * potential) / 2
                  + state.u @ (gradient + r.held_gradient - r.ionic_force - r.support_force)
                  + state.u @ r.curvature @ state.u / 2 + np.sum(state.p * state.p / r.masses) / 2)
        return MeanEvaluation(self, state, density, gradient, delta, potential,
                              hframes, rhs, force / r.masses, complex(energy))

    def tangent_moments(self, state, tangent):
        """Physical complex-linear density and positive gradient of delta rho."""
        r = self.reference
        density = np.zeros(r.nbar.shape, complex)
        gradient = np.zeros(r.ndof, complex)
        for s, psi, left, right in zip(self.sectors, state.frames, tangent.ket, tangent.bra):
            n1, g1 = s.coupling.moments(left, psi)
            n2, g2 = s.coupling.moments(psi, right.conj().T)
            density[:, s.density_component] += s.weight * (n1 + n2)
            gradient += s.weight * (g1 + g2)
        return density, gradient

    def tangent_rhs(self, evaluation, tangent):
        """Derivative of legs INCLUDING the moving-mean-frame connection.

Reconstructing the physical derivative also includes L Psi_dot^dagger and
Psi_dot R. This gives -i([H,delta rho]+[delta H,rho]) exactly.
A shared MeanEvaluation avoids recomputing mean H0/K/moments for each source.
"""
        if evaluation.model is not self:
            raise ValueError('mean evaluation belongs to a different model')
        self.tangent(tangent.ket, tangent.bra, tangent.u, tangent.p)
        density, gradient = self.tangent_moments(evaluation.state, tangent)
        field = self._kernel(density)
        ket, bra = [], []
        for index, (s, psi, left, right) in enumerate(zip(
                self.sectors, evaluation.state.frames, tangent.ket, tangent.bra)):
            drive = (s.coupling.potential_action(field[:, s.density_component], psi)
                     + s.coupling.ionic_action(tangent.u, psi))
            # Complex source: delta H is NOT Hermitian. Its independently
            # reconstructed adjoint uses conjugate field/direction, not drive.
            drive_adjoint = (s.coupling.potential_action(field[:, s.density_component].conj(), psi)
                             + s.coupling.ionic_action(tangent.u.conj(), psi))
            ket.append(-1j * (self._h_action(index, left, evaluation) + drive @ s.d))
            bra.append(1j * (self._h_action(index, right.conj().T, evaluation).conj().T
                             + s.d @ drive_adjoint.conj().T))
        return FactorTangent(tuple(ket), tuple(bra), tangent.p / self.reference.masses,
                             -self.reference.curvature @ tangent.u - gradient)

    def traces(self, state, tangent=None):
        """Raw per-sector traces; weights/number/spin remain explicit."""
        if tangent is None:
            return np.asarray([np.trace(s.d @ (psi.conj().T @ psi))
                               for s, psi in zip(self.sectors, state.frames)])
        return np.asarray([np.trace(psi.conj().T @ left) + np.trace(right @ psi)
                           for psi, left, right in zip(state.frames, tangent.ket, tangent.bra)])
