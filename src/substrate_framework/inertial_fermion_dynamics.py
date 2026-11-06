"""Finite held-reference electron/ion tangent, in Ry time (hbar=1).

The action has physical +g*u, independent supplied Cfrozen, and 1/2 n.K.n.
State order is (u, p, original-basis delta_rho for each explicit spin channel).
No numerical eta, fitted counterterm, population relaxation or material mapping.

On a unitary orbit rho=exp(-i*S) rho0 exp(i*S), its second-order
one-body energy is -Tr(H0 [S,[S,rho0]])/2. Entropy and occupation
constraints are constant on this orbit. Together with ionic kinetic/frozen,
u.gradient and density-interaction energies this is the conserved quadratic
energy of the tangent, not an integrated-work definition. Equal-population
blocks are outside that tangent; they are retained by the full EOM, but the
orbit-energy API rejects their nonzero initial perturbations. A zero energy
gap with unequal populations is allowed (possibly indefinite/zero energy).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm, solve


def _array(value, shape, name):
    value = np.asarray(value, dtype=complex)
    if value.shape != shape or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return value


def _hermitian(value, name):
    if not np.allclose(value, value.swapaxes(-1, -2).conj(), atol=1e-12, rtol=1e-12):
        raise ValueError(f"{name} must be Hermitian")


def _real_symmetric(value, name):
    if np.max(np.abs(value.imag), initial=0) > 1e-12:
        raise ValueError(f"{name} must be real")
    _hermitian(value, name)
    return value.real


def commutator(a, b):
    return a @ b - b @ a


def ionic_masses_ry(amass_species, ityp_one_based, amu_ry, cartesian=3):
    """QE M=amass(ityp)*AMU_RY, AMU_RY=AMU_SI/(2*m_e); no extra 2."""
    species = np.asarray(amass_species, dtype=float)
    indices = np.asarray(ityp_one_based)
    if (species.ndim != 1 or not len(species) or not np.all(np.isfinite(species))
            or np.any(species <= 0) or indices.ndim != 1 or not len(indices)
            or not np.issubdtype(indices.dtype, np.integer)
            or np.any(indices < 1) or np.any(indices > len(species))
            or not np.isfinite(amu_ry) or amu_ry <= 0
            or not isinstance(cartesian, int) or cartesian <= 0):
        raise ValueError("positive species masses, valid one-based ityp and AMU_RY required")
    return np.repeat(species[indices - 1] * amu_ry, cartesian)


@dataclass(frozen=True)
class FermionChannel:
    """One separately conserved spin sector, with explicit observable weight.

    Occupations are RAW probabilities, not quadrature-weighted. Columns of
    orbitals give their original-basis frame. Exact equal occupations identify
    null orbit sectors; no occupation cutoff is applied. H0 must commute with
    rho0. Density operators already identify their global density/spin slots.
    """
    h0: np.ndarray
    occupations: np.ndarray
    orbitals: np.ndarray
    g: np.ndarray
    densities: np.ndarray
    weight: float = 1.0
    spin: float = 0.0

    def __post_init__(self):
        f = np.asarray(self.occupations, dtype=float)
        if (f.ndim != 1 or not len(f) or not np.all(np.isfinite(f))
                or np.any(f < 0) or np.any(f > 1)):
            raise ValueError("occupations must be finite probabilities")
        d = len(f)
        h = _array(self.h0, (d, d), "H0")
        v = _array(self.orbitals, (d, d), "orbitals")
        g = np.asarray(self.g, dtype=complex)
        density = np.asarray(self.densities, dtype=complex)
        if g.ndim != 3 or density.ndim != 3:
            raise ValueError("g and densities must be arrays of operators")
        g = _array(g, (len(g), d, d), "g")
        density = _array(density, (len(density), d, d), "densities")
        for value, name in ((h, "H0"), (g, "g"), (density, "densities")):
            _hermitian(value, name)
        if not np.allclose(v.conj().T @ v, np.eye(d), atol=1e-12, rtol=1e-12):
            raise ValueError("orbitals must be unitary")
        if not np.isfinite(self.weight) or self.weight <= 0 or not np.isfinite(self.spin):
            raise ValueError("positive finite weight and finite spin required")
        rho = (v * f) @ v.conj().T
        if np.linalg.norm(commutator(h, rho)) > 1e-12 * max(1, np.linalg.norm(h)):
            raise ValueError("reference must be stationary under H0")
        for name, value in (("h0", h), ("occupations", f), ("orbitals", v),
                            ("g", g), ("densities", density)):
            value = value.copy()
            value.flags.writeable = False
            object.__setattr__(self, name, value)

    @property
    def rho0(self):
        return (self.orbitals * self.occupations) @ self.orbitals.conj().T


class InertialFermionDynamics:
    """Original-basis complex-linear EOM; physical states are Hermitian/real."""

    def __init__(self, channels, interaction, frozen_curvature, masses):
        self.channels = tuple(channels)
        if not self.channels:
            raise ValueError("at least one explicit channel required")
        self.nion = len(self.channels[0].g)
        self.ndensity = len(self.channels[0].densities)
        if not self.nion or not self.ndensity:
            raise ValueError("nonempty ionic and density spaces required")
        if any(len(c.g) != self.nion or len(c.densities) != self.ndensity
               for c in self.channels):
            raise ValueError("channel ionic/density dimensions must agree")
        self.interaction = _real_symmetric(_array(
            interaction, (self.ndensity, self.ndensity), "K"), "K").copy()
        self.frozen_curvature = _real_symmetric(_array(
            frozen_curvature, (self.nion, self.nion), "Cfrozen"), "Cfrozen").copy()
        mass = np.asarray(masses, dtype=float)
        if mass.shape != (self.nion,) or not np.all(np.isfinite(mass)) or np.any(mass <= 0):
            raise ValueError("one positive Ry mass per ionic coordinate required")
        self.masses = mass.copy()
        self.rho_slices = []
        cursor = 2 * self.nion
        for channel in self.channels:
            size = len(channel.occupations) ** 2
            self.rho_slices.append(slice(cursor, cursor + size))
            cursor += size
        self.state_size = cursor
        self._references = tuple(c.rho0 for c in self.channels)

    def zero_state(self):
        return np.zeros(self.state_size, dtype=complex)

    def unpack(self, state):
        state = _array(state, (self.state_size,), "state")
        matrices = tuple(state[s].reshape(c.h0.shape)
                         for s, c in zip(self.rho_slices, self.channels))
        return state[:self.nion], state[self.nion:2*self.nion], matrices

    def pack(self, displacement, momentum, matrices):
        if len(matrices) != len(self.channels):
            raise ValueError("one density matrix per channel required")
        parts = [_array(displacement, (self.nion,), "u"),
                 _array(momentum, (self.nion,), "p")]
        parts.extend(_array(rho, c.h0.shape, "delta_rho").ravel()
                     for rho, c in zip(matrices, self.channels))
        return np.concatenate(parts)

    def observables(self, state):
        """Return physical density and +ionic gradient (negative force)."""
        _, _, matrices = self.unpack(state)
        density = np.zeros(self.ndensity, dtype=complex)
        gradient = np.zeros(self.nion, dtype=complex)
        for channel, rho in zip(self.channels, matrices):
            density += channel.weight * np.einsum("aij,ji->a", channel.densities, rho)
            gradient += channel.weight * np.einsum("aij,ji->a", channel.g, rho)
        return density, gradient

    def conserved_populations(self, state):
        """Per-channel traces, weighted total particle number and spin."""
        _, _, matrices = self.unpack(state)
        traces = np.asarray([np.trace(rho) for rho in matrices])
        weights = np.asarray([c.weight for c in self.channels])
        spins = np.asarray([c.spin for c in self.channels])
        return traces, weights @ traces, (weights * spins) @ traces

    def derivative(self, state, external_force=None):
        u, p, matrices = self.unpack(state)
        density, gradient = self.observables(state)
        potential = self.interaction @ density
        electronic = []
        for c, rho0, delta in zip(self.channels, self._references, matrices):
            perturbation = (np.einsum("a,aij->ij", u, c.g)
                            + np.einsum("s,sij->ij", potential, c.densities))
            electronic.append(-1j * (commutator(c.h0, delta)
                                    + commutator(perturbation, rho0)))
        force = -self.frozen_curvature @ u - gradient
        if external_force is not None:
            force = force + _array(external_force, (self.nion,), "external_force")
        return self.pack(p / self.masses, force, electronic)

    def constant_generator(self):
        """Dense finite generator; not intended for full plane-wave production."""
        basis = np.eye(self.state_size, dtype=complex)
        return np.column_stack([self.derivative(column) for column in basis.T])

    def propagator(self, tau):
        if not np.isfinite(tau):
            raise ValueError("finite time required")
        return expm(tau * self.constant_generator())

    def propagate(self, state, tau):
        self.unpack(state)
        return self.propagator(tau) @ state

    def driven_state(self, z, force):
        """exp(-i*z*tau) response with zero independent initial electronic source."""
        if not np.isfinite(z):
            raise ValueError("finite frequency required")
        source = self.zero_state()
        source[self.nion:2*self.nion] = _array(force, (self.nion,), "force")
        return solve(-1j*z*np.eye(self.state_size) - self.constant_generator(), source)

    def electronic_curvature(self, z):
        generator = self.constant_generator()
        start = 2*self.nion
        response = solve(-1j*z*np.eye(self.state_size-start) - generator[start:, start:],
                         generator[start:, :self.nion])
        return -generator[self.nion:start, start:] @ response

    def frequency_inverse(self, z):
        return self.frozen_curvature + self.electronic_curvature(z) - z*z*np.diag(self.masses)

    def energy_terms(self, state, tolerance=1e-10):
        """Conserved source-action quadratic on the real unitary-orbit tangent.

        Non-passive occupations and unstable supplied curvature are allowed;
        conservation does NOT assert positive definiteness or equilibrium.
        Equal-population blocks (including both-zero/one blocks and populations)
        have no orbit coordinate. Reject, rather than pseudo-invert, their data.
        """
        u, p, matrices = self.unpack(state)
        if np.max(np.abs(np.r_[u.imag, p.imag]), initial=0) > tolerance:
            raise ValueError("orbit energy requires real ionic coordinates")
        bare = 0.0
        for c, rho0, delta in zip(self.channels, self._references, matrices):
            _hermitian(delta, "physical delta_rho")
            rotated = c.orbitals.conj().T @ delta @ c.orbitals
            difference = c.occupations[None, :] - c.occupations[:, None]
            if np.linalg.norm(rotated[difference == 0]) > tolerance:
                raise ValueError("equal-population perturbation is not a unitary-orbit tangent")
            s = np.divide(1j*rotated, difference, out=np.zeros_like(rotated), where=difference != 0)
            s = c.orbitals @ s @ c.orbitals.conj().T
            bare += (-0.5*c.weight*np.trace(c.h0 @ commutator(s, commutator(s, rho0)))).real
        density, gradient = self.observables(state)
        terms = {
            "ionic_kinetic": float(0.5*np.sum(p.real**2 / self.masses)),
            "ionic_frozen": float(0.5*u.real @ self.frozen_curvature @ u.real),
            "electronic_orbit": float(bare),
            "electron_ion": float((u @ gradient).real),
            "density_interaction": float((0.5*density @ self.interaction @ density).real),
        }
        terms["total"] = sum(terms.values())
        return terms

    def tangent_basis(self):
        """Real coordinates: ionic u,p then Hermitian X,Y for every f_i != f_j."""
        columns = []
        for index in range(2*self.nion):
            state = self.zero_state()
            state[index] = 1
            columns.append(state)
        for c, section in zip(self.channels, self.rho_slices):
            dimension = len(c.occupations)
            for i in range(dimension):
                for j in range(i+1, dimension):
                    if c.occupations[i] == c.occupations[j]:
                        continue
                    for phase in (1, 1j):
                        matrix = np.zeros_like(c.h0)
                        matrix[i, j], matrix[j, i] = phase/np.sqrt(2), np.conj(phase)/np.sqrt(2)
                        state = self.zero_state()
                        state[section] = (c.orbitals @ matrix @ c.orbitals.conj().T).ravel()
                        columns.append(state)
        return np.column_stack(columns)

    def prepared_covariance(self, coordinate_covariance):
        """Classical Gaussian on real tangent coordinates, NOT coupled quantum Gibbs.

        Caller supplies the entire real covariance (including cross terms).
        No f(1-f), thermal FDT, spin degeneracy or independent noise is invented.
        Equal-population/zero-population null-orbit blocks are absent by design;
        quantum quasifree occupation/degenerate noise is NOT represented.
        """
        basis = self.tangent_basis()
        covariance = _real_symmetric(_array(coordinate_covariance,
            (basis.shape[1], basis.shape[1]), "coordinate covariance"), "coordinate covariance")
        if np.linalg.eigvalsh(covariance).min() < -1e-12:
            raise ValueError("coordinate covariance must be positive semidefinite")
        return basis @ covariance @ basis.conj().T

    def propagate_covariance(self, covariance, tau):
        covariance = _array(covariance, (self.state_size, self.state_size), "covariance")
        _hermitian(covariance, "covariance")
        if np.linalg.eigvalsh(covariance).min() < -1e-12:
            raise ValueError("covariance must be positive semidefinite")
        propagator = self.propagator(tau)
        return propagator @ covariance @ propagator.conj().T
