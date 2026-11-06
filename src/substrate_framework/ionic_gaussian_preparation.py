"""Conditional local-trap vacuum at t=0, never a repair of dynamical C(q).

Periodic real-space canonical coordinates use the existing ordered IonicPreparation
API. Orthonormal q coordinates pair q with -q without a dagger; per-cell mean
coordinates have CCR i/N and therefore are deliberately NOT IonicPreparation.
The Hamiltonian and mathematical checks are preregistered in A28's
ionic-gaussian-registration.json. This is not a measured cooling protocol.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .inertial_fermion_dynamics import ionic_masses_ry
from .material_modes import AMU_RY
from .material_force_lattice import validated_force_grid
from .physical_bilinear_preparation import IonicPreparation, CrossPreparation


class MissingTrapPreparation(ValueError):
    """Actual clamped stiffness does not define the declared positive trap."""


@dataclass(frozen=True)
class IonicGaussianPreparation:
    """One identical independent cell vacuum, with exact periodic replication."""
    local: IonicPreparation
    grid_shape: tuple[int, int, int]
    masses_ry: np.ndarray
    stiffness_raw: np.ndarray
    omega_ry: np.ndarray

    @property
    def ncells(self):
        return int(np.prod(self.grid_shape))

    @property
    def ndof(self):
        return len(self.masses_ry)

    @property
    def symmetric(self):
        return (self.local.ordered + self.local.ordered.T)/2

    @property
    def mean_ordered(self):
        """Covariance of bar_u and bar_p, with CCR i/N (not canonical i)."""
        return self.local.ordered/self.ncells

    def _index(self, index):
        index = np.asarray(index)
        if index.shape != (3,) or not np.issubdtype(index.dtype, np.integer):
            raise ValueError('three integer periodic cell/q indices required')
        return index % np.asarray(self.grid_shape)

    def real_pair(self, cell, other):
        """<(u,p)_cell (u,p)_other>, no conjugation of the second operator."""
        if np.array_equal(self._index(cell), self._index(other)):
            return self.local.ordered.copy()
        return np.zeros_like(self.local.ordered)

    def q_pair(self, q, other, *, dagger=False):
        """Orthonormal FFT: delta(q+q') W, or delta(q-q') W with dagger."""
        left, right = self._index(q), self._index(other)
        difference = left-right if dagger else left+right
        if np.all(difference % np.asarray(self.grid_shape) == 0):
            return self.local.ordered.copy()
        return np.zeros_like(self.local.ordered)

    def product_cross(self, channels):
        """Reuse the existing explicit zero cross API, not independent noise."""
        return CrossPreparation.product(channels, self.ndof,
            self.local.state_authority+'; electronic/ionic tensor product at t=0')

    def archive_arrays(self):
        """Factorized exact 64-cell covariance; no dense 14208-square allocation.

        q_minus_index and the local block specify every q,q' covariance exactly.
        In real space the omitted different-cell blocks are identically zero.
        """
        indices = np.indices(self.grid_shape).reshape(3, -1).T
        w, v, d = self.local.ordered, self.symmetric, self.ndof
        return dict(masses_ry=self.masses_ry, stiffness_raw_ry_bohr2=self.stiffness_raw,
            stiffness_preparation_ry_bohr2=self.stiffness_raw.real,
            omega_ry=self.omega_ry, mean_u_bohr=np.zeros(d), mean_p_canonical=np.zeros(d),
            local_ordered=w, local_symmetrized=v, local_commutator=w-w.T,
            ordered_uu=w[:d,:d], ordered_pp=w[d:,d:], ordered_up=w[:d,d:], ordered_pu=w[d:,:d],
            symmetrized_uu=v[:d,:d], symmetrized_pp=v[d:,d:], symmetrized_up=v[:d,d:],
            mean_ordered=self.mean_ordered, mean_symmetrized=v/self.ncells,
            mean_commutator=(w-w.T)/self.ncells, q_index=indices,
            q_minus_index=(-indices) % np.asarray(self.grid_shape))


def prepare_local_trap_vacuum(force_lattice, atomic_masses_amu, *, state_authority,
                             reality_relative_tolerance=1e-10):
    """Ground states of positive real R=0 same-coordinate clamped stiffness.

    Only the preparation coefficient's insignificant imaginary residual may be
    omitted, under the stated relative criterion. Raw complex data survive.
    No C(q) is diagonalized, clipped, symmetrized or modified, and no electronic
    spin/k measure is applied. All cells and Cartesian traps are independent.
    """
    force = validated_force_grid(force_lattice)
    nat = force.shape[3]
    mass = np.asarray(atomic_masses_amu)
    if mass.shape != (nat,) or np.iscomplexobj(mass):
        raise ValueError('one original real atomic mass per atom required')
    if not np.isfinite(reality_relative_tolerance) or not 0 <= reality_relative_tolerance < 1:
        raise ValueError('explicit nonnegative relative reality tolerance below one required')
    masses = ionic_masses_ry(mass, np.arange(1, nat+1), AMU_RY)
    raw = np.diag(force[0,0,0].reshape(3*nat, 3*nat)).copy()
    bad = (raw.real <= 0) | (np.abs(raw.imag) > reality_relative_tolerance*np.abs(raw.real))
    if np.any(bad):
        details = [(int(a//3), int(a%3), float(raw[a].real), float(raw[a].imag))
                   for a in np.flatnonzero(bad)]
        raise MissingTrapPreparation('missing positive real clamped-site preparation stiffness '
                                     '(zero-based atom,xyz,real,imag): '+str(details))
    omega = np.sqrt(raw.real/masses)
    widths = masses*omega
    n = len(masses)
    ordered = np.zeros((2*n, 2*n), complex)
    index = np.arange(n)
    ordered[index,index] = 1/(2*widths)
    ordered[n+index,n+index] = widths/2
    ordered[index,n+index] = .5j
    ordered[n+index,index] = -.5j
    local = IonicPreparation(ordered, state_authority)
    for value in (raw, masses, omega):
        value.flags.writeable = False
    return IonicGaussianPreparation(local, tuple(force.shape[:3]), masses, raw, omega)
