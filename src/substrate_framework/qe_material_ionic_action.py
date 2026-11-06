"""Full-complex Gamma bare NC ionic derivative, in Ry/bohr.

QE 7.5 PHonon/PH/compute_dvloc.f90:91-106 and dvqpsi_us_only.f90:
135-141,199-238 fix the local phase and BOTH projector derivatives.
The saved physical reciprocal basis is already bohr^-1 (pw_restart_new).
For dimensionless u[nat,3], differentiate tau+epsilon*u, epsilon in bohr.
Complex u is extended linearly, never conjugated. No screening, H0, or
saved-band closure is implied. Only requested fields/active atoms are cached.
"""
from dataclasses import dataclass
import numpy as np
from .qe_nc_projectors import centered_projectors, place_projectors


def fft_indices(miller, shape):
    miller = np.asarray(miller)
    if miller.ndim != 2 or miller.shape[0] != 3 or miller.dtype.kind not in 'iu':
        raise ValueError('integer Miller columns [3,n] required')
    indices = tuple(miller[i] % shape[i] for i in range(3))
    if len(np.unique(np.ravel_multi_index(indices, shape))) != miller.shape[1]:
        raise ValueError('physical FFT aliases distinct Miller rows')
    return indices


@dataclass(frozen=True)
class DirectionalAction:
    """One requested full nat*3 direction; field Ry/bohr, factors physical."""
    owner: object
    direction: np.ndarray
    field: np.ndarray
    terms: tuple

    def components(self, psi):
        """Return local, derivative-of-bra, derivative-of-ket, each [npw,n]."""
        psi = np.asarray(psi, dtype=np.complex128)
        if psi.ndim != 2 or psi.shape[0] != len(self.owner.p) or not 1 <= psi.shape[1] <= 16:
            raise ValueError('full PW columns in blocks of 1..16 required')
        local = self.owner.local_action(self.field, psi)
        bra, ket = np.zeros_like(local), np.zeros_like(local)
        for pu, bd, adjoint in self.terms:
            bra += 1j * (bd @ (adjoint @ (pu[:, None] * psi)))
            ket -= 1j * pu[:, None] * (bd @ (adjoint @ psi))
        return local, bra, ket

    def apply(self, psi):
        local, bra, ket = self.components(psi)
        local += bra
        local += ket
        return local


class GammaIonicAction:
    """Reuse species tables, centered B, placed Bdagger/BD, and physical FFT.

    Inputs: physical fft_shape[3], integer pw_miller[3,npw] and
    local_miller[3,ngm], reciprocal_bohr_inv[3,3] columns,
    positions_bohr[3,nat], original atom species sequence, and species->
    NCRadialTables. Caller explicitly establishes k=q=0 and full-complex PW.
    Local Fourier support is the saved density cutoff, not every FFT bin.
    """
    def __init__(self, fft_shape, pw_miller, local_miller, reciprocal_bohr_inv,
                 positions_bohr, species, tables):
        self.shape = tuple(int(x) for x in fft_shape)
        if len(self.shape) != 3 or min(self.shape) < 1:
            raise ValueError('positive physical FFT dimensions required')
        self.pw_indices = fft_indices(pw_miller, self.shape)
        self.local_indices = fft_indices(local_miller, self.shape)
        reciprocal = np.asarray(reciprocal_bohr_inv, dtype=float)
        self.p = (reciprocal @ pw_miller).T
        self.g = (reciprocal @ local_miller).T
        self.positions = np.asarray(positions_bohr, dtype=float)
        self.species = tuple(species)
        if self.positions.shape != (3, len(self.species)) or not np.isfinite(self.positions).all():
            raise ValueError('original full atom geometry required')
        self.tables = tables
        self.centered = {}
        self.atoms = {}
        self.vloc = {}
        self.nnr = int(np.prod(self.shape))

    def atom_factors(self, atom):
        if atom not in self.atoms:
            symbol = self.species[atom]
            if symbol not in self.centered:
                self.centered[symbol] = centered_projectors(self.tables[symbol], self.p)
            b = place_projectors(self.centered[symbol], self.p, self.positions[:, atom])
            self.atoms[atom] = (b @ self.tables[symbol].upf.expanded_d, b.conj().T.copy())
        return self.atoms[atom]

    def local_action(self, field, psi):
        """Unitary PW/grid transforms; 1/sqrt(volume) cancels in V*psi."""
        out = np.empty_like(psi, dtype=np.complex128)
        grid = np.zeros(self.shape, dtype=np.complex128)
        for band in range(psi.shape[1]):
            grid.fill(0)
            grid[self.pw_indices] = psi[:, band]
            wave = np.fft.ifftn(grid, norm='ortho')
            wave *= field
            out[:, band] = np.fft.fftn(wave, norm='ortho')[self.pw_indices]
        return out

    def prepare(self, direction):
        u = np.array(direction, dtype=np.complex128, copy=True)
        if u.shape != (len(self.species), 3) or not np.isfinite(u).all():
            raise ValueError('finite dimensionless direction [nat,3] required')
        coefficients = np.zeros(len(self.g), dtype=np.complex128)
        terms = []
        for atom in np.flatnonzero(np.any(u != 0, axis=1)):
            symbol = self.species[atom]
            if symbol not in self.vloc:
                self.vloc[symbol] = self.tables[symbol].vloc(np.linalg.norm(self.g, axis=1))
            coefficients += (-1j * (self.g @ u[atom]) *
                             np.exp(-1j * (self.g @ self.positions[:, atom])) * self.vloc[symbol])
            bd, adjoint = self.atom_factors(int(atom))
            terms.append((self.p @ u[atom], bd, adjoint))
        grid = np.zeros(self.shape, dtype=np.complex128)
        grid[self.local_indices] = coefficients
        field = np.fft.ifftn(grid) * self.nnr
        u.setflags(write=False)
        field.setflags(write=False)
        return DirectionalAction(self, u, field, tuple(terms))


def saved_band_overlap(saved_adjoint, full_pw_action):
    """Explicit saved-band bilinear, NOT a full-space response Hamiltonian."""
    return saved_adjoint @ full_pw_action
