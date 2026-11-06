"""Unweighted full-PW Gamma dyad maps for prepared effective models.

X=left @ coefficient @ right.conj().T need not be Hermitian or a state.
Density is bohr^-3; moments return the POSITIVE ionic gradient in Ry/bohr.
The caller, not this channel, supplies spin/sector weights and force signs.
Only ordinary-NC, no-NLCC, full-complex k=q=Gamma inputs are supported; these
maps neither choose a preparation nor recalibrate PBE away from its expansion.
"""
from __future__ import annotations

import numpy as np

from .qe_material_gradient import GammaGradientReduction
from .qe_material_hamiltonian import PWFactors
from .qe_material_ionic_action import GammaIonicAction
from .qe_nc_projectors import centered_projectors, place_projectors


class GammaDyadCoupling:
    """Density/potential and all-coordinate gradient/action transpose pairs.

    Inputs remain the original bound factor/ionic objects and must not be
    mutated after construction. Identity checks include the placed projectors
    at every original atom; matching dimensions alone is insufficient. This
    one-time linear-size check does not cache duplicate placed projectors.
    All transient grid/PW factor arrays are bounded by block_size columns.
    """

    def __init__(self, factors: PWFactors, ionic: GammaIonicAction, block_size=8):
        if type(block_size) is not int or not 1 <= block_size <= 16:
            raise ValueError("block_size must be an integer in 1..16")
        self.factors = factors
        self.ionic = ionic
        self.block_size = block_size
        fft = factors.fft
        if (fft.physical != ionic.shape or fft.nnr != ionic.nnr
                or np.any(factors.k_bohr_inv != 0) or np.any(fft.q_bohr_inv != 0)
                or not np.array_equal(fft.wavevectors.T, ionic.p)
                or not np.array_equal(fft.nl, np.ravel_multi_index(
                    ionic.pw_indices, ionic.shape, order="F"))):
            raise ValueError("original Gamma PW rows, momenta and physical grid must match")
        self.npw = fft.miller.shape[1]
        self.nr = fft.nnr
        self.ndof = 3 * len(ionic.species)
        self.cell_volume_bohr3 = fft.omega
        self.grid_measure_bohr3 = fft.omega / fft.nnr
        covered = []
        names = set()
        for species in factors.species:
            if species.name in names or species.name not in ionic.tables:
                raise ValueError("one projector group per original species required")
            names.add(species.name)
            tables = ionic.tables[species.name]
            atoms = tuple(i + 1 for i, name in enumerate(ionic.species) if name == species.name)
            nh = tables.upf.nh
            if (species.atom_indices != atoms or not atoms
                    or tables.omega_bohr3 != fft.omega
                    or species.beta.shape != (self.npw, len(atoms) * nh)
                    or not np.array_equal(species.active_d, tables.upf.expanded_d)):
                raise ValueError("original species, atom order, volume and NC factors must match")
            if species.name not in ionic.centered:
                ionic.centered[species.name] = centered_projectors(tables, ionic.p)
            for slot, atom in enumerate(atoms):
                placed = place_projectors(ionic.centered[species.name], ionic.p,
                                          ionic.positions[:, atom - 1])
                if not np.array_equal(species.beta[:, slot*nh:(slot+1)*nh], placed):
                    raise ValueError("placed projectors disagree with original ionic geometry")
            covered.extend(atoms)
        if sorted(covered) != list(range(1, len(ionic.species) + 1)):
            raise ValueError("every original atom must be represented exactly once")
        self._gradient = GammaGradientReduction(factors, ionic)

    def _matrix(self, value, label):
        value = np.asarray(value)
        if value.ndim != 2 or value.shape[0] != self.npw or not np.isfinite(value).all():
            raise ValueError(f"{label} must be finite full-PW columns (npw,nvector)")
        return value

    def moments(self, left, right, coefficient=None):
        """Return (density[nr], positive_gradient[ndof]) for L C R†.

        None is the identity only when column counts agree; no identity matrix
        is allocated. Rectangular C is absorbed into blocks on the smaller leg
        without projecting the original PW rows or storing a PW-square dyad.
        """
        left, right = self._matrix(left, "left"), self._matrix(right, "right")
        nl, nr = left.shape[1], right.shape[1]
        if coefficient is None:
            if nl != nr:
                raise ValueError("identity coefficient requires equal column counts")
        else:
            coefficient = np.asarray(coefficient)
            if coefficient.shape != (nl, nr) or not np.isfinite(coefficient).all():
                raise ValueError("finite coefficient (nleft,nright) required")
        density = np.zeros(self.nr, dtype=np.complex128)
        gradient = np.zeros(self.ndof, dtype=np.complex128)
        fft = self.factors.fft
        for start in range(0, min(nl, nr), self.block_size):
            stop = min(start + self.block_size, min(nl, nr))
            if coefficient is None:
                lblock, rblock = left[:, start:stop], right[:, start:stop]
            elif nr <= nl:
                lblock, rblock = left @ coefficient[:, start:stop], right[:, start:stop]
            else:
                lblock = left[:, start:stop]
                rblock = right @ coefficient[start:stop, :].conj().T
            lwave, rwave = fft.from_g(lblock), fft.from_g(rblock)
            density += np.sum(lwave * rwave.conj(), axis=1) / self.cell_volume_bohr3
            gradient += self._gradient.nonlocal_block(lblock, rblock)
        gradient += self._gradient.local(density)
        return density, gradient

    def potential_action(self, field, vectors):
        """Apply a complex Ry field in the unchanged physical Fourier convention."""
        field = np.asarray(field)
        vectors = self._matrix(vectors, "vectors")
        if field.shape != (self.nr,) or not np.isfinite(field).all():
            raise ValueError("finite physical potential field (nr,) required")
        out = np.empty(vectors.shape, dtype=np.complex128)
        fft = self.factors.fft
        for start in range(0, vectors.shape[1], self.block_size):
            stop = min(start + self.block_size, vectors.shape[1])
            out[:, start:stop] = fft.to_g(field[:, None] * fft.from_g(vectors[:, start:stop]))
        return out

    def ionic_action(self, direction, vectors):
        """Apply positive g(direction), complex-linear in all atom-major entries."""
        direction = np.asarray(direction)
        vectors = self._matrix(vectors, "vectors")
        if direction.shape != (self.ndof,) or not np.isfinite(direction).all():
            raise ValueError("finite atom-major ionic direction (ndof,) required")
        prepared = self.ionic.prepare(direction.reshape(-1, 3))
        out = np.empty(vectors.shape, dtype=np.complex128)
        for start in range(0, vectors.shape[1], self.block_size):
            stop = min(start + self.block_size, vectors.shape[1])
            out[:, start:stop] = prepared.apply(vectors[:, start:stop])
        return out
