"""Matrix-free full-PW saved-density H0, with reusable ordinary-NC factors.

No dense PW-square Hamiltonian, dense projector-space D, diagonalization or
occupation model. Cached species factors share radial/angular work over atoms
and spin sectors with identical physical PW order. All energies are Ry.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_material_grid import PhysicalRhoFFT
from .qe_nc_projectors import centered_projectors, place_projectors


@dataclass(frozen=True)
class SpeciesProjectors:
    name: str
    atom_indices: tuple[int, ...]
    beta: np.ndarray             # npw x (nat_species * nh), atom-major
    active_d: np.ndarray         # nh x nh, shared; NOT global block diagonal

    def apply(self, block):
        # BLAS conjugate-transpose flag avoids copying the large cached beta.
        overlaps = zgemm(1., self.beta, block, trans_a=2).reshape(
            len(self.atom_indices), self.active_d.shape[0], block.shape[1])
        weighted = np.einsum("ij,ajb->aib", self.active_d, overlaps, optimize=True)
        return zgemm(1., self.beta, weighted.reshape(-1, block.shape[1]))


@dataclass(frozen=True)
class PWFactors:
    fft: PhysicalRhoFFT
    kinetic_ry: np.ndarray
    species: tuple[SpeciesProjectors, ...]
    k_bohr_inv: np.ndarray
    ik: int | None                   # None for reconstructed full-grid basis
    sources: tuple

    def matches(self, wfc):
        return (wfc.ik == self.ik and np.array_equal(wfc.miller, self.fft.miller)
                and np.array_equal(wfc.k_cart_bohr_inv, self.k_bohr_inv)
                and np.array_equal(wfc.reciprocal_bohr_inv, self.fft.reciprocal_bohr_inv))


def build_pw_factors(inputs, wfc) -> PWFactors:
    """Build once per distinct k/PW order; reuse across explicit spin actions."""
    ref = inputs.reference
    point, energies, _ = ref.select(wfc.ik, wfc.ispin)
    if (wfc.gamma_only or wfc.npol != 1 or wfc.scalef != 1 or point.npw != wfc.igwx
            or energies.size != wfc.nbnd
            or not np.allclose(wfc.k_cart_bohr_inv, point.k_cart_bohr_inv, atol=2e-13, rtol=0)
            or not np.allclose(wfc.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, atol=2e-13, rtol=0)):
        raise ValueError("WFC/reference identity, row count or physical units disagree")
    return _build_basis_factors(inputs, wfc.miller, wfc.k_cart_bohr_inv,
                               wfc.reciprocal_bohr_inv, wfc.ik, dict(wfc.source))


def build_pw_factors_for_basis(inputs, basis) -> PWFactors:
    """Checked reconstructed basis; never invent a saved XML ik record."""
    from .qe_bloch_frame import BlochBasis
    if not isinstance(basis, BlochBasis):
        raise ValueError("a checked BlochBasis is required")
    basis.validate(inputs.reference)
    provenance = dict(kind="unitary-reconstructed-basis", identity=basis.identity,
                      source=dict(basis.slot.source_binding),
                      xml_source=dict(basis.reference_source),
                      ispin=basis.slot.ispin, mesh_index=basis.slot.mesh_index,
                      source_ik=basis.slot.source_ik,
                      operation=basis.slot.operation.index,
                      gauge_policy=basis.slot.gauge_policy)
    return _build_basis_factors(inputs, basis.miller, basis.k_cart_bohr_inv,
                               basis.reciprocal_bohr_inv, None, provenance)


def _build_basis_factors(inputs, miller, k, reciprocal, ik, provenance):
    """Shared physical k+G algebra, independent of saved/reconstructed identity."""
    ref = inputs.reference
    volume = abs(float(np.linalg.det(ref.cell_bohr)))
    fft = PhysicalRhoFFT.create(ref.fft_smooth, miller, reciprocal, volume)
    momenta = (reciprocal @ miller + k[:, None]).T
    factors = []
    for species in ref.species:
        tables = inputs.tables[species.name]
        centered = centered_projectors(tables, momenta)
        atoms = tuple(i for i, name in enumerate(ref.atom_names) if name == species.name)
        beta = np.empty((miller.shape[1], len(atoms)*tables.upf.nh), dtype=np.complex128, order="F")
        for slot, atom in enumerate(atoms):
            beta[:, slot*tables.upf.nh:(slot+1)*tables.upf.nh] = place_projectors(
                centered, momenta, ref.positions_bohr[:, atom])
        factors.append(SpeciesProjectors(species.name, tuple(a+1 for a in atoms), beta, tables.upf.expanded_d))
    return PWFactors(fft, np.sum(momenta**2, axis=1), tuple(factors),
                     k, ik, inputs.sources + (provenance,))


@dataclass(frozen=True)
class MaterialHamiltonian:
    background: object
    factors: PWFactors
    ispin: int
    block_size: int = 8

    def __post_init__(self):
        if self.ispin not in (1, 2) or type(self.block_size) is not int or self.block_size < 1:
            raise ValueError("explicit one-based spin and positive block size required")
        if (self.factors.fft.physical != self.background.fft.physical
                or self.factors.sources[:-1] != self.background.inputs.sources):
            raise ValueError("H0 factors and background must share exactly bound saved inputs")

    def apply(self, vectors):
        """Return H0 vectors in unchanged PW row order; allocate only column blocks."""
        vectors = np.asarray(vectors)
        single = vectors.ndim == 1
        block_input = vectors[:, None] if single else vectors
        if (block_input.ndim != 2 or block_input.shape[0] != self.factors.kinetic_ry.size
                or not np.isfinite(block_input).all()):
            raise ValueError("finite vectors (npw,) or (npw,nvector) required")
        output = np.empty(block_input.shape, dtype=np.complex128)
        fft = self.factors.fft
        potential = self.background.vrs[:, self.ispin-1, None]
        for start in range(0, block_input.shape[1], self.block_size):
            stop = min(start+self.block_size, block_input.shape[1])
            block = block_input[:, start:stop]
            result = self.factors.kinetic_ry[:, None]*block
            result = result + fft.to_g(potential * fft.from_g(block))
            for species in self.factors.species:
                result += species.apply(block)
            output[:, start:stop] = result
        return output[:, 0] if single else output

    def saved_eigenpair_diagnostics(self, wfc):
        """Raw residuals are diagnostics of an output-density H0, not a gate."""
        if wfc.ispin != self.ispin or not self.factors.matches(wfc):
            raise ValueError("saved vectors must match this exact spin/k/PW order")
        _, energies, _ = self.background.inputs.reference.select(wfc.ik, wfc.ispin)
        action = self.apply(wfc.coefficients)
        residual = action - wfc.coefficients*energies[None, :]
        gram = wfc.coefficients.conj().T @ wfc.coefficients
        return dict(action_ry=action, residual_ry=residual, gram=gram,
                    projected_h0_ry=wfc.coefficients.conj().T @ action,
                    residual_norm_ry=np.linalg.norm(residual, axis=0),
                    eigenvalues_ry=energies, source=dict(wfc.source),
                    qualification="Saved eigenpairs precede output-density potential rebuild; no convergence gate or vector repair.")
