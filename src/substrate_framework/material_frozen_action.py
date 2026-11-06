"""Exact frozen same-k H0 + V + g(u), with independently sized FFT/GEMM panels.

Reuse checked material factors and prepared ionic derivatives; no radial DFT
reconstruction, saved-band closure, small-term cutoff or Hermitian repair.
The caller owns H0 work reservations and must keep shared factors immutable.
See P254 A28/fused-action-registration.json for the source-first contract.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg.blas import zgemm
from threadpoolctl import threadpool_limits

from .qe_material_hamiltonian import MaterialHamiltonian, SpeciesProjectors


class FrozenMaterialAction:
    """One immutable stage snapshot on an existing checked full-PW Hamiltonian.

    With neither potential nor direction this is ordinary H0 (including for
    energy contractions). Potential is a physical F-order field (nr,), direction
    is the original atom-major (ndof,) vector, both allowed to be complex.
    Coupling is required only for direction and must share h0.factors by identity.
    The returned full frame is the only full-frame allocation; transient arrays
    have at most projector_block_size PW columns or fft_block_size grid columns.
    Original factor/prepared-direction caches are additional sector-local data.
    """

    def __init__(self, h0, coupling=None, *, potential=None, direction=None,
                 fft_block_size=8, projector_block_size=64):
        if not isinstance(h0, MaterialHamiltonian):
            raise TypeError("checked MaterialHamiltonian required, not a callable wrapper")
        for name, value in (("fft_block_size", fft_block_size),
                            ("projector_block_size", projector_block_size)):
            if type(value) is not int or value < 1:
                raise ValueError(name + " must be a positive integer")
        if coupling is not None and coupling.factors is not h0.factors:
            raise ValueError("coupling and H0 must share the same checked PWFactors")
        self.h0 = h0
        self.fft_block_size = fft_block_size
        self.projector_block_size = projector_block_size
        fft = h0.factors.fft
        field = np.array(h0.background.vrs[:, h0.ispin-1], dtype=complex, copy=True)
        if field.shape != (fft.nnr,) or not np.isfinite(field).all():
            raise ValueError("finite original physical saved potential required")
        if potential is not None:
            potential = np.asarray(potential)
            if potential.shape != field.shape or not np.isfinite(potential).all():
                raise ValueError("finite physical stage field (nr,) required")
            field += potential
        self.direction = None
        self.terms = ()
        if direction is not None:
            if coupling is None:
                raise ValueError("checked ionic coupling required for direction")
            u = np.array(direction, dtype=complex, copy=True)
            if u.shape != (coupling.ndof,) or not np.isfinite(u).all():
                raise ValueError("finite original atom-major direction (ndof,) required")
            u.setflags(write=False)
            self.direction = u
            # This is an exact zero test, not a small-displacement cutoff.
            if np.any(u != 0):
                with threadpool_limits(limits=1):
                    prepared = coupling.ionic.prepare(u.reshape(-1, 3))
                ionic_field = prepared.field
                # GammaIonicAction exposes a physical 3D grid; Bloch exposes its
                # F-order flattening. The checked coupling fixes this identity.
                if ionic_field.shape == fft.physical:
                    ionic_field = ionic_field.reshape(-1, order="F")
                if ionic_field.shape != field.shape:
                    raise ValueError("ionic and Hamiltonian physical grids differ")
                field += ionic_field
                self.terms = prepared.terms
        field.setflags(write=False)
        self.field = field
        self.has_local = bool(np.any(field != 0))
        # Only small species D matrices are conjugated/copied, never full beta.
        self.adjoint_species = tuple(
            SpeciesProjectors(s.name, s.atom_indices, s.beta, s.active_d.conj().T)
            for s in h0.factors.species)

    def apply(self, vectors):
        """Apply the frozen operator; no implicit work-count reservation."""
        with threadpool_limits(limits=1):
            return self._apply(vectors, adjoint=False)

    def adjoint_apply(self, vectors):
        """Exact Hilbert adjoint, also for complex fields/directions/raw factors."""
        with threadpool_limits(limits=1):
            return self._apply(vectors, adjoint=True)

    def _apply(self, vectors, *, adjoint):
        vectors = np.asarray(vectors)
        single = vectors.ndim == 1
        x = vectors[:, None] if single else vectors
        factors = self.h0.factors
        if x.ndim != 2 or x.shape[0] != factors.kinetic_ry.size:
            raise ValueError("finite full-PW vectors (npw,) or (npw,ncolumn) required")
        output = np.empty(x.shape, dtype=complex)
        kinetic = factors.kinetic_ry.conj() if adjoint else factors.kinetic_ry
        field = self.field.conj() if adjoint else self.field
        fft = factors.fft
        for start in range(0, x.shape[1], self.fft_block_size):
            sl = slice(start, start+self.fft_block_size)
            block = x[:, sl]
            if not np.isfinite(block).all():
                raise ValueError("finite full-PW vectors required")
            output[:, sl] = kinetic[:, None]*block
            if self.has_local:
                grid = fft.from_g(block)
                grid *= field[:, None]
                output[:, sl] += fft.to_g(grid)
        species = self.adjoint_species if adjoint else factors.species
        # Wide species overlaps are deliberately OUTSIDE FFT column blocks.
        # Reuse the admitted SpeciesProjectors contraction, including raw D.
        for start in range(0, x.shape[1], self.projector_block_size):
            sl = slice(start, start+self.projector_block_size)
            block = x[:, sl]
            for item in species:
                output[:, sl] += item.apply(block)
            for pu, bd, bra in self.terms:
                # Prepared DirectionalAction's two independent derivative legs.
                # BLAS transpose flags implement the true adjoint without copying
                # PW-sized factors or assuming that active D is Hermitian.
                if adjoint:
                    p = pu.conj()
                    overlap = zgemm(1., bd, p[:, None]*block, trans_a=2)
                    output[:, sl] += zgemm(1j, bra, overlap, trans_a=2)
                    overlap = zgemm(1., bd, block, trans_a=2)
                    leg = zgemm(-1j, bra, overlap, trans_a=2)
                    leg *= p[:, None]
                else:
                    overlap = zgemm(1., bra, pu[:, None]*block)
                    output[:, sl] += zgemm(1j, bd, overlap)
                    overlap = zgemm(1., bra, block)
                    leg = zgemm(-1j, bd, overlap)
                    leg *= pu[:, None]
                output[:, sl] += leg
        return output[:, 0] if single else output
