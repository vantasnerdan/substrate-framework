"""Checked same-spin edges of a declared complete Bloch grid.

The edge fold ell obeys kt=ks+q+B ell. It is NOT either source-star
reciprocal_shift. q is kept in the caller's representative, so reverse()
uses exactly -q and -ell rather than independently refolding its support.
"""
from dataclasses import dataclass

import numpy as np

from .qe_bloch_frame import BlochBasis
from .qe_bloch_grid import qe_fold, readonly


@dataclass(frozen=True)
class BlochEdge:
    source_basis: BlochBasis
    target_basis: BlochBasis
    q_fractional: np.ndarray
    q_cart_bohr_inv: np.ndarray
    reciprocal_fold: np.ndarray
    q_mesh_index: tuple[int, int, int]
    mesh: tuple[int, int, int]
    ispin: int

    @classmethod
    def create(cls, reference, source_basis, target_basis, q_fractional):
        """Bind two admitted full-PW bases, including generic commensurate q."""
        if not isinstance(source_basis, BlochBasis) or not isinstance(target_basis, BlochBasis):
            raise ValueError('two admitted BlochBasis objects required')
        source_basis.validate(reference)
        target_basis.validate(reference)
        source, target = source_basis.slot, target_basis.slot
        if source.ispin != target.ispin or source.ispin not in (1, 2):
            raise ValueError('directed scalar edge must preserve explicit spin')
        mesh = tuple(reference.monkhorst_pack[f'nk{i}'] for i in (1, 2, 3))
        if any(type(n) not in (int, np.int64, np.int32) or n < 1 for n in mesh):
            raise ValueError('positive integer declared complete-grid dimensions required')
        shift = np.array([reference.monkhorst_pack[f'k{i}'] for i in (1, 2, 3)])
        if not np.isin(shift, (0, 1)).all():
            raise ValueError('declared Monkhorst-Pack shifts must be zero or one')
        for slot in (source, target):
            index = np.asarray(slot.mesh_index)
            if (index.shape != (3,) or index.dtype.kind not in 'iu'
                    or np.any(index < 0) or np.any(index >= mesh)
                    or not np.allclose(slot.k_fractional, qe_fold((index+shift/2)/mesh),
                                       atol=2e-10, rtol=0)):
                raise ValueError('basis slot does not lie on declared complete grid')
        q = np.asarray(q_fractional, dtype=float)
        if q.shape != (3,) or not np.isfinite(q).all():
            raise ValueError('finite reciprocal-fractional q required')
        step = np.rint(q*np.asarray(mesh))
        if not np.allclose(q*np.asarray(mesh), step, atol=2e-10, rtol=0):
            raise ValueError('q must be commensurate with the declared complete grid')
        index = tuple(int(v) for v in step.astype(np.int64) % mesh)
        expected = (np.asarray(source.mesh_index)+np.asarray(index)) % mesh
        if not np.array_equal(expected, target.mesh_index):
            raise ValueError('target mesh slot is not source plus q')
        delta = target.k_fractional-source.k_fractional-q
        fold = np.rint(delta).astype(np.int64)
        if not np.allclose(delta, fold, atol=2e-10, rtol=0):
            raise ValueError('edge cannot be connected by an integer reciprocal fold')
        return cls(source_basis, target_basis, readonly(q),
                   readonly(reference.reciprocal_bohr_inv @ q), readonly(fold),
                   index, mesh, source.ispin)

    def validate(self, reference):
        expected = type(self).create(reference, self.source_basis, self.target_basis, self.q_fractional)
        if (self.mesh != expected.mesh or self.q_mesh_index != expected.q_mesh_index
                or self.ispin != expected.ispin
                or not np.array_equal(self.reciprocal_fold, expected.reciprocal_fold)
                or not np.array_equal(self.q_cart_bohr_inv, expected.q_cart_bohr_inv)):
            raise ValueError('edge identity/fold/units disagree with its checked bases')
        return self

    @property
    def source_dimension(self):
        return self.source_basis.miller.shape[1]

    @property
    def target_dimension(self):
        return self.target_basis.miller.shape[1]

    def reverse(self, reference):
        """Reverse geometry only; an operator adjoint must retain its own support."""
        return type(self).create(reference, self.target_basis, self.source_basis, -self.q_fractional)
