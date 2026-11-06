"""Bound WFC admission and inherited-order unitary PW frame transport.

All source rows and raw band columns survive. Cutoff boundaries are recorded,
not repaired: strict-interior completeness is required, boundary membership
is inherited from the original native file. No native target ordering claim.
"""
from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType
import numpy as np
from .qe_bloch_grid import readonly
from .qe_material_grid import PhysicalRhoFFT
from .qe_saved_native import load_wavefunctions

CUTOFF_BOUNDARY_RY = 1e-8


def _digest(*arrays):
    h = sha256()
    for a in arrays:
        a = np.asarray(a)
        h.update(str((a.shape, a.dtype.str)).encode())
        h.update(a.tobytes(order='C'))
    return h.hexdigest()


def _cutoff_admission(reference, source):
    """Enumerate small integer bounding box, never a PW-square object."""
    b, k, m = source.reciprocal_bohr_inv, source.k_cart_bohr_inv, source.miller
    energy = np.sum((b @ m + k[:, None])**2, axis=0)
    cutoff = 2*reference.ecutwfc_ha
    if np.any(energy > cutoff + CUTOFF_BOUNDARY_RY):
        raise ValueError('original source rows exceed declared PW cutoff')
    rows = {tuple(row) for row in m.T}
    if len(rows) != source.igwx:
        raise ValueError('duplicate original PW rows')
    inverse = np.linalg.inv(b)
    center = -(inverse @ k)
    radius = np.linalg.norm(inverse, axis=1)*np.sqrt(cutoff+CUTOFF_BOUNDARY_RY)
    low, high = np.ceil(center-radius).astype(int), np.floor(center+radius).astype(int)
    missing_boundary = 0
    for i in range(low[0], high[0]+1):
        for j in range(low[1], high[1]+1):
            z = np.arange(low[2], high[2]+1)
            candidates = np.array([np.full(z.size, i), np.full(z.size, j), z])
            e = np.sum((b @ candidates + k[:, None])**2, axis=0)
            for index in np.flatnonzero(e <= cutoff+CUTOFF_BOUNDARY_RY):
                if tuple(candidates[:, index]) not in rows:
                    if e[index] < cutoff-CUTOFF_BOUNDARY_RY:
                        raise ValueError('source file omits a strict-interior PW row')
                    missing_boundary += 1
    return dict(cutoff_ry=cutoff, boundary_tolerance_ry=CUTOFF_BOUNDARY_RY,
                included_boundary_rows=int(np.count_nonzero(np.abs(energy-cutoff) <= CUTOFF_BOUNDARY_RY)),
                omitted_boundary_rows=missing_boundary,
                row_policy='complete strict interior; original native cutoff boundary retained')


def admit_source(reference, slot, source):
    """Match each actual spin header, all rows and bands independently."""
    point, energies, _ = reference.select(slot.source_ik, slot.ispin)
    if (source.ik != slot.source_ik or source.ispin != slot.ispin
            or dict(source.source) != dict(slot.source_binding)
            or source.gamma_only or source.npol != 1 or source.scalef != 1
            or source.igwx != point.npw or source.nbnd != energies.size
            or source.coefficients.shape != (source.igwx, source.nbnd)
            or not np.allclose(source.k_cart_bohr_inv, point.k_cart_bohr_inv, atol=2e-13, rtol=0)
            or not np.allclose(source.reciprocal_bohr_inv, reference.reciprocal_bohr_inv, atol=2e-13, rtol=0)):
        raise ValueError('original WFC spin/k/header/binding does not match XML sector')
    PhysicalRhoFFT.create(reference.fft_smooth, source.miller, source.reciprocal_bohr_inv,
                          abs(np.linalg.det(reference.cell_bohr)))
    audit = _cutoff_admission(reference, source)
    audit.update(source_sha256=source.source['sha256'], ispin=source.ispin, ik=source.ik,
                 source_rows_sha256=_digest(source.miller))
    return MappingProxyType(audit)


@dataclass(frozen=True)
class BlochBasis:
    miller: np.ndarray
    k_cart_bohr_inv: np.ndarray
    reciprocal_bohr_inv: np.ndarray
    row_map: np.ndarray
    source_miller: np.ndarray
    slot: object
    reference_source: object
    source_admission: object
    identity: str

    def validate(self, reference):
        slot = self.slot
        if dict(self.reference_source) != dict(reference.source):
            raise ValueError('reconstructed basis has a different XML authority')
        if (self.source_admission['source_sha256'] != slot.source_binding['sha256']
                or self.source_admission['ispin'] != slot.ispin
                or self.source_admission['ik'] != slot.source_ik
                or self.source_admission['source_rows_sha256'] != _digest(self.source_miller)):
            raise ValueError('reconstructed basis source admission mismatch')
        point, _, _ = reference.select(slot.source_ik, slot.ispin)
        saved = reference.saved_symmetries
        if saved is None or not 1 <= slot.operation.index <= saved.nsym:
            raise ValueError('basis operation is absent from saved crystal group')
        original = saved.operations[slot.operation.index-1]
        if (original.effective_time_reversal is not False
                or not np.array_equal(original.rotation, slot.operation.rotation)
                or not np.array_equal(original.fractional_translation, slot.operation.fractional_translation)
                or not np.array_equal(original.equivalent_atoms, slot.operation.equivalent_atoms)
                or slot.expected_rows != point.npw
                or self.source_miller.shape != (3, point.npw)
                or not np.allclose(reference.reciprocal_bohr_inv @ slot.source_k_fractional,
                                   point.k_cart_bohr_inv, atol=2e-13, rtol=0)
                or not np.allclose(slot.k_fractional, original.rotation @ slot.source_k_fractional
                                   + slot.reciprocal_shift, atol=2e-10, rtol=0)):
            raise ValueError('basis operation/source/fold is not the saved unitary route')
        expected = slot.operation.rotation @ self.source_miller - slot.reciprocal_shift[:, None]
        if (not np.array_equal(self.miller, expected)
                or not np.array_equal(self.row_map, np.arange(slot.expected_rows))
                or not np.array_equal(self.k_cart_bohr_inv, slot.k_cart_bohr_inv)
                or not np.allclose(self.reciprocal_bohr_inv, reference.reciprocal_bohr_inv, atol=2e-13, rtol=0)
                or self.identity != _basis_identity(self.miller, self.reciprocal_bohr_inv, slot)):
            raise ValueError('reconstructed row/k/gauge identity mismatch')
        return self


def _basis_identity(miller, reciprocal, slot):
    return sha256((_digest(miller, reciprocal, slot.k_cart_bohr_inv, slot.operation.rotation,
                           slot.operation.fractional_translation, slot.reciprocal_shift)
                   + slot.source_binding['sha256'] + slot.gauge_policy
                   + str((slot.ispin, slot.mesh_index, slot.source_ik, slot.operation.index))).encode()).hexdigest()


@dataclass(frozen=True)
class BlochFrame:
    basis: BlochBasis
    coefficients: np.ndarray
    eigenvalues_ry: np.ndarray
    raw_occupations: np.ndarray
    gauge_identity: str

    @property
    def ispin(self):
        return self.basis.slot.ispin


def transport_bloch_frame(reference, slot, source, *, source_admission=None):
    """Transport raw columns without diagonalization, normalization or band gauge."""
    audit = admit_source(reference, slot, source) if source_admission is None else source_admission
    if (audit['source_sha256'] != source.source['sha256'] or audit['ispin'] != slot.ispin
            or audit['ik'] != slot.source_ik):
        raise ValueError('cached source admission belongs to another spin/source')
    m = readonly(slot.operation.rotation @ source.miller - slot.reciprocal_shift[:, None])
    # Zero translations share the original immutable coefficients, no star-sized copy.
    f = slot.operation.fractional_translation
    coefficients = source.coefficients
    if np.any(f):
        phase = np.exp(-2j*np.pi*((source.miller.T + slot.source_k_fractional) @ f))
        coefficients = phase[:, None]*source.coefficients
        coefficients.setflags(write=False)
    identity = _basis_identity(m, source.reciprocal_bohr_inv, slot)
    basis = BlochBasis(m, slot.k_cart_bohr_inv, source.reciprocal_bohr_inv,
                       readonly(np.arange(source.igwx)), source.miller, slot,
                       reference.source, audit, identity)
    basis.validate(reference)
    PhysicalRhoFFT.create(reference.fft_smooth, m, basis.reciprocal_bohr_inv,
                          abs(np.linalg.det(reference.cell_bohr)))
    return BlochFrame(basis, coefficients, slot.eigenvalues_ry, slot.raw_occupations,
                      identity + ':raw-saved-columns')


def load_bloch_frame(reference, slot):
    return transport_bloch_frame(reference, slot, load_wavefunctions(slot.source_binding))
