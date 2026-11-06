"""Three physical Cartesian L traces of one borrowed low-rank density dyad.

For X = weight * left @ coefficient @ right.H, return Tr(L_j X).
Only bounded PW-column temporaries and projector/column overlaps are formed.
The original UPF channel matrix and the occupation coefficient are distinct;
neither is diagonalized, inverted, conjugated, or silently made Hermitian.
"""
import numpy as np
from scipy.linalg.blas import zgemm

from substrate_framework.material_mixed_field_ionic import first_cartesian_jets
from substrate_framework.material_physical_operators import (
    PhysicalFieldVertices, _vectors, _weighted_overlap)
from substrate_framework.paired_field_moments import _place_first
from substrate_framework.qe_material_field_derivatives import UPFRadial


def _overlap(bank, frame, block_size):
    """B.H @ frame without a whole-frame BLAS layout-conversion copy."""
    result = np.empty((bank.shape[1], frame.shape[1]), complex)
    for start in range(0, frame.shape[1], block_size):
        sl = slice(start, start + block_size)
        result[:, sl] = zgemm(1., bank, frame[:, sl], trans_a=2)
    return result


def _kinetic_moments(momenta, left, right, coefficient, block_size):
    diagonal = np.zeros(len(momenta), complex)
    # Absorb the independent coefficient only into bounded column blocks.
    if right.shape[1] <= left.shape[1]:
        for start in range(0, right.shape[1], block_size):
            sl = slice(start, start + block_size)
            weighted = left @ coefficient[:, sl]
            diagonal += np.einsum('pi,pi->p', right[:, sl].conj(), weighted)
    else:
        for start in range(0, left.shape[1], block_size):
            sl = slice(start, start + block_size)
            weighted = right @ coefficient[sl, :].conj().T
            diagonal += np.einsum('pi,pi->p', weighted.conj(), left[:, sl])
    return 2 * (momenta.T @ diagonal)


def _species_moments(vertex, species, left, right, coefficient):
    size = vertex.block_size
    base_left = _overlap(species.beta, left, size)
    base_right = base_left if right is left else _overlap(species.beta, right, size)
    weighted_base = _weighted_overlap(species, base_left, left.shape[1]) @ coefficient
    inputs = vertex.hamiltonian.background.inputs
    jets = first_cartesian_jets(UPFRadial(inputs.tables[species.name]), vertex.momenta)
    derivative = np.empty_like(species.beta, order='F')
    result = np.empty(3, complex)
    for axis, center in enumerate(jets[1:]):
        _place_first(vertex, species, center, derivative)
        first_left = _overlap(derivative, left, size)
        first_right = first_left if right is left else _overlap(derivative, right, size)
        weighted_first = _weighted_overlap(species, first_left, left.shape[1]) @ coefficient
        result[axis] = (np.vdot(first_right, weighted_base)
                        + np.vdot(base_right, weighted_first))
    return result


def physical_field_dyad_moments(vertex, left, right, coefficient, *, weight):
    """Return complex Cartesian ``(3,)`` physical L moments, weight ONCE.

    ``vertex`` is the actual endpoint PhysicalFieldVertices, with qualified
    direct-physical UPF banks. ``coefficient`` may be singular, complex, and
    rectangular; its axes are left columns then right columns. No field
    direction is applied and no full-PW field-action frame is allocated.
    """
    if not isinstance(vertex, PhysicalFieldVertices):
        raise TypeError('original PhysicalFieldVertices required')
    l, _ = _vectors(left, len(vertex.momenta))
    r = l if right is left else _vectors(right, len(vertex.momenta))[0]
    d = np.asarray(coefficient)
    if (d.shape != (l.shape[1], r.shape[1]) or min(d.shape) < 1
            or not np.isfinite(d).all()):
        raise ValueError('finite original coefficient[left columns,right columns] required')
    if not np.isscalar(weight) or not np.isreal(weight) or not np.isfinite(weight):
        raise ValueError('explicit finite real original source weight required')
    result = _kinetic_moments(vertex.momenta, l, r, d, vertex.block_size)
    for species in vertex.hamiltonian.factors.species:
        result += _species_moments(vertex, species, l, r, d)
    return weight * result
