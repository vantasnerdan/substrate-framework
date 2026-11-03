"""Paired Gamma field traces and complex-linear physical first-field drives.

The two bilinear pairs have the exact ordering of StreamedJointResponse:
<L_j right, left> + <L_j right2, left2>. Only the traces carry weight.
No sqrt(N), dual sign, volume factor, FFT replacement or field contact is added.
"""
import numpy as np
from scipy.linalg.blas import zgemm

from .material_mixed_field_ionic import first_cartesian_jets
from .material_physical_operators import PhysicalFieldVertices, _vectors, _weighted_overlap
from .qe_material_field_derivatives import UPFRadial
from .qe_nc_projectors import place_projectors


def _paired_vectors(vertex, left, right, left2, right2):
    if not isinstance(vertex, PhysicalFieldVertices):
        raise TypeError('original PhysicalFieldVertices required')
    rows = len(vertex.momenta)
    l, _ = _vectors(left, rows)
    r, single = _vectors(right, rows)
    l2, single2 = _vectors(left2, rows)
    r2, _ = _vectors(right2, rows)
    if l.shape != r.shape or l2.shape != r2.shape:
        raise ValueError('matching widths within each independent bilinear pair required')
    return (l, r, l2, r2), single, single2


def _kinetic_response(momenta, pairs, direction, block_size):
    left, right, left2, right2 = pairs
    diagonal = np.zeros(len(momenta), complex)
    for l, r in ((left, right), (left2, right2)):
        for start in range(0, r.shape[1], block_size):
            sl = slice(start, start+block_size)
            diagonal += np.einsum('pi,pi->p', r[:, sl].conj(), l[:, sl])
    ell = 2*(momenta.T @ diagonal)
    kinetic = 2*(momenta @ direction)
    drive = kinetic[:, None]*right
    adjoint_drive = kinetic.conj()[:, None]*left2
    return ell, drive, adjoint_drive


def _place_first(vertex, species, center, out):
    nh = species.active_d.shape[0]
    positions = vertex.hamiltonian.background.inputs.reference.positions_bohr
    for slot, atom in enumerate(species.atom_indices):
        out[:, slot*nh:(slot+1)*nh] = place_projectors(
            center, vertex.momenta, positions[:, atom-1])


def _overlaps(bank, pairs):
    return tuple(zgemm(1., bank, value, trans_a=2) for value in pairs)


def _add_drive(species, derivative, weighted_base, first_overlap, out, block_size):
    weighted_first = _weighted_overlap(species, first_overlap, out.shape[1])
    for start in range(0, out.shape[1], block_size):
        sl = slice(start, start+block_size)
        out[:, sl] += zgemm(1., derivative, weighted_base[:, sl])
        out[:, sl] += zgemm(1., species.beta, weighted_first[:, sl])


def _add_species_response(vertex, species, pairs, direction, ell, drive, adjoint_drive):
    inputs = vertex.hamiltonian.background.inputs
    jets = first_cartesian_jets(UPFRadial(inputs.tables[species.name]), vertex.momenta)
    base = _overlaps(species.beta, pairs)
    weighted = tuple(_weighted_overlap(species, value, value.shape[1]) for value in base)
    first_right = np.zeros_like(base[1])
    first_left2 = np.zeros_like(base[2])
    derivative = np.empty_like(species.beta, order='F')
    forward_center = np.zeros_like(jets[1])
    adjoint_center = np.zeros_like(jets[1])
    for axis, center in enumerate(jets[1:]):
        _place_first(vertex, species, center, derivative)
        first_left, first_r, first_l2, first_right2 = _overlaps(derivative, pairs)
        # D is the original Hermitian, possibly complex/off-diagonal channel D.
        # This is <B_j D B^H r + B D B_j^H r, l>, not its conjugate.
        ell[axis] += (np.vdot(first_r, weighted[0])
                      + np.vdot(weighted[1], first_left)
                      + np.vdot(first_right2, weighted[2])
                      + np.vdot(weighted[3], first_l2))
        first_right += direction[axis]*first_r
        first_left2 += direction[axis].conjugate()*first_l2
        forward_center += direction[axis]*center
        adjoint_center += direction[axis].conjugate()*center
    # For complex a the second projector term uses sum(a_j B_j^H),
    # NOT (sum(a_j B_j))^H. Independent overlaps above preserve this.
    _place_first(vertex, species, forward_center, derivative)
    _add_drive(species, derivative, weighted[1], first_right, drive, vertex.block_size)
    _place_first(vertex, species, adjoint_center, derivative)
    _add_drive(species, derivative, weighted[2], first_left2, adjoint_drive, vertex.block_size)


def paired_field_moments_and_drive(vertex, left, right, left2, right2,
                                   field_direction, *, weight):
    """Return ``(ell, L(a) right, L(a)^H left2)`` without three action frames.

    ``ell[j] = weight * (vdot(L_j right, left) +
    vdot(L_j right2, left2))``. Pair widths may differ. Vector drives retain
    vector shape; frame drives retain their original full-PW column count.
    Cartesian radial first jets and placed projector workspace are shared by
    both dyads and both drives, species by species. Endpoint Hermiticity never
    implies Hermiticity of an arbitrary complex tangent field direction.
    """
    pairs, single, single2 = _paired_vectors(vertex, left, right, left2, right2)
    direction = np.asarray(field_direction, dtype=complex)
    if direction.shape != (3,) or not np.isfinite(direction).all():
        raise ValueError('finite complex Cartesian field direction (3,) required')
    w = np.asarray(weight)
    if w.ndim != 0 or not np.isfinite(w):
        raise ValueError('finite scalar original sector weight required')
    ell, drive, adjoint_drive = _kinetic_response(
        vertex.momenta, pairs, direction, vertex.block_size)
    for species in vertex.hamiltonian.factors.species:
        _add_species_response(vertex, species, pairs, direction, ell, drive, adjoint_drive)
    return (weight*ell, drive[:, 0] if single else drive,
            adjoint_drive[:, 0] if single2 else adjoint_drive)


class CartesianPhysicalFieldAxes:
    """Indexable original Cartesian actions, prepared only when requested."""
    def __init__(self, vertex):
        self.vertex = vertex
        self._axes = [None, None, None]

    def __len__(self):
        return 3

    def __getitem__(self, axis):
        value = self._axes[axis]
        if value is None:
            value = self.vertex.prepare(np.eye(3)[axis])
            self._axes[axis] = value
        return value

    def __iter__(self):
        return (self[axis] for axis in range(3))
