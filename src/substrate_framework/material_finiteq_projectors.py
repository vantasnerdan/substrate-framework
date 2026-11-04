"""Stream the two projector legs of the declared straight-line Peierls vertex.

For A(r)=a exp(i q.r), differentiate exp[-i integral_(r'->r) A.dl].
At line parameter t the Fourier kernel is V(pt-tq, ps+(1-t)q),
and its derivative is (a.d_pt+a.d_ps)V. The common atomic placement
phases cancel under that derivative, but their q-dependent value is retained.
The right derivative bank uses conjugate(a), since the bank is later daggered.
No continuity reconstruction or all-electron transverse claim is involved.
"""
from types import SimpleNamespace

import numpy as np
from scipy.linalg.blas import zgemm

from .qe_material_field_derivatives import UPFRadial, centered_derivatives
from .qe_nc_projectors import place_projectors


def directional_bank(radial, momenta, direction):
    """Centered value and complex-linear directional derivative; one radial sum."""
    values = radial.evaluate(np.einsum('pi,pi->p', momenta, momenta))
    cached = SimpleNamespace(tables=radial.tables, evaluate=lambda s: values)
    base, real, _, _ = centered_derivatives(cached, momenta, direction.real)
    if np.any(direction.imag):
        imag = centered_derivatives(cached, momenta, direction.imag)[1]
        real = real+1j*imag
    return base, real


def placed_blocks(inputs, species, pt, ps, q, direction, t, atom_block_size):
    """Yield (Bt, Bt_a, Bs, Bs_conj(a), D), bounded by atom block width."""
    radial = UPFRadial(inputs.tables[species.name])
    u, v = pt-t*q, ps+(1-t)*q
    target = directional_bank(radial, u, direction)
    source = directional_bank(radial, v, direction.conj())
    nh = species.active_d.shape[0]
    atoms = species.atom_indices
    for start in range(0, len(atoms), atom_block_size):
        selected = atoms[start:start+atom_block_size]
        banks = [np.empty((len(p), len(selected)*nh), complex, order='F')
                 for p in (u, u, v, v)]
        for slot, atom in enumerate(selected):
            tau = inputs.reference.positions_bohr[:, atom-1]
            sl = slice(slot*nh, (slot+1)*nh)
            for bank, centered, p in zip(banks, (*target, *source), (u, u, v, v)):
                bank[:, sl] = place_projectors(centered, p, tau)
        yield (*banks, species.active_d)


def weighted_overlap(bank, vectors, d):
    width, nh = vectors.shape[1], d.shape[0]
    overlap = zgemm(1., bank, vectors, trans_a=2).reshape(-1, nh, width)
    return np.einsum('ij,ajb->aib', d, overlap, optimize=True).reshape(-1, width)


def forward_block(bt, dt, bs, ds, d, vectors):
    """Both nonlocal motion legs, forward source -> target."""
    return (zgemm(1., dt, weighted_overlap(bs, vectors, d))
            + zgemm(1., bt, weighted_overlap(ds, vectors, d)))


def adjoint_block(bt, dt, bs, ds, d, vectors):
    """Independently encoded dagger; never substitute a refolded reverse edge."""
    dh = d.conj().T
    return (zgemm(1., bs, weighted_overlap(dt, vectors, dh))
            + zgemm(1., ds, weighted_overlap(bt, vectors, dh)))
