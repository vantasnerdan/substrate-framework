"""Original coefficient-absorbed ionic blocks with request-local overlap reuse."""
import numpy as np
from scipy.linalg.blas import zgemm


def _overlap(cache, bank, block, key, *, momenta=None, axis=None):
    identity = (id(bank), key, block.shape[1], id(momenta), axis)
    if identity not in cache:
        vectors = block if momenta is None else momenta[:, axis, None]*block
        cache[identity] = zgemm(1., bank, vectors, trans_a=2)
    return cache[identity]


def ionic_gradient_block(coupling, left, right, left_key, right_key, cache):
    """Match BlochFiniteQIonic._nonlocal_gradient_block's arithmetic order.

    D has already been absorbed into the SAME bounded PW factors used by the
    orbital FFTs. Applying it after projection is algebraically equivalent but
    changes cancellation-sensitive material gradients. No coefficient, UPF
    channel, momentum derivative or original ordered term is inferred away.
    """
    result = np.zeros((coupling.nat, 3), dtype=np.complex128)
    width = left.shape[1]
    for source, target in coupling.species:
        atoms = np.asarray(source.atom_indices, dtype=int)-1
        shape = (len(atoms), source.active_d.shape[0], width)
        bl = _overlap(cache, target.beta, left, left_key).reshape(shape)
        br = _overlap(cache, source.beta, right, right_key).reshape(shape)
        channel_key = ('channel', id(source.active_d))
        if channel_key not in cache:
            cache[channel_key] = source.active_d.conj().T
        channel = cache[channel_key]
        weighted_key = ('weighted', id(target.beta), id(source.active_d), left_key, width)
        if weighted_key not in cache:
            cache[weighted_key] = np.einsum('ij,ajb->aib', channel, bl, optimize=True)
        dbl = cache[weighted_key]
        for axis in range(3):
            dl = _overlap(cache, target.beta, left, left_key,
                          momenta=coupling.pt, axis=axis).reshape(shape)
            dr = _overlap(cache, source.beta, right, right_key,
                          momenta=coupling.ps, axis=axis).reshape(shape)
            derivative_key = ('weighted-first', id(target.beta), id(source.active_d),
                              left_key, width, id(coupling.pt), axis)
            if derivative_key not in cache:
                cache[derivative_key] = np.einsum('ij,ajb->aib', channel, dl, optimize=True)
            result[atoms, axis] = 1j*(
                np.einsum('aib,aib->a', br.conj(), cache[derivative_key])
                - np.einsum('aib,aib->a', dr.conj(), dbl))
    return result
