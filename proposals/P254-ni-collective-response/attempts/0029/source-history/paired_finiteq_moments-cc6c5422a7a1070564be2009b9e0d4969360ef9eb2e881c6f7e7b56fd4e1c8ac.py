"""Paired finite-q ionic moments and independent directional drives.

Only original source/target real-space transforms are shared: the Hamiltonian,
physical FFT backend, finite-q support and nonlocal gradient primitive are not
changed. The two dyads may have different ranks and endpoint PW dimensions.
"""

import numpy as np
from scipy.linalg.blas import zgemm


def _add_nonlocal_direction(coupling, direction, vectors, result, *, adjoint):
    """Add the original ionic nonlocal action, without repeating local FFTs."""
    u = direction.conj() if adjoint else direction
    pin, pout = ((coupling.pt, coupling.ps) if adjoint
                 else (coupling.ps, coupling.pt))
    width = vectors.shape[1]
    for source, target in coupling.species:
        incoming, outgoing = (target, source) if adjoint else (source, target)
        atoms = np.asarray(incoming.atom_indices, dtype=int)-1
        shape = (len(atoms), incoming.active_d.shape[0], width)
        d = incoming.active_d.conj().T if adjoint else incoming.active_d
        base = zgemm(1., incoming.beta, vectors, trans_a=2).reshape(shape)
        dbase = np.einsum('ij,ajb->aib', d, base, optimize=True)
        for axis in range(3):
            ua = u[atoms, axis, None, None]
            if not np.any(ua):
                continue
            derivative = zgemm(1., incoming.beta, pin[:, axis, None]*vectors,
                               trans_a=2).reshape(shape)
            weighted = np.einsum('ij,ajb->aib', d, derivative, optimize=True)*ua
            result += 1j*zgemm(1., outgoing.beta, weighted.reshape(-1, width))
            result -= 1j*pout[:, axis, None]*zgemm(
                1., outgoing.beta, (dbase*ua).reshape(-1, width))


def _accumulate_pair(coupling, left, right, action, channel, gradient, *, adjoint):
    """Accumulate one dyad, reusing its driven leg's inverse transform."""
    bridge = coupling.bridge
    outrows = (coupling.edge.source_dimension if adjoint
               else coupling.edge.target_dimension)
    drive = np.empty((outrows, left.shape[1]), dtype=np.complex128)
    multiplier = bridge.multiplication_phase*action.field
    if adjoint:
        multiplier = multiplier.conj()
    for start in range(0, left.shape[1], bridge.block_size):
        sl = slice(start, start+bridge.block_size)
        lblock, rblock = left[:, sl], right[:, sl]
        lwave = bridge.target_fft.from_g(lblock)
        rwave = bridge.source_fft.from_g(rblock)
        channel += np.sum(lwave*rwave.conj(), axis=1)
        gradient += coupling._nonlocal_gradient_block(lblock, rblock)
        vectors, wave, fft = ((lblock, lwave, bridge.source_fft) if adjoint
                              else (rblock, rwave, bridge.target_fft))
        result = fft.to_g(multiplier[:, None]*wave)
        _add_nonlocal_direction(coupling, action.direction, vectors, result,
                                adjoint=adjoint)
        drive[:, sl] = result
    return drive


def _finish_moments(coupling, density, gradient, weight):
    """Apply the original phase, volume, spin, weight and positive transpose."""
    bridge = coupling.bridge
    channel = density[:, coupling.edge.ispin-1]
    channel *= weight*bridge.multiplication_phase.conj()/bridge.source_fft.omega
    gradient *= weight
    sampled = bridge.density_fft.to_g(channel)
    ref = coupling.inputs.reference
    for atom, name in enumerate(ref.atom_names):
        weighted = (sampled*coupling.vloc[name].conj()
                    * np.exp(1j*(coupling.gq @ ref.positions_bohr[:, atom])))
        gradient[atom] += 1j*bridge.source_fft.omega*(weighted @ coupling.gq)
    return density, gradient.reshape(-1)


def paired_moments_and_drive(coupling, left, right, left2, right2, direction, *, weight):
    """Return summed ``(density, gradient, drive, adjoint_drive)``.

    Moments equal ``coupling.moments(left, right, weight=weight)`` plus
    ``coupling.moments(left2, right2, weight=weight)``. Drives equal
    ``coupling.prepare(direction).apply(right)`` and its independent
    ``adjoint_apply(left2)``. Weight applies only to moments. Each pair has
    equal internal column counts; their ranks need not match one another.
    """
    bridge = coupling.bridge
    weight = bridge._weight(weight)
    left = bridge._matrix(left, coupling.edge.target_dimension)
    right = bridge._matrix(right, coupling.edge.source_dimension)
    left2 = bridge._matrix(left2, coupling.edge.target_dimension)
    right2 = bridge._matrix(right2, coupling.edge.source_dimension)
    if left.shape[1] != right.shape[1] or left2.shape[1] != right2.shape[1]:
        raise ValueError('identity coefficient requires equal column counts within each dyad')
    action = coupling.prepare(direction)
    density = np.zeros((bridge.source_fft.nnr, 2), dtype=np.complex128)
    gradient = np.zeros((coupling.nat, 3), dtype=np.complex128)
    channel = density[:, coupling.edge.ispin-1]
    drive = _accumulate_pair(coupling, left, right, action, channel, gradient,
                             adjoint=False)
    adjoint_drive = _accumulate_pair(coupling, left2, right2, action, channel,
                                     gradient, adjoint=True)
    density, gradient = _finish_moments(coupling, density, gradient, weight)
    return density, gradient, drive, adjoint_drive
