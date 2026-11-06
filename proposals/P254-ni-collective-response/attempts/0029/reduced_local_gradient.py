"""Exact all-coordinate local Fourier gradient matrix, built once per q/support."""
from time import perf_counter
import numpy as np


def local_gradient_matrix(inputs, support, vloc):
    """M[G,3*atom+axis] = i Omega conj(vloc) exp(+i gq.R) gq_axis.

    This is the physical positive-gradient transpose, not a projected adjoint
    imposed on data. The forward local action uses its conjugated entries.
    """
    started = perf_counter()
    gq = support.wavevectors.T
    ref = inputs.reference
    nat = len(ref.atom_names)
    matrix = np.empty((len(gq), 3*nat), complex, order='F')
    for atom, name in enumerate(ref.atom_names):
        weighted = 1j*support.omega*vloc[name].conj()*np.exp(1j*(gq @ ref.positions_bohr[:, atom]))
        matrix[:, 3*atom:3*atom+3] = weighted[:, None]*gq
    matrix.setflags(write=False)
    return matrix, dict(local_gradient_build_seconds=perf_counter()-started,
                        local_gradient_bytes=matrix.nbytes,
                        local_gradient_coordinates=3*nat,
                        local_gradient_support=len(gq))
