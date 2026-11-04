"""Preparation-time one-body operator from the existing Gamma coupled dual.

The dual is already pulled back: do not propagate this operator again or
multiply by another Fourier/source/spin factor. Basis columns are the actual
saved active columns expressed in the consumer's declared Hilbert coordinates.
"""
import numpy as np


def homogeneous_electronic_operators(initial_dual, active_basis, *, mechanical_size,
                                     source_count, source_slot):
    """Return [probe,row,row] operators; ordinary transpose of each dual block.

    This is the single-sector Gamma adapter. General-q edge embeddings are not
    this contract. Physical adjoints are computed from the resulting explicit
    matrices, not by conjugating the nonconjugating dual pairing.
    """
    basis = np.asarray(active_basis, dtype=np.complex128)
    dual = np.asarray(initial_dual, dtype=np.complex128)
    if dual.ndim == 1:
        dual = dual[:, None]
    if basis.ndim != 2 or dual.ndim != 2:
        raise ValueError('active basis and coupled dual must be matrices')
    active_count = basis.shape[1]
    if (source_count <= 0 or not 0 <= source_slot < source_count
            or mechanical_size < 0
            or dual.shape[0] != mechanical_size+source_count*active_count**2):
        raise ValueError('exact coupled layout and saved source slot required')
    coefficients = dual[mechanical_size:].reshape(
        source_count,active_count,active_count,dual.shape[1],order='C')[source_slot]
    return np.stack([basis@coefficients[:,:,probe].T@basis.conj().T
                     for probe in range(dual.shape[1])])
