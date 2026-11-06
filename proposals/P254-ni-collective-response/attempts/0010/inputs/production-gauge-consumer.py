"""Explicit ordinary-NC screened operator views; never infer a chemical-potential shift."""
import numpy as np


def screened_native_matrix(record):
    """Return the untouched same-gauge ordinary/EPA matrix, without copying it."""
    if record['kind'] != 'SCREENED_NATIVE':
        raise ValueError('A screened gauge view requires SCREENED_NATIVE data')
    return record['matrix']


def _shift(header, record):
    value = header['fermi_shift_ry_per_bohr'][record['perturbation'] - 1]
    if not header['fermi_shift_applied'] and value != 0:
        raise ValueError('Fermi shift outside its actual native applied branch')
    return value


def screened_hamiltonian_matrix(header, record):
    """Reconstruct absolute dH only on demand; never modify the native matrix."""
    native = screened_native_matrix(record)
    shift = _shift(header, record)
    if shift == 0:
        return native
    absolute = native.copy()
    absolute.flat[::absolute.shape[0] + 1] += shift
    return absolute


def screened_hamiltonian_diagonal(header, record):
    """Reconstruct the absolute diagonal without allocating another full matrix."""
    return np.diag(screened_native_matrix(record)) + _shift(header, record)
