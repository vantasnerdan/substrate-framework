"""Durable full complex128 payloads, with typed metadata for exact zeros.

Only exactly zero freshly computed data becomes an exact-zero descriptor. Every
nonzero value, however tiny, retains its actual durable NPY bytes and binding.
"""
from pathlib import Path

import numpy as np
from substrate_framework.partitioned_q_storage import sync_directory, sync_file


def exact_zero_payload(shape):
    """Describe all original rows/columns without inventing a file or byte cost."""
    if (not isinstance(shape, (list, tuple)) or len(shape) != 2
            or any(isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer))
                   or n < 1 for n in shape)):
        raise ValueError('positive full-row/source-column integer shape required')
    return dict(kind='exact-zero', shape=[int(n) for n in shape], dtype='complex128')


def is_exact_zero_payload(entry):
    """Return whether entry is a complete, strictly typed exact-zero descriptor."""
    if (not isinstance(entry, dict) or set(entry) != {'kind', 'shape', 'dtype'}
            or entry['kind'] != 'exact-zero' or entry['dtype'] != 'complex128'):
        return False
    try:
        exact_zero_payload(entry['shape'])
    except ValueError:
        return False
    return True


def write_payload(path, value, budget):
    """Publish finite full complex128 data; exact zeros create no target file.

    Existing targets are always refused, including for an exact-zero result.
    Metadata publication is separately charged by the checkpoint/journal writer.
    """
    if (not isinstance(value, np.ndarray) or value.ndim != 2
            or min(value.shape) < 1 or value.dtype != np.complex128
            or not np.isfinite(value).all()):
        raise ValueError('actual finite full-row complex128 matrix required')
    budget.check()
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    if not np.any(value):
        return exact_zero_payload(value.shape)
    binding = budget.array(path, value)
    sync_file(binding['path'])
    sync_directory(Path(binding['path']).parent)
    return binding
