"""Proposed raw-current to global canonical-field generalized-force map.

J_Ha=-L_Ry/2 and Pi-dot_el=-N e L give the source-vector map M=2e.
For the existing per-cell sum(w C_J), the proposed extensive covariance is
N M C_J M.T, not N**2 M C_J M.T. Independent physical verification is
required; these functions perform no state admission or physical computation.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
from substrate_framework.prepared_reference import array


def current_force_matrix(polarizations):
    """Return real M=2e, using the original Maxwell polarization convention."""
    e = np.asarray(polarizations)
    if e.ndim != 2 or e.shape[1] != 3 or not 1 <= len(e) <= 3:
        raise ValueError('one to three explicitly chosen Cartesian polarizations required')
    e = array(e, e.shape, 'polarizations', real=True)
    if not np.allclose(e @ e.T, np.eye(len(e)), rtol=1e-13, atol=1e-13):
        raise ValueError('orthonormal polarizations required; no normalization repair')
    return 2 * e


def canonical_force_covariance(current_covariance, ncells, polarizations):
    """Map one ordered raw complex 3x3 matrix without symmetrizing or repair.

    The input is sum over explicit spin/internal-k sectors of w*C_J, with
    quadrature already applied once. Output is covariance of the electronic
    generalized force entering global canonical Pi-dot, not Pi-state noise
    or a mass-normalized acceleration. Partial sums are not renormalized.
    """
    if type(ncells) is not int or ncells < 1:
        raise ValueError('positive integer number of physical cells required')
    covariance = array(current_covariance, (3, 3), 'raw current covariance')
    matrix = current_force_matrix(polarizations)
    return ncells * (matrix @ covariance @ matrix.T)
