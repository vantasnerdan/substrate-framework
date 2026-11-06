"""Two-coordinate centered scalar/force curvature at a self-conjugate q.

Inputs are per-primitive-cell free energies and conjugate forces, with amplitudes
in bohr. The background need not be ionically stationary: odd gradients cancel.
No ASR, mass weighting, clipping or post-hoc symmetry correction is performed.
"""
from decimal import Decimal

import numpy as np

DIRECTIONS = {'e1': (1, 0), 'e2': (0, 1), 'sum': (1, 1)}


def scalar_hessian(energies, step):
    """Six signed directional energies plus an independent center; keep decimals."""
    h2 = Decimal(str(step)) ** 2
    center = Decimal(str(energies['center']))
    paired = {name: Decimal(str(energies[name + '+']))
              + Decimal(str(energies[name + '-'])) for name in DIRECTIONS}
    h11 = (paired['e1'] - 2 * center) / h2
    h22 = (paired['e2'] - 2 * center) / h2
    h12 = (paired['sum'] - paired['e1'] - paired['e2'] + 2 * center) / (2 * h2)
    return np.array([[float(h11), float(h12)], [float(h12), float(h22)]])


def force_hessian(forces, step):
    """Differentiate both independently observed conjugate forces; no symmetrizing."""
    return np.column_stack([
        -(np.asarray(forces[name + '+']) - np.asarray(forces[name + '-'])) / (2 * step)
        for name in ('e1', 'e2')])


def conjugate_forces(cartesian_forces, signs):
    """Ni-x and H-y on ordered (Ni,H) copies; average per primitive cell."""
    forces = np.asarray(cartesian_forces).reshape(len(signs), 2, 3)
    return np.array([np.dot(signs, forces[:, 0, 0]),
                     np.dot(signs, forces[:, 1, 1])]) / len(signs)
