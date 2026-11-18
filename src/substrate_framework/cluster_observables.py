"""Real-space observables for source-bound material-cluster calculations."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates


def sample_periodic_density(density, cell_A, point_A) -> float:
    """Trilinear sample of a node-centred periodic grid, including skew cells.

    Grid index (i,j,k) represents (i/Nx,j/Ny,k/Nz) @ cell_A.
    Density units are preserved. This is interpolation, not a density model.
    """
    array = np.asarray(density)
    if array.ndim != 3:
        raise ValueError('Density must be a three-dimensional periodic grid')
    fractional = np.linalg.solve(np.asarray(cell_A).T, np.asarray(point_A)) % 1.0
    coordinate = fractional * array.shape
    return float(map_coordinates(array, coordinate[:, None], order=1,
                                 mode='grid-wrap', prefilter=False)[0])
