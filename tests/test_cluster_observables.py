from itertools import product

import numpy as np
import pytest

from substrate_framework.cluster_observables import sample_periodic_density


def test_density_sampling_preserves_skew_cell_axes_and_periodic_images():
    shape = np.array([20, 24, 30])
    i, j, k = np.ogrid[:20, :24, :30]
    density = 10 + np.cos(2*np.pi*i/20) + 2*np.sin(2*np.pi*j/24) + 3*np.cos(4*np.pi*k/30)
    cell = np.array([[3.2, .2, .4], [1.1, 4.3, .1], [.3, .8, 5.1]])
    grid_coordinate = np.array([3.25, 6.5, 12.75])
    lower = grid_coordinate.astype(int)
    fraction = grid_coordinate-lower
    expected = 0.0
    for corner in product((0, 1), repeat=3):
        corner = np.array(corner)
        weight = np.prod(np.where(corner, fraction, 1-fraction))
        index = (lower+corner) % shape
        value = (10 + np.cos(2*np.pi*index[0]/20) + 2*np.sin(2*np.pi*index[1]/24)
                 + 3*np.cos(4*np.pi*index[2]/30))
        expected += weight*value
    point = (grid_coordinate/shape) @ cell
    assert sample_periodic_density(density, cell, point) == pytest.approx(expected, abs=1e-12)
    translated = point+np.array([4, -3, 2]) @ cell
    assert sample_periodic_density(density, cell, translated) == pytest.approx(expected, abs=1e-12)
