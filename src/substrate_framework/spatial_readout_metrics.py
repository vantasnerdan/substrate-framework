"""Selected scalar-observable weights, not energy or full-state localization."""
from itertools import product

import numpy as np


def periodic_distances(cells, other, mesh, cell_bohr, *, offset_bohr=(0., 0., 0.)):
    """Exact closest-image distances; cell_bohr has lattice vectors as COLUMNS.

    Enumerate a provably sufficient integer box from the smallest singular
    value of the BvK supercell. Unlike rounding fractional coordinates alone,
    this also works for skew cells. offset_bohr is the selected probe position
    minus the anchor position inside their original primitive cells.
    """
    cells, other, mesh = np.asarray(cells), np.asarray(other), np.asarray(mesh)
    lattice = np.asarray(cell_bohr, float)
    offset = np.asarray(offset_bohr, float)
    if (lattice.shape != (3, 3) or not np.isfinite(lattice).all()
            or mesh.shape != (3,) or np.any(mesh <= 0)
            or offset.shape != (3,) or not np.isfinite(offset).all()):
        raise ValueError('finite nonsingular original geometry and positive mesh required')
    supercell = lattice @ np.diag(mesh)
    sigma = float(np.linalg.svd(supercell, compute_uv=False)[-1])
    if sigma <= 0:
        raise ValueError('singular original cell')
    fractional_offset = np.linalg.solve(supercell, offset)
    distances = []
    for delta in cells-other:
        fractional = delta/mesh+fractional_offset
        nearest = np.rint(fractional)
        upper = float(np.linalg.norm(supercell @ (fractional-nearest)))
        radius = upper/sigma
        # One extra integer protects the bound against floating roundoff.
        ranges = [range(int(np.floor(x-radius))-1, int(np.ceil(x+radius))+2)
                  for x in fractional]
        best = upper
        for image in product(*ranges):
            best = min(best, float(np.linalg.norm(supercell @ (fractional-image))))
        distances.append(best)
    return np.asarray(distances)


def selected_weight(values, distances, *, value_units):
    """P=(sum |z|^2)^2/sum |z|^4 effective cells; RMS distance in bohr.

    Arrays are one selected scalar per cell, NOT all atoms or a state norm.
    Per-cell distance/weight rows resolve distance without arbitrary binning.
    Zero signal has undefined normalized weights, participation and extent.
    """
    z, r = np.asarray(values, complex), np.asarray(distances, float)
    if (z.ndim != 1 or not len(z) or r.shape != z.shape
            or not np.isfinite(z).all() or not np.isfinite(r).all() or np.any(r < 0)):
        raise ValueError('one finite scalar and nonnegative distance per cell required')
    scale = float(np.max(np.abs(z)))
    weights = np.abs(z)**2
    base = dict(definition='w_R=|z_R|^2; P=1/sum((w_R/sum(w))^2); rms=sqrt(sum(w_R*r_R^2)/sum(w))',
                value_units=value_units, weight_units='('+value_units+')^2',
                squared_weight_by_cell=weights, distance_bohr_by_cell=r,
                total_squared_weight=float(np.sum(weights)),
                maximum_sampled_distance_bohr=float(np.max(r)))
    if scale == 0:
        return dict(base, defined=False, reason='exact zero selected-observable signal',
                    normalized_weight_by_cell=None, participation_cells=None,
                    participation_fraction=None, rms_distance_bohr=None,
                    outermost_distance_weight_fraction=None)
    scaled = np.abs(z/scale)**2
    probability = scaled/np.sum(scaled)
    participation = float(1/np.sum(probability**2))
    outer = np.isclose(r, np.max(r), rtol=1e-12, atol=1e-12)
    return dict(base, defined=True, normalized_weight_by_cell=probability,
                participation_cells=participation, participation_fraction=participation/len(z),
                rms_distance_bohr=float(np.sqrt(probability @ r**2)),
                outermost_distance_weight_fraction=float(np.sum(probability[outer])))
