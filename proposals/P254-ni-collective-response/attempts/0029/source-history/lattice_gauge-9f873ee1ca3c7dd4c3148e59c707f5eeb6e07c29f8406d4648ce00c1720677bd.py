"""Generic observables for square periodic two-dimensional lattice links.

Authority status: conditional, unpromoted numerical infrastructure linked to
open goal issue #28 (vantasnerdan/substrate-framework). The functions validate
finite unitary link matrices and compute only the objects named by their APIs.
They do not claim that an input represents a nontrivial bundle, prescribed
total flux, constant curvature, a minimum-action background, or a continuum
field configuration.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

LinkField = tuple[np.ndarray, np.ndarray]


def _validated_links(links: Sequence[np.ndarray]) -> LinkField:
    """Return complex link arrays after enforcing the public input contract."""

    if len(links) != 2:
        raise ValueError(f"links must contain two lattice directions; got {len(links)}")
    first = np.asarray(links[0], dtype=complex)
    second = np.asarray(links[1], dtype=complex)
    if first.ndim != 4 or second.ndim != 4:
        raise ValueError("each link array must have shape (N, N, d, d)")
    if first.shape != second.shape:
        raise ValueError(
            f"link directions must have identical shapes; got {first.shape} and {second.shape}"
        )
    n_first, n_second, rows, columns = first.shape
    if n_first != n_second or n_first < 2:
        raise ValueError(f"links require a square periodic lattice with N >= 2; got {first.shape}")
    if rows != columns or rows < 1:
        raise ValueError(f"link matrices must be nonempty and square; got {first.shape[-2:]}")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("link matrices must contain only finite values")

    identity = np.eye(rows, dtype=complex)
    for direction, field in enumerate((first, second)):
        products = field @ np.swapaxes(field.conj(), -1, -2)
        if not np.allclose(products, identity, rtol=1e-10, atol=1e-10):
            raise ValueError(f"links in direction {direction} must be unitary")
    return first, second


def plaquette_holonomies(links: Sequence[np.ndarray]) -> np.ndarray:
    """Return the oriented elementary plaquette holonomy at every site."""

    first, second = _validated_links(links)
    n_side = first.shape[0]
    plaquettes = np.empty_like(first)
    for n1 in range(n_side):
        for n2 in range(n_side):
            plaquettes[n1, n2] = (
                first[n1, n2]
                @ second[(n1 + 1) % n_side, n2]
                @ first[n1, (n2 + 1) % n_side].conj().T
                @ second[n1, n2].conj().T
            )
    return plaquettes


def cycle_holonomy(
    links: Sequence[np.ndarray],
    direction: int,
    fixed_coordinate: int,
) -> np.ndarray:
    """Return a Wilson holonomy around one periodic coordinate cycle.

    ``direction=0`` traverses the first coordinate at fixed second coordinate;
    ``direction=1`` traverses the second coordinate at fixed first coordinate.
    """

    first, second = _validated_links(links)
    n_side = first.shape[0]
    dimension = first.shape[-1]
    if direction not in (0, 1):
        raise ValueError(f"direction must be 0 or 1; got {direction!r}")
    if fixed_coordinate not in range(n_side):
        raise ValueError(
            f"fixed_coordinate must be in range({n_side}); got {fixed_coordinate!r}"
        )

    holonomy = np.eye(dimension, dtype=complex)
    if direction == 0:
        for n1 in range(n_side):
            holonomy = holonomy @ first[n1, fixed_coordinate]
    else:
        for n2 in range(n_side):
            holonomy = holonomy @ second[fixed_coordinate, n2]
    return holonomy


def covariant_laplacian_spectrum(
    links: Sequence[np.ndarray],
    *,
    lattice_spacing: float = 1.0,
) -> np.ndarray:
    """Return the Hermitian nearest-neighbor covariant-Laplacian spectrum.

    The input is a square periodic ``N x N`` link field of unitary ``d x d``
    matrices. Eigenvalues carry units of ``lattice_spacing**-2``. This is a
    finite-lattice operator; continuum or bundle interpretations require
    separate construction-level evidence.
    """

    first, second = _validated_links(links)
    spacing = float(lattice_spacing)
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError(f"lattice_spacing must be finite and positive; got {spacing!r}")

    n_side = first.shape[0]
    dimension = first.shape[-1]
    size = n_side * n_side * dimension
    laplacian = np.zeros((size, size), dtype=complex)

    def flattened_index(n1: int, n2: int, component: int) -> int:
        return (n1 * n_side + n2) * dimension + component

    for n1 in range(n_side):
        for n2 in range(n_side):
            for direction, field in ((0, first), (1, second)):
                forward = (
                    ((n1 + 1) % n_side, n2)
                    if direction == 0
                    else (n1, (n2 + 1) % n_side)
                )
                link = field[n1, n2]
                for row in range(dimension):
                    site_index = flattened_index(n1, n2, row)
                    forward_index = flattened_index(forward[0], forward[1], row)
                    laplacian[site_index, site_index] += 1
                    laplacian[forward_index, forward_index] += 1
                    for column in range(dimension):
                        neighbor_index = flattened_index(
                            forward[0], forward[1], column
                        )
                        laplacian[site_index, neighbor_index] -= link[row, column]
                        laplacian[neighbor_index, site_index] -= np.conj(
                            link[row, column]
                        )

    eigenvalues = np.linalg.eigvalsh(laplacian).real / spacing**2
    return np.sort(eigenvalues)
