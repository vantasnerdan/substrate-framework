"""Exact octahedral moment projectors for neutral Euler-cell candidates.

The projectors encode physical rotation-orbit cancellation.  They do not
construct an invariant Euler mode, a lattice band, or a particle.
"""

from itertools import permutations, product

import sympy as sp


def proper_octahedral_rotations():
    """Return the 24 signed permutation matrices with determinant ``+1``."""

    rotations = []
    for permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            matrix = sp.zeros(3)
            for row, column in enumerate(permutation):
                matrix[row, column] = signs[row]
            if matrix.det() == 1:
                rotations.append(sp.ImmutableMatrix(matrix))
    return tuple(rotations)


def octahedral_vector_orbit_sum(vector):
    """Sum a polar vector over the proper octahedral rotation orbit."""

    value = sp.ImmutableMatrix(vector)
    if value.shape != (3, 1):
        raise ValueError("vector must have three components")
    total = sum(
        (rotation * value for rotation in proper_octahedral_rotations()),
        sp.zeros(3, 1),
    )
    return sp.ImmutableMatrix(total)


def octahedral_rank2_orbit_sum(tensor):
    """Sum a rank-two tensor over the rotation orbit."""

    value = sp.ImmutableMatrix(tensor)
    if value.shape != (3, 3):
        raise ValueError("tensor must be 3 by 3")
    total = sp.zeros(3)
    for rotation in proper_octahedral_rotations():
        total += rotation * value * rotation.T
    return sp.ImmutableMatrix(total)


def octahedral_rank3_orbit_sum(tensor):
    """Sum a rank-three tensor symmetric in its first two indices."""

    value = sp.ImmutableDenseNDimArray(tensor)
    if value.shape != (3, 3, 3):
        raise ValueError("tensor must have shape 3 by 3 by 3")
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if value[i, j, k] != value[j, i, k]:
                    raise ValueError("tensor must be symmetric in its first two indices")
    total = sp.MutableDenseNDimArray.zeros(3, 3, 3)
    for rotation in proper_octahedral_rotations():
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    total[i, j, k] += sum(
                        rotation[i, a] * rotation[j, b] * rotation[k, c] * value[a, b, c]
                        for a in range(3)
                        for b in range(3)
                        for c in range(3)
                    )
    return sp.ImmutableDenseNDimArray(total)
