import pytest
import sympy as sp

from substrate_framework.euler_neutral_cell import (
    octahedral_rank2_orbit_sum,
    octahedral_rank3_orbit_sum,
    octahedral_vector_orbit_sum,
    proper_octahedral_rotations,
)


def test_proper_octahedral_group_has_24_exact_rotations():
    rotations = proper_octahedral_rotations()
    assert len(rotations) == 24
    assert len(set(rotations)) == 24
    assert all(rotation.det() == 1 for rotation in rotations)
    assert all(rotation.T * rotation == sp.eye(3) for rotation in rotations)


def test_vector_and_rank_two_moment_orbits_remove_low_anisotropy():
    a, b, c, d, e, f = sp.symbols("a b c d e f", real=True)
    vector = sp.Matrix([a, b, c])
    matrix = sp.Matrix([[a, d, e], [d, b, f], [e, f, c]])
    assert octahedral_vector_orbit_sum(vector) == sp.zeros(3, 1)
    assert octahedral_rank2_orbit_sum(matrix) == 8 * sp.trace(matrix) * sp.eye(3)
    antisymmetric = sp.Matrix([[0, a, b], [-a, 0, c], [-b, -c, 0]])
    assert octahedral_rank2_orbit_sum(antisymmetric) == sp.zeros(3)


def test_symmetric_rank_three_orbit_has_no_octahedral_invariant():
    values = sp.MutableDenseNDimArray.zeros(3, 3, 3)
    counter = 1
    for i in range(3):
        for j in range(i, 3):
            for k in range(3):
                values[i, j, k] = counter
                values[j, i, k] = counter
                counter += 1
    assert octahedral_rank3_orbit_sum(values) == sp.ImmutableDenseNDimArray.zeros(3, 3, 3)


def test_invalid_moment_domains_are_rejected():
    with pytest.raises(ValueError):
        octahedral_vector_orbit_sum([1, 2])
    with pytest.raises(ValueError):
        octahedral_rank2_orbit_sum([[0, 1], [0, 0]])
    bad = sp.MutableDenseNDimArray.zeros(3, 3, 3)
    bad[0, 1, 2] = 1
    with pytest.raises(ValueError):
        octahedral_rank3_orbit_sum(bad)
