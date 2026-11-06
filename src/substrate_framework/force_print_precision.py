"""Exact component-box arithmetic for declared force-print precision.

No material, unit conversion or numerical acceptance threshold is implicit.
"""
from itertools import product

import sympy as sp


def component_box_max_squared(center, half_quantum):
    """Maximum squared Euclidean norm for independent component intervals."""
    radius = sp.Rational(str(half_quantum))
    return sum((abs(sp.Rational(str(value))) + radius) ** 2 for value in center)


def corner_squared_norms(center, half_quantum):
    """Independently enumerate vertices instead of selecting their signs."""
    values = tuple(sp.Rational(str(value)) for value in center)
    radius = sp.Rational(str(half_quantum))
    return tuple(sum((value + sign * radius) ** 2
                     for value, sign in zip(values, signs))
                 for signs in product((-1, 1), repeat=len(values)))


def barycentric_box_identity():
    """Squared-norm interpolation gap on a three-component box.

    The caller checks the polynomial identities and the declared 0<=t<=1
    domain. Its nonnegative variance proves interior enclosure by vertices.
    """
    center = sp.symbols('f0:3', real=True)
    parameters = sp.symbols('t0:3', real=True)
    radius = sp.Symbol('h', nonnegative=True)
    total_weight, weighted_norm = 0, 0
    weighted_position = [0, 0, 0]
    for bits in product((0, 1), repeat=3):
        weight = sp.prod(t if bit else 1 - t
                         for t, bit in zip(parameters, bits))
        corner = [f + radius * (2 * bit - 1) for f, bit in zip(center, bits)]
        total_weight += weight
        weighted_norm += weight * sum(value ** 2 for value in corner)
        for axis in range(3):
            weighted_position[axis] += weight * corner[axis]
    point = [f + radius * (2 * t - 1) for f, t in zip(center, parameters)]
    return {
        'weight_sum': sp.expand(total_weight),
        'position_gaps': tuple(sp.expand(a - b) for a, b in zip(weighted_position, point)),
        'variance': sp.expand(weighted_norm - sum(value ** 2 for value in point)),
        'expected_variance': 4 * radius ** 2 * sum(t * (1 - t) for t in parameters),
        'parameters': parameters,
        'radius': radius,
    }
