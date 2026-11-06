"""Source-parameterized arithmetic for a fixed-cell raw-force boundary.

No engine, material constants, publication schema, or acceptance gate is implicit.
"""
from itertools import product

import numpy as np
import sympy as sp
from scipy.linalg import norm


def atomic_units_2006(c, mu0, h, charge, electron_mass):
    """Derive Angstrom/Bohr and eV/Ry from SI defining relations."""
    eps0 = 1 / (mu0 * c**2)
    hbar = h / (2 * sp.pi)
    coulomb = charge**2 / (4 * sp.pi * eps0)
    bohr_m = hbar**2 / (electron_mass * coulomb)
    rydberg_joule = electron_mass * coulomb**2 / (2 * hbar**2)
    return sp.simplify(bohr_m * 10**10), sp.simplify(rydberg_joule / charge)


def squared_norms_decimal(vectors):
    """Exact decimal squares retain every sign/component before reduction."""
    return [sum(sp.Rational(str(x))**2 for x in row) for row in vectors]


def raw_gate(vectors, gate):
    squares = squared_norms_decimal(vectors)
    return bool(squares and max(squares) <= sp.Rational(str(gate))**2)


def force_rounding_upper(vectors, component_half_quantum):
    """Triangle inequality encloses norms of rounded three-component forces."""
    return np.array([norm(row) for row in vectors]) + np.sqrt(3) * component_half_quantum


def com_project(forces, masses):
    """Orthogonal force projection onto the mass-defined COM tangent space."""
    f, m = np.asarray(forces, float), np.asarray(masses, float)
    return f - m[:, None] * (m @ f / (m @ m))


def cubic_pair_distance(first, second, length):
    """Exact minimum-image distance for an orthogonal equal-edge periodic cell."""
    delta = [sp.Rational(str(a)) - sp.Rational(str(b)) for a, b in zip(first, second)]
    length = sp.Rational(str(length))
    wrapped = [d - length * sp.floor(d / length + sp.Rational(1, 2)) for d in delta]
    return sp.sqrt(sum(d*d for d in wrapped))


def image_enumeration_distance(first, second, cell):
    """Independent finite-image Euclidean search for points inside one cell."""
    first, second, cell = np.asarray(first), np.asarray(second), np.asarray(cell)
    return min(norm(first-second + np.asarray(n) @ cell) for n in product(range(-1, 2), repeat=3))
