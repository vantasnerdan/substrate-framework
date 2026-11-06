"""Finite ordinary-3D Ewald sums in physical bohr/Ry, without force processing.

Internal kernels for qe_ionic_ewald. Cell vectors are columns, atom rows are
Cartesian. QE rgen's R-dtau sign is reversed here so pair forces point along
Ri-Rj-R. All signed reciprocal vectors are retained (no Gamma factor).
"""
from itertools import product
from math import pi, sqrt

import numpy as np
from scipy.special import erfc


def _integer_box(bounds):
    return np.asarray(list(product(*(range(-int(n), int(n) + 1) for n in bounds))), dtype=int)


def reciprocal_vectors(cell, cutoff):
    """Complete spherical physical reciprocal set, excluding G=0."""
    basis = 2 * pi * np.linalg.inv(cell).T
    bounds = np.ceil(cutoff * np.linalg.norm(cell, axis=0) / (2 * pi)).astype(int)
    miller = _integer_box(bounds)
    vectors = miller @ basis.T
    squared = np.einsum('ij,ij->i', vectors, vectors)
    keep = (squared > 0) & (squared <= cutoff * cutoff)
    return vectors[keep], squared[keep]


def _real_pair(displacement, cell, inverse, translations, radius, guard):
    fractional = inverse @ displacement
    folded = cell @ (fractional - np.rint(fractional))
    vectors = folded - translations
    radii = np.linalg.norm(vectors, axis=1)
    gap = float(np.min(np.abs(radii - radius)))
    if gap <= guard:
        raise ValueError('real cutoff crossing: the finite scalar derivative is not defined here')
    return vectors, radii, gap


def real_sum(cell, positions, charges, alpha, cutoff):
    """Pair-counted scalar and its local derivative; exact self only is excluded."""
    inverse = np.linalg.inv(cell)
    bounds = np.ceil(cutoff * np.linalg.norm(inverse, axis=1) + .5).astype(int)
    translations = _integer_box(bounds) @ cell.T
    guard = 64 * np.finfo(float).eps * max(cutoff, float(np.linalg.norm(cell)))
    energy, gap, count = 0., float('inf'), 0
    forces = np.zeros_like(positions)
    root = sqrt(alpha)
    for i in range(len(charges)):
        for j in range(i, len(charges)):
            vectors, radii, pair_gap = _real_pair(
                positions[i] - positions[j], cell, inverse, translations, cutoff, guard)
            gap = min(gap, pair_gap)
            if i != j and np.any(radii <= guard):
                raise ValueError('coincident periodic ions: Coulomb scalar is singular')
            keep = (radii <= cutoff) & (radii > 0)
            r, d = radii[keep], vectors[keep]
            count += len(r) * (1 if i == j else 2)
            screened = erfc(root * r)
            pair_charge = 2 * charges[i] * charges[j]  # QE e2=2 in Ry
            energy += pair_charge * float(np.sum(screened / r)) * (.5 if i == j else 1.)
            if i != j:
                radial = screened / r**3 + 2 * root / sqrt(pi) * np.exp(-alpha * r*r) / r**2
                pair_force = pair_charge * (radial @ d)
                forces[i] += pair_force
                forces[j] -= pair_force
    return energy, forces, count, gap


def reciprocal_sum(cell, positions, charges, alpha, cutoff):
    """Structure-factor scalar and analytic Cartesian force using one G set."""
    vectors, squared = reciprocal_vectors(cell, cutoff)
    volume = abs(float(np.linalg.det(cell)))
    weights = np.exp(-squared / (4 * alpha)) / squared
    phase = np.exp(1j * (positions @ vectors.T))
    structure = charges @ phase
    energy = 4 * pi / volume * float(np.sum(weights * np.abs(structure)**2))
    force = (8 * pi / volume * charges[:, None]
             * ((weights[None, :] * (phase * structure.conj()).imag) @ vectors))
    return energy, force, len(vectors)
