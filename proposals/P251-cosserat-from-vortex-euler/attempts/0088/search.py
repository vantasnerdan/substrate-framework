"""Hypothesis generation only; print exact rational Rayleigh witnesses."""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

from itertools import product

import numpy as np
from scipy.linalg import eigh


def modes(bound):
    return [k for k in product(range(-bound, bound+1), repeat=3)
            if k != (0, 0, 0) and next(x for x in k if x) > 0]


def basis(bound):
    result = []
    for k in modes(bound):
        vector = np.array(k)
        axis = int(np.argmin(np.abs(vector)))
        p = np.cross(vector, np.eye(3, dtype=int)[axis])
        second = np.cross(vector, p)
        for polarization in (p, second):
            for sine in (False, True):
                result.append((k, polarization, sine))
    return result


def background(n):
    # Twice the Fourier coefficients of the equal-amplitude ABC field.
    return {(n, 0, 0): np.array([0, -1j, 1]),
            (-n, 0, 0): np.array([0, 1j, 1]),
            (0, n, 0): np.array([1, 0, -1j]),
            (0, -n, 0): np.array([1, 0, 1j]),
            (0, 0, n): np.array([-1j, 1, 0]),
            (0, 0, -n): np.array([1j, 1, 0])}


def run(n, bound, supplied=None, name="ABC"):
    fields = basis(bound)
    u = background(n) if supplied is None else supplied
    for wave, value in u.items():
        assert np.array_equal(1j*np.cross(wave, value), n*value)
    forces = []
    frequencies = set()
    for wave, polarization, sine in fields:
        field = {}
        for sign in (-1, 1):
            phase = -sign*1j if sine else 1
            for carrier, value in u.items():
                target = tuple(sign*np.array(wave)+np.array(carrier))
                field[target] = field.get(target, np.zeros(3, complex)) \
                    + phase*np.cross(polarization, value)
        frequencies.update(field)
        forces.append(field)
    frequencies = sorted(frequencies)
    index = {wave: i for i, wave in enumerate(frequencies)}
    f = np.zeros((len(frequencies), 3, len(fields)), complex)
    w = np.zeros_like(f)
    for column, force in enumerate(forces):
        for wave, value in force.items():
            f[index[wave], :, column] = value
            w[index[wave], :, column] = 1j*np.cross(wave, value)
    f = f.reshape(-1, len(fields))
    w = w.reshape(-1, len(fields))
    matrix = (n*f.conj().T@w-w.conj().T@w).real
    assert np.array_equal(matrix, matrix.T)
    assert np.array_equal(matrix, np.rint(matrix))
    masses = np.array([np.dot(p, p)/2 for _, p, _ in fields])
    normalized = matrix/np.sqrt(masses[:, None]*masses[None, :])/16
    values, vectors = eigh(normalized, subset_by_index=(len(fields)-2, len(fields)-1))
    print("SEARCH_ONLY", name, n, bound, "dimension", len(fields), "top", values, flush=True)
    if values[-1] <= 1e-8:
        return
    direction = vectors[:, -1]/np.sqrt(masses)
    direction /= np.max(np.abs(direction))
    integers = np.rint(100*direction).astype(int)
    # Python integer accumulation is exact, including cancellation.
    exact = sum(int(integers[i])*int(matrix[i, j])*int(integers[j])
                for i in range(len(fields)) for j in range(len(fields)))
    mass_twice = sum(int(q)**2*int(np.dot(p, p))
                     for q, (_, p, _) in zip(integers, fields, strict=True))
    print("RATIONAL_CANDIDATE K numerator /16:", exact,
          "mass numerator /2:", mass_twice, flush=True)
    if exact > 0:
        for q, (wave, p, sine) in zip(integers, fields, strict=True):
            if q:
                print(int(q), wave, tuple(map(int, p)), "sin" if sine else "cos")


if __name__ == "__main__":
    for carrier in (1, 2, 3):
        for cutoff in (1, 2, 3):
            run(carrier, cutoff)
