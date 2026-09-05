"""Exact finite-field minors for the frozen analytic compact-velocity probe."""

from itertools import product
from math import comb

import numpy as np

ORDER = 6
PRIME = 101


def index_set(order):
    return sorted((a for a in product(range(order + 1), repeat=3) if sum(a) <= order),
                  key=lambda a: (sum(a), a))


def rat(num, den=1):
    return num * pow(den, -1, PRIME) % PRIME


def cross(a, b):
    return [(a[1] * b[2] - a[2] * b[1]) % PRIME,
            (a[2] * b[0] - a[0] * b[2]) % PRIME,
            (a[0] * b[1] - a[1] * b[0]) % PRIME]


def modular_rank(matrix):
    """Exact Gaussian elimination over F_101; products fit signed int64."""
    matrix = np.array(matrix, dtype=np.int64, copy=True) % PRIME
    pivot = 0
    for col in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[pivot:, col])
        if not len(candidates):
            continue
        chosen = pivot + int(candidates[0])
        matrix[[pivot, chosen]] = matrix[[chosen, pivot]]
        matrix[pivot] = matrix[pivot] * pow(int(matrix[pivot, col]), -1, PRIME) % PRIME
        if pivot + 1 < matrix.shape[0]:
            factor = matrix[pivot + 1:, col].copy()
            matrix[pivot + 1:] = (matrix[pivot + 1:]
                                   - factor[:, None] * matrix[pivot]) % PRIME
        pivot += 1
        if pivot == matrix.shape[0]:
            break
    return pivot


waves = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
         (rat(3, 5), rat(4, 5), 0), (0, rat(3, 5), rat(4, 5)),
         (rat(4, 5), 0, rat(3, 5)),
         (rat(1, 3), rat(2, 3), rat(2, 3)),
         (rat(2, 3), rat(-2, 3), rat(1, 3))]
amplitudes = []
for number, wave in enumerate(waves):
    reference = [number + 2, 2 * number + 1, 3 - number]
    dot = sum(reference[i] * wave[i] for i in range(3)) % PRIME
    cosine = [(reference[i] - dot * wave[i]) % PRIME for i in range(3)]
    sine = [(-v) % PRIME for v in cross(wave, cosine)]
    amplitudes.append((cosine, sine))

indices = index_set(ORDER)
outputs = [a for a in index_set(ORDER + 1) if sum(a)]
output_index = {a: i for i, a in enumerate(outputs)}
jets = {}
for alpha in index_set(ORDER):
    degree = sum(alpha)
    value = [0, 0, 0]
    for wave, (cosine, sine) in zip(waves, amplitudes, strict=True):
        monomial = np.prod([pow(wave[i], alpha[i], PRIME) for i in range(3)]) % PRIME
        for j in range(3):
            trigonometric = cosine[j] * [1, 0, -1, 0][degree % 4]
            trigonometric += sine[j] * [0, 1, 0, -1][degree % 4]
            value[j] = (value[j] + int(monomial) * trigonometric) % PRIME
    jets[alpha] = value

matrix = np.zeros((2 * len(outputs), 3 * len(indices)), dtype=np.int64)
angular = np.zeros((3, 3 * len(indices)), dtype=np.int64)
for j in range(3):
    for column, alpha in enumerate(indices):
        col = j * len(indices) + column
        beta = list(alpha)
        beta[j] += 1
        matrix[output_index[tuple(beta)], col] += 1
        for gamma in product(*(range(a + 1) for a in alpha)):
            factor = (-1)**sum(gamma)
            for a, g in zip(alpha, gamma, strict=True):
                factor *= comb(a, g)
            omega_jet = jets[gamma]
            basis = [int(k == j) for k in range(3)]
            force_jet = cross(basis, omega_jet)
            for i in range(3):
                beta = [alpha[k] - gamma[k] + int(k == i) for k in range(3)]
                matrix[len(outputs) + output_index[tuple(beta)], col] += factor * force_jet[i]
        for axis in range(3):
            value = 0
            if j == axis:
                for m in range(3):
                    if alpha[m]:
                        lower = list(alpha)
                        lower[m] -= 1
                        value += alpha[m] * jets[tuple(lower)][m]
            if alpha[j]:
                lower = list(alpha)
                lower[j] -= 1
                value -= alpha[j] * jets[tuple(lower)][axis]
            angular[axis, col] = (-1)**sum(alpha) * value
matrix %= PRIME
angular %= PRIME
rank = modular_rank(matrix)
print(f"N={ORDER}; prime={PRIME}; unknowns={matrix.shape[1]}; equations={matrix.shape[0]}")
print(f"constraint modular rank={rank}; proved row-count upper bound={matrix.shape[0]}")
for axis in range(3):
    augmented = modular_rank(np.vstack([matrix, angular[axis]]))
    print(f"spin axis {axis}: augmented modular rank={augmented}")
print(f"all-spin augmented modular rank={modular_rank(np.vstack([matrix, angular]))}")
if rank != matrix.shape[0]:
    print("SCOPE: unaugmented modular rank below row-count bound; no angular verdict from minors alone")
