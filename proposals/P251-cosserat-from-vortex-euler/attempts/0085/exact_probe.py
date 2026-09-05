"""Characteristic-zero continuation of the frozen finite-field probe.

The x->15x, omega->225omega rescaling makes all wavevectors and amplitudes
integer and preserves all constraint/angular row ranks. No new prototype
is selected from the first comparison.
"""

from itertools import product
from math import comb, prod

from flint import fmpz_mat

ORDER = 6


def index_set(order):
    return sorted((a for a in product(range(order + 1), repeat=3) if sum(a) <= order),
                  key=lambda a: (sum(a), a))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


waves = [(15, 0, 0), (0, 15, 0), (0, 0, 15), (9, 12, 0),
         (0, 9, 12), (12, 0, 9), (5, 10, 10), (10, -10, 5)]
amplitudes = []
for number, wave in enumerate(waves):
    reference = [number + 2, 2 * number + 1, 3 - number]
    dot = sum(reference[i] * wave[i] for i in range(3))
    cosine = [225 * reference[i] - dot * wave[i] for i in range(3)]
    raw_sine = cross(wave, cosine)
    assert all(v % 15 == 0 for v in raw_sine)
    sine = [-v // 15 for v in raw_sine]
    amplitudes.append((cosine, sine))

indices = index_set(ORDER)
outputs = [a for a in index_set(ORDER + 1) if sum(a)]
output_index = {a: i for i, a in enumerate(outputs)}
jets = {}
for alpha in index_set(ORDER):
    degree = sum(alpha)
    value = [0, 0, 0]
    for wave, (cosine, sine) in zip(waves, amplitudes, strict=True):
        monomial = prod(wave[i]**alpha[i] for i in range(3))
        for j in range(3):
            trigonometric = cosine[j] * [1, 0, -1, 0][degree % 4]
            trigonometric += sine[j] * [0, 1, 0, -1][degree % 4]
            value[j] += monomial * trigonometric
    jets[alpha] = value

matrix = [[0] * (3 * len(indices)) for _ in range(2 * len(outputs))]
angular = [[0] * (3 * len(indices)) for _ in range(3)]
for j in range(3):
    for column, alpha in enumerate(indices):
        col = j * len(indices) + column
        beta = list(alpha)
        beta[j] += 1
        matrix[output_index[tuple(beta)]][col] += 1
        for gamma in product(*(range(a + 1) for a in alpha)):
            factor = (-1)**sum(gamma)
            factor *= prod(comb(a, g) for a, g in zip(alpha, gamma, strict=True))
            force_jet = cross([int(k == j) for k in range(3)], jets[gamma])
            for i in range(3):
                beta = tuple(alpha[k] - gamma[k] + int(k == i) for k in range(3))
                matrix[len(outputs) + output_index[beta]][col] += factor * force_jet[i]
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
            angular[axis][col] = (-1)**sum(alpha) * value

print(f"N={ORDER}; exact integer prototype; lambda=15; columns={len(matrix[0])}", flush=True)
rank = fmpz_mat(matrix).rank()
print(f"constraint exact rational rank={rank}", flush=True)
for axis in range(3):
    augmented = fmpz_mat(matrix + [angular[axis]]).rank()
    print(f"spin axis {axis}: augmented exact rational rank={augmented}", flush=True)
print(f"all-spin augmented exact rational rank={fmpz_mat(matrix + angular).rank()}", flush=True)
