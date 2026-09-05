"""Append-only exact diagnostic of the three finite-jet row dependencies."""

import numpy as np

from jet_probe import PRIME, matrix, outputs

work = matrix.T.copy()
pivots = []
row = 0
for col in range(work.shape[1]):
    candidates = np.flatnonzero(work[row:, col])
    if not len(candidates):
        continue
    chosen = row + int(candidates[0])
    work[[row, chosen]] = work[[chosen, row]]
    work[row] = work[row] * pow(int(work[row, col]), -1, PRIME) % PRIME
    factors = work[:, col].copy()
    factors[row] = 0
    work = (work - factors[:, None] * work[row]) % PRIME
    pivots.append(col)
    row += 1
    if row == work.shape[0]:
        break
free = [col for col in range(work.shape[1]) if col not in pivots]
for col in free:
    null = np.zeros(work.shape[1], dtype=np.int64)
    null[col] = 1
    for index, pivot in enumerate(pivots):
        null[pivot] = -work[index, col] % PRIME
    assert np.all(matrix.T @ null % PRIME == 0)
    support = {}
    for idx in np.flatnonzero(null):
        component = int(idx) // len(outputs)
        degree = sum(outputs[int(idx) % len(outputs)])
        support[component, degree] = support.get((component, degree), 0) + 1
    print(f"left-kernel free coordinate {col}; support by (component,degree)={support}")
