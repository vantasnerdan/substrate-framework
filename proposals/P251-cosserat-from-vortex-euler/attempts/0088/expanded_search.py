"""Different exact Beltrami wave geometry after the ABC search."""

from search import run

import numpy as np


def helical(waves, eigenvalue):
    field = {}
    for index, wave in enumerate(waves):
        k = np.array(wave)
        p = np.cross(k, np.eye(3, dtype=int)[np.argmin(abs(k))])
        value = eigenvalue*p+1j*np.cross(k, p)
        phase = (1, 1j, -1, -1j)[index % 4]
        field[wave] = phase*value
        field[tuple(-k)] = (phase*value).conjugate()
    return field


if __name__ == "__main__":
    candidates = {
        "close_pair": [(4, 3, 0), (4, -3, 0)],
        "tilted_triad": [(4, 3, 0), (4, -3, 0), (0, 4, 3)],
        "six_waves": [(4, 3, 0), (4, -3, 0), (0, 4, 3),
                      (0, 4, -3), (3, 0, 4), (-3, 0, 4)],
    }
    for name, directions in candidates.items():
        for cutoff in (2, 3):
            run(5, cutoff, helical(directions, 5), name)
