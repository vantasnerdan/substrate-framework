"""Original Gamma canonical coordinates and normalized common-functional legs.

Only a component evaluator: no mean, contact, preparation or XC replacement.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in ('0045','0042','0041','0030','0029')]
import numpy as np
from electron_ion_functional import EdgePair, evaluate_components


def coordinate_indices(c):
    return np.r_[np.arange(c.d), np.arange(2*c.d, 2*c.d+c.f)]


def vertices(c, block=None):
    """Hamiltonian derivatives and independently stored gradient contractions.

    Both returned arrays [slot,coordinate,target,source] include normalization.
    Reverse remains a contraction coefficient, NOT an inferred dagger.
    """
    block = c.q(0) if block is None else block
    return (np.concatenate((block.ionic/c.sqrt_n, c.field_modes), axis=1),
            np.concatenate((block.ionic_reverse/c.sqrt_n,
                            c.field_reverse_modes), axis=1))


def ordered_coordinates(reference, a, b):
    c = reference.c
    sq, sm = reference.transitions(0)
    indices = coordinate_indices(c)
    w = sq[a] @ c.mechanical.w0[0] @ sm[b].T
    reverse = sm[b] @ c.mechanical.w0[0] @ sq[a].T
    return -1j*w[np.ix_(indices, indices)], -1j*reverse[np.ix_(indices, indices)].T


def common_pair(c, reference, a, b, *, ionic_only=False):
    """Full128 explicit CAR sectors, all114 coordinates, NO further 1/N.

    Restricting to the first d coordinates is the exactly same ionic component
    as A41's unscaled-vertex / N edge convention, not a separate approximation.
    """
    forward, reverse = vertices(c)
    dg, dl = ordered_coordinates(reference, a, b)
    if ionic_only:
        forward, reverse = forward[:, :c.d], reverse[:, :c.d]
        dg, dl = dg[:c.d, :c.d], dl[:c.d, :c.d]
    gg, gl = reference.electron(a, b)
    bg, bl = reference.electron(b, a)
    total = {}
    for source in range(c.ns):
        target = c.target[0, source]
        left, right = forward[source], reverse[source].transpose(0, 2, 1)
        pair = EdgePair(left, right, left, right, gg[source], gl[source],
                        gg[target], gl[target], bg[source], bl[source], dg, dl)
        values = evaluate_components(pair, n=1, t=reference.times[a], s=reference.times[b])
        for name, value in values.items():
            if name.startswith('sigma'):
                if name not in total:
                    total[name] = np.zeros((c.ns, c.m, c.m), complex)
                total[name][target] += value
            else:
                total[name] = total.get(name, 0) + value
    return total
