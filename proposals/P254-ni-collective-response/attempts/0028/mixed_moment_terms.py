"""Borrowed original-coefficient terms shared by mixed moment contractions."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MomentTerm:
    group: str
    left: np.ndarray
    D: np.ndarray
    right: np.ndarray
    weight: float


def moment_terms(coupling, groups):
    """Keep original columns, coefficients, ordering and per-source weights.

    Repeated immutable arrays retain one normalized object so request-local
    arithmetic caches can share factors. Nothing is copied, inverted, repaired
    or made Hermitian beyond the existing bridge's required matrix conversion.
    """
    matrices, coefficients, result = {}, {}, []
    bridge = coupling.bridge
    for group, block in groups.items():
        for dyad in block.dyads:
            values = []
            for frame, rows in ((dyad.left, coupling.edge.target_dimension),
                                (dyad.right, coupling.edge.source_dimension)):
                key = (id(frame), rows)
                if key not in matrices:
                    matrices[key] = bridge._matrix(frame, rows)
                values.append(matrices[key])
            left, right = values
            key = (id(dyad.D), left.shape[1], right.shape[1])
            if key not in coefficients:
                coefficient = np.asarray(dyad.D)
                if (coefficient.shape != (left.shape[1], right.shape[1])
                        or not np.isfinite(coefficient).all()):
                    raise ValueError('finite original rectangular dyad coefficient required')
                coefficients[key] = coefficient
            result.append(MomentTerm(group, left, coefficients[key], right,
                                     bridge._weight(dyad.weight)))
    return result
