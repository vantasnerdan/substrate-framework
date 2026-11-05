"""Unrestricted spatial/spin/momentum coherences for the same A41 diagram.

This is a coordinate adapter, not a propagator or second approximation. Electron
axes may include ALL band, spin and k indices without block-diagonal assumptions.
Mode axes are independent Hermitian Cartesian coordinates (e.g. all original
real-space cell/ion coordinates, or a complete real Fourier transform retaining
both sine/cosine sectors). Supply the actual Hamiltonian derivative matrices in
those coordinates with ALL normalization already included: there is no extra N.
No momentum/spin symmetry is assumed and no paired mode is discarded. Complex
traveling-wave u_q are NOT independent Hermitian coordinates: their caller must
retain the full paired-coordinate metric when converting to/from this API.
The bounded material consumer uses the equivalent original directed-q edge sum,
not a prohibitively large dense all-cell vertex tensor.
"""
from dataclasses import dataclass

import numpy as np

from electron_ion_functional import EdgePair, evaluate, jvp, transpose


@dataclass(frozen=True)
class GlobalPair:
    vertex_t: np.ndarray
    vertex_s: np.ndarray
    g_greater: np.ndarray
    g_lesser: np.ndarray
    back_greater: np.ndarray
    back_lesser: np.ndarray
    d_greater: np.ndarray
    d_lesser: np.ndarray


def edge_operands(pair):
    """All equivalent appearances remain shared, not independent new variables."""
    return EdgePair(pair.vertex_t,pair.vertex_t,pair.vertex_s,pair.vertex_s,
                    pair.g_greater,pair.g_lesser,pair.g_greater,pair.g_lesser,
                    pair.back_greater,pair.back_lesser,pair.d_greater,pair.d_lesser)


def evaluate_global(pair, *, t, s):
    return evaluate(edge_operands(pair),n=1,t=t,s=s)


def jvp_global(pair, direction, *, t, s):
    return jvp(edge_operands(pair),edge_operands(direction),n=1,t=t,s=s)


def transpose_global(pair, dual, *, t, s):
    """Scatter-add EVERY shared occurrence under the nonconjugating pairing."""
    cotangent = transpose(edge_operands(pair),dual,n=1,t=t,s=s)
    return GlobalPair(cotangent.left_t+cotangent.right_t,
        cotangent.left_s+cotangent.right_s,
        cotangent.gs_greater+cotangent.gt_greater,
        cotangent.gs_lesser+cotangent.gt_lesser,
        cotangent.back_greater,cotangent.back_lesser,
        cotangent.d_greater,cotangent.d_lesser)
