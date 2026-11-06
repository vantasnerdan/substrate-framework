"""Attempt0013 inputs frozen before execution; arbitrary dimensionless units.

Four total orbitals: two up and two down, not four orbitals per spin.
No parameters are extracted from Ni, SG15, PBE, or observed pilot results.
"""
from __future__ import annotations

import numpy as np
import sympy as sp


def _blocks(diagonal, up=0, down=0):
    matrix = np.diag(np.asarray(diagonal, dtype=float))
    matrix[0, 1] = matrix[1, 0] = up
    matrix[2, 3] = matrix[3, 2] = down
    return matrix


def symbolic_reference():
    """Exact reference for symbolic response helpers; no floating zero gaps."""
    return ([sp.S.Zero, sp.log(2), sp.log(3), sp.log(5)],
            [sp.Rational(1, 2), sp.Rational(1, 3),
             sp.Rational(1, 4), sp.Rational(1, 6)])


def make_fixture():
    """Return fresh arrays; B[a], D[s], W[a,b] act on the same 4D space.

    d_s=Tr(D_s gamma), Z=diag(+,+,-,-), h(u)=h0-u.B+u.W.u/2.
    The reference has N0=5/4 and M0=5/12 at beta=1, mu=0.
    U is positive definite with a deliberately nonzero cross-spin channel.
    """
    energies, occupations = symbolic_reference()
    h0 = np.diag(np.asarray(energies, dtype=float))
    gamma0 = np.diag(np.asarray(occupations, dtype=float))
    B = np.array([
        _blocks([1/2, -1/3, 1/4, 1/5], 1/4, 1/3),
        _blocks([-1/4, 2/5, 1/3, -1/2], -1/6, 1/5),
    ])
    D = np.array([
        _blocks([1, 2/3, 0, 0], 1/5, 0),
        _blocks([0, 0, 3/4, 1/2], 0, 1/4),
    ])
    W11 = _blocks([1/3, 1/5, -1/7, 1/4], 1/11, -1/13)
    W12 = _blocks([1/7, -1/8, 1/9, 1/10], 1/12, 1/14)
    W22 = _blocks([-1/6, 1/4, 1/5, -1/9], -1/10, 1/8)
    return {
        "h0": h0, "gamma0": gamma0, "beta": 1., "mu0": 0.,
        "N0": 5/4, "M0": 5/12,
        "d0": np.einsum("sij,ji->s", D, gamma0),
        "B": B, "D": D, "W": np.array([[W11, W12], [W12, W22]]),
        "U": np.array([[3/4, 1/4], [1/4, 2/3]]),
        "Kion": np.array([[2, 1/5], [1/5, 3/2]]),
        "Z": np.diag([1., 1., -1., -1.]), "spin": (0, 0, 1, 1),
    }
