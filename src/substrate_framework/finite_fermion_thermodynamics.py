"""Isothermal one-body Gibbs curvature and mean-number Legendre projection."""
from __future__ import annotations

import sympy as sp

from .finite_fermion_response import isolated_static, noise_lines


def grandcanonical_curvature(energies, occupations, operators, beta):
    """-Omega_xx for h(x)=h-sum x_a A_a, at fixed beta and mu.

    The zero-energy Kubo-Mori divided difference is beta f(1-f), including
    all matrix elements within degenerate blocks, not merely diagonals.
    """
    static = isolated_static(energies, occupations, operators)
    populations = noise_lines(energies, occupations, operators)[sp.S.Zero]
    return sp.ImmutableMatrix(static + beta * populations)


def mean_number_curvature(energies, occupations, operators, beta):
    """-F_xx at fixed mean N, F=Omega+mu*N, NOT the finite canonical ensemble."""
    all_operators = tuple(operators) + (sp.eye(len(energies)),)
    full = grandcanonical_curvature(energies, occupations, all_operators, beta)
    n = len(operators)
    if full[n, n] == 0:
        raise ValueError("mean-number Legendre projection needs nonzero compressibility")
    return sp.ImmutableMatrix(full[:n, :n] - full[:n, n:] * full[n:, :n] / full[n, n])


def calibrated_inverse(stiffness, response, calibration_curvature, mass, z):
    """Quadratic action D(z)=K_cal-M z²-C(z)+C_cal.

    This is an algebraic definition, not permission to set C_cal=chi_R(0).
    C_cal must be the curvature already included in K_cal for the same
    ensemble/action. A GC free-energy Hessian includes occupation covariance;
    an isolated positive-gap bath does not. Population dynamics, if needed,
    must be specified independently rather than invented by subtraction.
    """
    return sp.ImmutableMatrix(sp.Matrix(stiffness) - sp.Matrix(mass) * z**2
                              - sp.Matrix(response) + sp.Matrix(calibration_curvature))
