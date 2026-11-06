"""Connected finite-fermion two-point spectra and +i commutator response to -x B."""
from __future__ import annotations

import sympy as sp


def retarded(energies, occupations, operators, z, hbar=1):
    """chi_ab=sum_ij (f_i-f_j) A^a_ij A^b_ji/(e_j-e_i-hbar*z).

    Fourier convention is exp(+i z t), Im(z)>0. Operators need not be real
    or mutually commuting. No spin factor, implicit broadening or static term.
    """
    return sp.ImmutableMatrix(len(operators), len(operators), lambda a, b: sum(
        (occupations[i] - occupations[j]) * operators[a][i, j]
        * operators[b][j, i] / (energies[j] - energies[i] - hbar * z)
        for i in range(len(energies)) for j in range(len(energies))
        if energies[j] != energies[i]))


def noise_lines(energies, occupations, operators):
    """Energy-keyed weights of C_ab(t)=sum_D W_ab(D) exp(-i D t/hbar).

    S_ab(omega)=2pi sum_D W_ab(D) delta_Dirac(omega-D/hbar).
    Includes diagonal and all degenerate coherences at exactly D=0.
    """
    result = {}
    for i, ei in enumerate(energies):
        for j, ej in enumerate(energies):
            gap = sp.expand_log(sp.simplify(ej - ei), force=True)
            weight = sp.Matrix(len(operators), len(operators), lambda a, b:
                               occupations[i] * (1 - occupations[j])
                               * operators[a][i, j] * operators[b][j, i])
            result[gap] = result.get(gap, sp.zeros(len(operators))) + weight
    return {gap: sp.ImmutableMatrix(weight) for gap, weight in result.items()}


def symmetrized_lines(lines):
    """Weights of <{delta A_a(t),delta A_b(0)}> / 2; transpose, not adjoint."""
    return {gap: sp.ImmutableMatrix((weight + lines[-gap].T) / 2)
            for gap, weight in lines.items()}


def isolated_static(energies, occupations, operators):
    """Finite isolated retarded z->0 limit; exact zero gaps have zero commutator."""
    return retarded(energies, occupations, operators, 0)
