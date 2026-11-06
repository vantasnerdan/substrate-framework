"""Finite CAR representation and physical Gibbs ensembles, without smearing aliases."""
from __future__ import annotations

import sympy as sp


def annihilators(orbitals: int) -> tuple[sp.ImmutableMatrix, ...]:
    """Occupation bit i labels orbital i; lower occupied bits give the CAR sign."""
    result = []
    for i in range(orbitals):
        matrix = sp.zeros(1 << orbitals)
        for state in range(1 << orbitals):
            if state & (1 << i):
                matrix[state ^ (1 << i), state] = (-1) ** (
                    (state & ((1 << i) - 1)).bit_count())
        result.append(sp.ImmutableMatrix(matrix))
    return tuple(result)


def second_quantize(one_body) -> sp.ImmutableMatrix:
    """Return dGamma(A)=sum A_ij c_i† c_j in the full finite Fock space."""
    matrix = sp.Matrix(one_body)
    if matrix.rows != matrix.cols:
        raise ValueError("one-body operator must be square")
    c = annihilators(matrix.rows)
    result = sp.zeros(1 << matrix.rows)
    for i in range(matrix.rows):
        for j in range(matrix.rows):
            result += matrix[i, j] * c[i].H * c[j]
    return sp.ImmutableMatrix(result)


def fermi_dirac(energies, beta, mu=0) -> tuple:
    """Physical noninteracting grandcanonical occupations, beta=1/(k_B T)."""
    return tuple(sp.simplify(1 / (1 + sp.exp(beta * (e - mu))))
                 for e in energies)


def fock_gibbs(energies, beta, mu=0, number=None):
    """Exact diagonal Gibbs density and partition sum; number selects finite canonical."""
    dimension = 1 << len(energies)
    many_energies = tuple(sum(e for i, e in enumerate(energies) if s & (1 << i))
                          for s in range(dimension))
    weights = [sp.simplify(sp.exp(-beta * (e - mu * s.bit_count())))
               if number is None or s.bit_count() == number else sp.S.Zero
               for s, e in enumerate(many_energies)]
    partition = sp.simplify(sum(weights))
    rho = sp.diag(*(sp.simplify(w / partition) for w in weights))
    return sp.ImmutableMatrix(rho), partition, many_energies
