"""Finite mean-field force response; instantaneous interaction, dynamic screening.

The joint positive susceptibility uses h'= -sum u_a B_a + sum D_s U_st delta d_t.
B may be nonlocal Hermitian operators; D labels declared density channels.
This algebra does not identify the supplied kernel or vertices with material PBE.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve


def screening_fields(joint, n_external, interaction):
    """Return U delta d / u from the selfconsistent induced-density equation."""
    chi = np.asarray(joint, dtype=complex)
    interaction = np.asarray(interaction, dtype=complex)
    dd = chi[n_external:, n_external:]
    db = chi[n_external:, :n_external]
    return interaction @ solve(np.eye(len(dd)) + dd @ interaction, db)


def screened_response(joint, n_external, interaction):
    """Eliminate induced density at the SAME frequency as the supplied chi."""
    chi = np.asarray(joint, dtype=complex)
    return (chi[:n_external, :n_external]
            - chi[:n_external, n_external:]
            @ screening_fields(chi, n_external, interaction))


def static_screened_mixed(dynamic_joint, static_joint, n_external, interaction):
    """Bare/dynamically responding left side; frozen static screened right side.

    This is an explicitly different approximation, not selfconsistent screening
    of the same functional at finite frequency. Its static limit matches only
    when dynamic_joint is the same constrained static response used on the right.
    """
    chi = np.asarray(dynamic_joint, dtype=complex)
    return (chi[:n_external, :n_external]
            - chi[:n_external, n_external:]
            @ screening_fields(static_joint, n_external, interaction))


def linearized_density_response(energies, occupations, external, densities,
                                interaction, z):
    """Independent direct i*d(gamma)/dt=[h(gamma,u),gamma] linear solve.

    Returns Tr[B_a delta gamma_b] without using a precontracted susceptibility.
    Diagonal populations are stationary in this isolated number-conserving model.
    """
    energies = np.asarray(energies, dtype=float)
    occupations = np.asarray(occupations, dtype=float)
    external = np.asarray(external, dtype=complex)
    densities = np.asarray(densities, dtype=complex)
    interaction = np.asarray(interaction, dtype=complex)
    size = len(energies)
    shape = size * size
    equation = np.zeros((shape, shape), dtype=complex)
    rhs = np.zeros((shape, len(external)), dtype=complex)
    for i in range(size):
        for j in range(size):
            row = i * size + j
            weight = occupations[j] - occupations[i]
            equation[row, row] = energies[i] - energies[j] - z
            for k in range(size):
                for l in range(size):
                    equation[row, k * size + l] += weight * sum(
                        densities[s, i, j] * interaction[s, t]
                        * densities[t, l, k]
                        for s in range(len(densities))
                        for t in range(len(densities)))
            rhs[row] = weight * external[:, i, j]
    solutions = solve(equation, rhs)
    return np.array([[np.trace(operator @ solutions[:, b].reshape(size, size))
                      for b in range(len(external))] for operator in external])
