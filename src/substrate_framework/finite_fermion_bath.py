"""Positive-gap electron-hole to complex Gaussian oscillator two-point map only."""
from __future__ import annotations

import sympy as sp


def oscillator_map(energies, occupations, operators):
    """Return (Delta, g) for X_a=g_a b+g_a* b†, H=Delta b†b.

    Energies must be ordered and occupations passive. Exactly degenerate
    pairs are deliberately absent: their finite static occupation covariance
    is not a thermal oscillator at Delta=0.
    """
    modes = []
    for i, ei in enumerate(energies):
        for j in range(i + 1, len(energies)):
            gap = sp.expand_log(sp.simplify(energies[j] - ei), force=True)
            if gap == 0:
                continue
            if gap.is_positive is not True:
                raise ValueError("energies must be ordered")
            difference = sp.simplify(occupations[i] - occupations[j])
            if difference.is_nonnegative is not True:
                raise ValueError("oscillator map requires passive occupations")
            coupling = sp.Matrix([sp.sqrt(difference) * a[i, j] for a in operators])
            modes.append((gap, sp.ImmutableMatrix(coupling)))
    return tuple(modes)


def oscillator_noise(modes, beta, size):
    """Exact harmonic thermal two-point line weights, complex coupling retained."""
    lines = {}
    for gap, g in modes:
        bose = sp.simplify(1 / (sp.exp(beta * gap) - 1))
        for energy, weight in ((gap, (bose + 1) * g * g.H),
                               (-gap, bose * g.conjugate() * g.T)):
            lines[energy] = lines.get(energy, sp.zeros(size)) + weight
    return {gap: sp.ImmutableMatrix(w) for gap, w in lines.items()}


def oscillator_retarded(modes, z, hbar, size):
    """Both positive and negative poles of the +i commutator convention."""
    result = sp.zeros(size)
    for gap, g in modes:
        result += g * g.H / (gap - hbar * z)
        result += g.conjugate() * g.T / (gap + hbar * z)
    return sp.ImmutableMatrix(result)
