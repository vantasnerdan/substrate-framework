"""Exact sine-Gordon mass-tower induction of the Newton coupling.

This module closes the dimensionless coordinate that
:mod:`substrate_framework.induced_gravity` and claim ``C-GRV-001`` proved
unreachable by dimensional analysis alone.  Dimensional analysis forces only

``G = q * a**2 * c**3 / hbar``

with ``q`` dimensionless and free.  Here ``q`` is *derived* from the accepted
sine-Gordon spectrum: the one-loop vacuum polarization of each mass species
on the emergent geometry induces an Einstein-Hilbert term through the
Seeley-DeWitt ``a2`` heat-kernel coefficient, and the tower edge ``E = 8*pi``
of the accepted action lattice (``C-SG-007``) acts as the spectrum's own
ultraviolet boundary, so every level's logarithm ``log(Lambda**2/m**2)`` is
positive and finite with ``Lambda`` the tower edge.

The imported standard result (role explicit): a real scalar of mass ``m``
and nonminimal coupling ``xi`` contributes, at one loop,

``Gamma_ind ⊃ integral sqrt(-g) * R * (1/6 - xi) * m**2 / (2 * (4*pi)**2) * log(Lambda**2/m**2)``

in inverse-Compton-length units for ``m``.  Matching to the Einstein-Hilbert
normalization ``integral sqrt(-g) * R / (16*pi*G)`` gives the per-species
inverse coupling below.  The regulator is the tower-edge cutoff; the
logarithm's additive scheme constant is absorbed by defining ``Lambda`` as
the tower edge exactly (no floating scheme parameter).

This module does not derive the action lattice premise, the quantum
identification ``J_scale = hbar``, the Gordon-metric emergence of the
geometry, the breather-as-loop-species semiclassical premise, or any
empirical confrontation.  Those are declared premises of the composing
campaign, not outputs of this ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp

from .sine_gordon import breather_action_lattice_energy

TOWER_EDGE_ACTION = 8 * sp.pi
"""Accepted upper edge of the C-SG-007 breather action domain ``J < 8*pi``."""

PHONON_NORMALIZED_MASS = sp.Integer(1)
"""Normalized sG phonon gap of ``U_tt - U_xx + sin U``; the phonon is the
``n -> 0``-style light level of the tower and enters the loop sum as the
``n = 1`` breather level, never as a separate species."""

KINK_NORMALIZED_MASS = sp.Integer(8)
"""Normalized static-kink energy of the C-SG-002 Hamiltonian conventions."""


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) > 0.0:
            raise ValueError(f"{name} must be real and positive")
    return expression


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) >= 0.0:
            raise ValueError(f"{name} must be real and nonnegative")
    return expression


def induced_inverse_g_species(mass: Any, xi: Any, cutoff: Any) -> sp.Expr:
    """One species's contribution to ``1/G`` in inverse-length-squared units.

    With ``m``, ``Lambda`` inverse Compton lengths and ``xi`` the nonminimal
    coupling, the imported a2 coefficient matched to the Einstein-Hilbert
    normalization gives

    ``(1/G)_s = (1/6 - xi) * m**2 * log(Lambda**2 / m**2) / (2*pi)``.

    At the canonical-scalar premise ``xi = 0`` this is
    ``m**2 * log(Lambda**2 / m**2) / (12*pi)``.
    """

    m = _positive_exact(mass, "mass")
    lam = _positive_exact(cutoff, "cutoff")
    x = sp.sympify(xi)
    return (sp.Rational(1, 6) - x) * m**2 * sp.log(lam**2 / m**2) / (2 * sp.pi)


def tower_edge_cutoff() -> sp.Expr:
    """The spectrum's own ultraviolet boundary in normalized mass units.

    The C-SG-007 lattice terminates at action ``J = 8*pi`` where the breather
    energy ``16*sin(J/16)`` reaches ``16``; the heaviest admissible state of
    the accepted spectrum fixes the cutoff.  Expressed as a normalized mass
    the boundary is the edge action value ``8*pi`` itself: no state sits at
    or above it, so every species logarithm is strictly positive.
    """

    return TOWER_EDGE_ACTION


def maximum_lattice_level(action_quantum: Any) -> sp.Expr:
    """Largest admissible integer ``n`` with ``n*h < 8*pi`` (C-SG-007)."""

    h = _positive_exact(action_quantum, "action_quantum")
    return sp.floor(TOWER_EDGE_ACTION / h - sp.Integer(1) / 10**12)


@dataclass(frozen=True)
class TowerSpecies:
    """One loop species: normalized mass and declared multiplicity."""

    name: str
    normalized_mass: sp.Expr
    multiplicity: int


def sg_mass_tower(action_quantum: Any, kink_multiplicity: int = 2) -> tuple[TowerSpecies, ...]:
    """The accepted sine-Gordon spectrum as a loop-species tower.

    Species: the kink-antikink pair (normalized mass 8, declared
    multiplicity, default 2) and every admissible C-SG-007 breather level
    ``n`` with normalized mass ``E_n = 16*sin(n*h/16)``.  The phonon is not a
    separate species: it is the ``n = 1`` breather level (the accepted
    Dashen-Hasslacher-Neveu identification of the elementary excitation with
    the lightest breather), so listing it separately would double-count.
    The breather-as-loop-species premise and the kink multiplicity are
    declared campaign premises; both are exposed for mutation testing.
    """

    if not isinstance(kink_multiplicity, int) or kink_multiplicity < 0:
        raise ValueError("kink_multiplicity must be a nonnegative integer")
    h = _positive_exact(action_quantum, "action_quantum")
    species: list[TowerSpecies] = [
        TowerSpecies("kink_antikink", KINK_NORMALIZED_MASS, kink_multiplicity),
    ]
    n_max = maximum_lattice_level(h)
    if not n_max.is_number:
        raise ValueError("action_quantum must be numeric for tower enumeration")
    for n in range(1, int(n_max) + 1):
        energy = breather_action_lattice_energy(n, h)
        species.append(TowerSpecies(f"breather_n{n}", sp.sympify(energy), 1))
    return tuple(species)


def tower_inverse_g_normalized(
    action_quantum: Any, xi: Any = 0, kink_multiplicity: int = 2
) -> sp.Expr:
    """Exact normalized induced ``1/G`` from the full tower.

    The sum

    ``sum_s mult_s * (1/6 - xi) * m_s**2 * log((8*pi)**2 / m_s**2) / (2*pi)``

    in units of the squared phonon inverse Compton length ``(omega_0/c)**2``.
    """

    cutoff = tower_edge_cutoff()
    total = sp.Integer(0)
    for species in sg_mass_tower(action_quantum, kink_multiplicity):
        total += species.multiplicity * induced_inverse_g_species(
            species.normalized_mass, xi, cutoff
        )
    return total


def newton_q_normalized(action_quantum: Any, xi: Any = 0, kink_multiplicity: int = 2) -> sp.Expr:
    """The C-GRV-001 dimensionless coordinate ``q`` as a pure tower number.

    With the physical scale map ``m_s = E_s * omega_0 / c`` (inverse Compton
    lengths; ``E_s`` normalized tower masses) the induced relation reads

    ``1/G = (c**3/hbar) * (omega_0/c)**2 * S,``

    where ``S`` is :func:`tower_inverse_g_normalized`.  Matching to the
    C-GRV-001 monomial ``G = q * ell**2 * c**3 / hbar`` with the C-MED-003
    length ``ell = c/omega_0`` gives exactly ``q = 1/S``: a pure number of
    the upstream spectrum, the action quantum, and the declared
    multiplicities.
    """

    return 1 / tower_inverse_g_normalized(action_quantum, xi, kink_multiplicity)


def species_contribution_table(
    action_quantum: Any, xi: Any = 0, kink_multiplicity: int = 2
) -> tuple[dict[str, Any], ...]:
    """Per-species ledger rows for audit: mass, multiplicity, contribution."""

    cutoff = tower_edge_cutoff()
    rows = []
    for species in sg_mass_tower(action_quantum, kink_multiplicity):
        contribution = species.multiplicity * induced_inverse_g_species(
            species.normalized_mass, xi, cutoff
        )
        rows.append(
            {
                "name": species.name,
                "normalized_mass": species.normalized_mass,
                "multiplicity": species.multiplicity,
                "contribution": contribution,
            }
        )
    return tuple(rows)
