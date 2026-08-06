from __future__ import annotations

import sympy as sp

import substrate_framework as framework
from substrate_framework.sg_spectral_induction import (
    KINK_NORMALIZED_MASS,
    PHONON_NORMALIZED_MASS,
    TOWER_EDGE_ACTION,
    induced_inverse_g_species,
    maximum_lattice_level,
    newton_q_normalized,
    sg_mass_tower,
    tower_edge_cutoff,
    tower_inverse_g_normalized,
)


def test_sg_spectral_induction_api_is_exported_from_package() -> None:
    assert framework.tower_inverse_g_normalized is tower_inverse_g_normalized
    assert framework.newton_q_normalized is newton_q_normalized
    assert framework.sg_mass_tower is sg_mass_tower
    assert framework.induced_inverse_g_species is induced_inverse_g_species


def test_per_species_coefficient_matches_imported_a2_matching() -> None:
    m = sp.Symbol("m", positive=True)
    lam = sp.Symbol("Lambda", positive=True)
    exact = induced_inverse_g_species(m, 0, lam)
    assert sp.simplify(exact - m**2 * sp.log(lam**2 / m**2) / (12 * sp.pi)) == 0
    xi = sp.Symbol("xi", real=True)
    general = induced_inverse_g_species(m, xi, lam)
    assert sp.simplify(
        general - (sp.Rational(1, 6) - xi) * m**2 * sp.log(lam**2 / m**2) / (2 * sp.pi)
    ) == 0


def test_conformal_coupling_mutation_annihilates_contribution() -> None:
    # xi = 1/6 must kill the induced R term exactly: the verifier is sensitive
    # to the sign and value of the nonminimal coupling.
    assert induced_inverse_g_species(3, sp.Rational(1, 6), 10) == 0


def test_tower_edge_cutoff_is_the_accepted_action_edge() -> None:
    assert tower_edge_cutoff() == TOWER_EDGE_ACTION == 8 * sp.pi


def test_every_species_logarithm_is_strictly_positive() -> None:
    cutoff = tower_edge_cutoff()
    for species in sg_mass_tower(1):
        ratio = cutoff / species.normalized_mass
        assert float(ratio) > 1.0, species.name


def test_tower_composition_matches_declared_premises() -> None:
    tower = sg_mass_tower(1)
    assert tower[0].normalized_mass == KINK_NORMALIZED_MASS
    assert tower[0].multiplicity == 2
    breather_levels = [s for s in tower if s.name.startswith("breather_n")]
    # C-SG-007: n * h < 8*pi with h = 1 gives n = 1..25.
    assert len(breather_levels) == 25
    assert breather_levels[0].normalized_mass == 16 * sp.sin(sp.Rational(1, 16))


def test_phonon_is_not_double_counted() -> None:
    # The phonon is the n=1 breather level (DHN identification); a standalone
    # phonon species would double-count the lightest excitation.
    tower = sg_mass_tower(1)
    assert not any(s.name == "phonon" for s in tower)
    n1 = [s for s in tower if s.name == "breather_n1"][0]
    # The n=1 level sits at the phonon gap scale (16*sin(1/16) ~ 1 at h=1).
    assert abs(float(n1.normalized_mass) - float(PHONON_NORMALIZED_MASS)) < 1e-3


def test_maximum_lattice_level_obeys_the_open_domain() -> None:
    assert maximum_lattice_level(1) == 25
    assert maximum_lattice_level(2) == 12
    h = sp.Rational(8) * sp.pi / 10  # exactly ten levels would touch the edge
    assert maximum_lattice_level(h) == 9


def test_kink_multiplicity_mutation_moves_the_sum_by_the_kink_term() -> None:
    with_kinks = tower_inverse_g_normalized(1, 0, 2)
    without_kinks = tower_inverse_g_normalized(1, 0, 0)
    kink_term = 2 * induced_inverse_g_species(KINK_NORMALIZED_MASS, 0, tower_edge_cutoff())
    assert sp.simplify(with_kinks - without_kinks - kink_term) == 0
    assert float(with_kinks) > float(without_kinks)


def test_action_quantum_mutation_removes_levels() -> None:
    s_h1 = tower_inverse_g_normalized(1)
    s_h2 = tower_inverse_g_normalized(2)
    assert len(sg_mass_tower(2)) == 1 + 12
    assert float(s_h2) < float(s_h1)


def test_q_is_the_exact_reciprocal_of_the_tower_sum() -> None:
    s_total = tower_inverse_g_normalized(1)
    q = newton_q_normalized(1)
    assert sp.simplify(q * s_total - 1) == 0


def test_symbolic_sum_matches_direct_float_summation() -> None:
    # Independent numeric route: straight float loop over the same premises.
    import math

    cutoff = 8 * math.pi
    total = 2 * (8.0**2) * math.log(cutoff**2 / 8.0**2) / (12 * math.pi)
    for n in range(1, 26):
        e_n = 16 * math.sin(n / 16)
        total += e_n**2 * math.log(cutoff**2 / e_n**2) / (12 * math.pi)
    assert abs(float(tower_inverse_g_normalized(1)) - total) < 1e-10 * total


def test_h_tending_downward_grows_the_tower_without_bound() -> None:
    # The h -> 0 continuum limit must be diagnosed, not hidden: halving h
    # roughly doubles the level count and the sum grows accordingly.
    s_half = tower_inverse_g_normalized(sp.Rational(1, 2))
    assert len(sg_mass_tower(sp.Rational(1, 2))) == 1 + 50
    assert float(s_half) > 1.5 * float(tower_inverse_g_normalized(1))
