#!/usr/bin/env python3
"""Independent SymPy/raw-source oracle with explicit exposing mutations.

No precomputed Q outputs are embedded. Literal decimal input anchors below are
transcribed from acquired AME2020 records, not from production-parser outputs.
Run: python proposals/P254-ni-collective-response/attempts/0028/oracle_nuclear_ledger.py
"""
from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
import json
import re

import sympy as sp

from consume_nuclear_ledger import SOURCE_DIR, build_ledger, load_sources
from substrate_framework.ni_h_nuclear_ledger import (
    HE4, H1, NEUTRON, Nuclide, Reaction, daughter_network, decay_paths,
    exact_decimal, mass_difference, nickel_channels, reaction_energy,
    two_body_kinetic_energies,
)

# Original input anchors: AME lines 37,38,44,549,575,590,591,602,617,618,
# 629,643,644,655,669,670. These are input values, not a calculated result.
ANCHORS = {
    (1, 0): '8071.31806', (1, 1): '7288.971064', (4, 2): '2424.91587',
    (56, 28): '-53907.649', (58, 28): '-60228.871',
    (59, 29): '-56358.454', (59, 30): '-47215.678',
    (60, 28): '-64473.243', (61, 29): '-61984.063', (61, 30): '-56348.906',
    (62, 28): '-66746.440', (63, 29): '-65579.868', (63, 30): '-62213.429',
    (64, 28): '-67099.034', (65, 29): '-67263.677', (65, 30): '-65912.024',
}


def rational(value):
    return sp.Rational(value.numerator, value.denominator)


def rejects(function, error=ValueError):
    try:
        function()
    except error:
        return
    raise AssertionError('Exposing mutation was not rejected')


def main():
    masses, decays = load_sources()
    # Separate regex reader over original records, avoiding production slices
    # and exact_decimal. The field widths are the official Fortran specification.
    pattern = re.compile(r'^.(.{3})(.{5})(.{5})(.{5}) (.{3})(.{4}) (.{14})(.{12})')
    independent = {}
    for line in (SOURCE_DIR / 'mass_1.mas20').read_text(encoding='ascii').splitlines()[36:]:
        match = pattern.match(line)
        if not match:
            continue
        _, n, z, a, element, origin, mass, sigma = match.groups()
        assert int(n) + int(z) == int(a)
        if '*' not in mass:
            independent[(int(a), int(z))] = sp.Rational(mass.strip().replace('#', '.'))
    for key, text in ANCHORS.items():
        assert independent[key] == sp.Rational(text)
        assert rational(masses[Nuclide(*key)].excess_kev) == sp.Rational(text)
    assert exact_decimal('28667#') == Fraction(28667)
    assert exact_decimal('*') is None
    # A genuine extrapolated record must keep its flag rather than become data.
    assert masses[Nuclide(3, 3)].extrapolated
    assert masses[Nuclide(3, 3)].excess_raw.strip() == '28667#'

    def m(a, z):
        return independent[(a, z)]

    channels = nickel_channels(include_binary_capture=True)
    assert len(channels) == 16
    mutations = []
    symbolic_balances = []
    u, electron_mass = sp.symbols('u_c2 electron_rest_energy', positive=True)
    for reaction in channels:
        a = reaction.reactants[0].a
        if reaction.mechanism == 'proposed Ni+2p':
            if reaction.name.endswith('Cu+p'):
                expected = m(a, 28) + m(1, 1) - m(a + 1, 29)
                expected_products = (Nuclide(a + 1, 29), H1)
            elif reaction.name.endswith('Ni+alpha'):
                expected = m(a, 28) + 2*m(1, 1) - m(a - 2, 28) - m(4, 2)
                expected_products = (Nuclide(a - 2, 28), HE4)
            elif reaction.name.endswith('Zn+n'):
                expected = m(a, 28) + 2*m(1, 1) - m(a + 1, 30) - m(1, 0)
                expected_products = (Nuclide(a + 1, 30), NEUTRON)
            else:
                raise AssertionError(reaction.name)
        else:
            expected = m(a, 28) + m(1, 1) - m(a + 1, 29)
            expected_products = (Nuclide(a + 1, 29),)
            assert reaction.radiation == ('gamma/cascade',)
        assert reaction.products == expected_products
        actual = reaction_energy(reaction, masses)
        assert rational(actual.q_atomic_kev) == expected
        # Atomic masses include A*u and electron rest masses. Both cancel in
        # these reactions; total electronic binding does NOT cancel in general.
        full = sum(x.a*u + m(x.a, x.z) for x in reaction.reactants) - sum(x.a*u + m(x.a, x.z) for x in reaction.products)
        electron_shift = (sum(x.z for x in reaction.reactants) - sum(x.z for x in reaction.products))*electron_mass
        assert sp.expand(full - expected) == 0 and electron_shift == 0
        assert expected != 0 and -expected != expected  # reversed Q-sign mutation
        if reaction.name.endswith('Zn+n'):
            # Kim source's sign claim is independently tested, not imposed on Q.
            assert (expected > 0) == (a in (62, 64))
        first = reaction.products[0]
        rejects(lambda: reaction_energy(replace(reaction, products=(Nuclide(first.a + 1, first.z),) + reaction.products[1:]), masses))
        rejects(lambda: reaction_energy(replace(reaction, products=(Nuclide(first.a, first.z - 1),) + reaction.products[1:]), masses))
        # Sensitivity/refinement: exact finite perturbation agrees with the
        # isotope's stoichiometric coefficient; rounding cannot silently vanish.
        key, coefficient = actual.coefficients[0]
        delta = Fraction(1, 10**7)
        changed = dict(masses)
        changed[key] = replace(masses[key], excess_kev=masses[key].excess_kev + delta)
        assert reaction_energy(reaction, changed).q_atomic_kev - actual.q_atomic_kev == coefficient*delta
        errors = []
        for places in (0, 3, 6):
            scale = 10**places
            refined = dict(masses)
            for nuclide, _ in actual.coefficients:
                value = masses[nuclide].excess_kev
                refined[nuclide] = replace(masses[nuclide], excess_kev=Fraction(round(value*scale), scale))
            error = abs(reaction_energy(reaction, refined).q_atomic_kev - actual.q_atomic_kev)
            bound = sum(abs(c) for _, c in actual.coefficients)*Fraction(1, 2*scale)
            assert error <= bound
            errors.append({'decimal_places': places, 'q_error_kev': str(error), 'roundoff_bound_kev': str(bound)})
        symbolic_balances.append({'channel': reaction.name, 'Q_kev': str(expected), 'refinement': errors})
        mutations.append({'channel': reaction.name, 'sign_flip_exposed': True, 'wrong_A_rejected': True, 'wrong_Z_rejected': True, 'mass_perturbation_exposed': True})

    # Balanced but wrong exit channel cannot pass a conservation-only oracle.
    for a in (58, 60, 62, 64):
        group = [r for r in channels if r.reactants[0].a == a and r.mechanism == 'proposed Ni+2p']
        cu, alpha, neutron = group
        assert reaction_energy(replace(cu, products=alpha.products), masses).q_atomic_kev != reaction_energy(cu, masses).q_atomic_kev
        assert reaction_energy(replace(neutron, products=cu.products), masses).q_atomic_kev != reaction_energy(neutron, masses).q_atomic_kev
        # Spectator proton cancels in Q and its uncertainty coefficient, but the
        # reaction mechanism is still different from binary gamma capture.
        binary = next(r for r in channels if r.reactants[0].a == a and r.radiation)
        assert reaction_energy(cu, masses) == reaction_energy(binary, masses)
        assert dict(reaction_energy(cu, masses).coefficients)[H1] == 1
        assert cu.mechanism != binary.mechanism and cu.products != binary.products

    # Relativistic two-body oracle independently starts from E1=(W²+m1²-m2²)/2W.
    q, m1, m2 = sp.symbols('q m1 m2', positive=True)
    w = m1 + m2 + q
    first = sp.cancel((w*w + m1*m1 - m2*m2)/(2*w) - m1)
    second = sp.cancel((w*w + m2*m2 - m1*m1)/(2*w) - m2)
    assert sp.simplify(first + second - q) == 0
    assert sp.simplify(first - q*(q + 2*m2)/(2*w)) == 0
    assert sp.simplify(second - q*(q + 2*m1)/(2*w)) == 0
    for energy, rest1, rest2 in ((0, 3, 9), (7, 3, 9), (7, 0, 9), (Fraction(1, 10), 3, 9)):
        left, right = two_body_kinetic_energies(energy, rest1, rest2)
        values = {q: rational(Fraction(energy)), m1: rest1, m2: rest2}
        assert rational(left) == first.subs(values) and rational(right) == second.subs(values)
        assert left + right == energy
        assert two_body_kinetic_energies(energy, rest2, rest1) == (right, left)
    rejects(lambda: two_body_kinetic_energies(-1, 3, 9))
    rejects(lambda: two_body_kinetic_energies(1.0, 3, 9), TypeError)
    # Low-Q limit has inverse mass energy partition; the first correction is
    # retained in the exact kernel rather than assuming mass-number fractions.
    assert sp.simplify(sp.limit(first/q, q, 0) - m2/(m1 + m2)) == 0

    products = tuple(x for r in channels for x in r.products)
    network = daughter_network(products, decays)
    assert network[Nuclide(56, 28)][0].daughter == Nuclide(56, 27)
    assert network[Nuclide(56, 27)][0].daughter == Nuclide(56, 26)
    assert network[Nuclide(59, 29)][0].daughter == Nuclide(59, 28)
    assert network[Nuclide(59, 28)][0].daughter == Nuclide(59, 27)
    assert network[Nuclide(61, 29)][0].daughter == Nuclide(61, 28)
    assert not network[Nuclide(63, 29)] and not network[Nuclide(65, 29)]
    delayed = next(path for path in network[Nuclide(59, 30)] if path.mode == 'B+p')
    assert delayed.daughter == Nuclide(58, 28) and delayed.nuclear_ejectiles == (H1,)
    assert delayed.source_branch_raw == 'B+p=0.10 2' and delayed.inclusive_parent_mode == 'B+'
    assert decays[Nuclide(59, 29)].half_life_raw == '81.5'
    assert decays[Nuclide(61, 29)].half_life_raw == '3.343'
    assert decays[Nuclide(59, 28)].half_life_unit == 'ky'
    assert any('?' in path.source_branch_raw for path in network[Nuclide(58, 28)])
    for parent, paths in network.items():
        for path in paths:
            assert parent.a == path.daughter.a + sum(x.a for x in path.nuclear_ejectiles)
            dz = parent.z - path.daughter.z - sum(x.z for x in path.nuclear_ejectiles)
            assert dz == {'B+': 1, 'B+p': 1, 'B-': -1, '2B+': 2}[path.mode]
            independent_q = m(parent.a, parent.z) - m(path.daughter.a, path.daughter.z) - sum(m(x.a, x.z) for x in path.nuclear_ejectiles)
            actual_q = mass_difference((parent,), (path.daughter,) + path.nuclear_ejectiles, masses)
            assert rational(actual_q.q_atomic_kev) == independent_q
    rejects(lambda: decay_paths(replace(decays[Nuclide(59, 30)], branches_raw='SF=100')))
    ledger = build_ledger()
    assert len(ledger['reaction_rows']) == 16
    assert ledger['calculated_fusion_rate'] is None and ledger['calculated_heat_rate'] is None
    assert all(row['energy_disposal']['deposited_energy'] is None for row in ledger['reaction_rows'])
    print(json.dumps({'status': 'passed', 'source_bindings': ledger['source_bindings'], 'independent_balances': symbolic_balances, 'exposing_mutations': mutations, 'decay_topology_checked': True, 'kinematics_and_refinement_checked': True}, indent=2))


if __name__ == '__main__':
    main()
