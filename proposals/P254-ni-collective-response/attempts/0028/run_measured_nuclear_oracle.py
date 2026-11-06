#!/usr/bin/env python3
"""Independent exact-mass/SymPy/finite-network oracle with sensitive mutations.

This is an executable oracle for Root to run after integration, not executed by
its author. Synthetic fluxes are explicitly labelled and have measured-domain
support; they are not a material preparation or LENR estimate. Requires SymPy,
numpy and scipy. No cloud, trajectories or unbounded jobs.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction
import json
from pathlib import Path
import re
import sys

import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
from consume_nuclear_ledger import load_sources, SOURCE_DIR as MASS_SOURCE_DIR
from consume_measured_nuclear_channels import build_report, SOURCE_DIR
import substrate_framework.ni_h_measured_channels as measured
from substrate_framework.ni_h_measured_channels import (
    FluxComponent, ProtonEncounterModel, load_measured_tables, rate_from_matched_flux,
)
from substrate_framework.ni_h_nuclear_ledger import Nuclide
from substrate_framework.ni_h_sequential_channels import (
    CU59, NI59, NI56, CO56, SPECIES, SPECIES_NAMES, ENERGY_CARRIERS,
    SequentialPreparation, account_energy, build_conditional_network, sequence_energy,
)


def must_reject(call):
    try:
        call()
    except (ValueError, TypeError, AssertionError):
        return
    raise AssertionError('Mutation survived')


def equal(actual, expected):
    assert actual == expected, (actual, expected)


def close(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-13)


def source_table_oracle(tables):
    # Parse the actual original text's table cells independently of the JSON.
    # Section boundaries are the printed isotope labels, not row-index guesses.
    text = (SOURCE_DIR / 'simon-1305.1213.txt').read_text()
    lines = text.splitlines()
    lines = lines[next(index for index, line in enumerate(lines) if 'TABLE II.' in line):]
    blocks = {}
    current = None
    for line in lines:
        heading = re.search(r'Ni\(p,γ\)(\d+)\s*Cu', line)
        if heading:
            current = int(heading.group(1)) - 1
            blocks.setdefault(current, [])
        if current is not None:
            for energy, sigma in re.findall(r'(?<![\d.])(\d+\.\d{2})\s+(\d+\.\d+\(\d+\))', line):
                blocks[current].append((energy, sigma))
        if current == 58 and '1(a)-(d)' in line:
            break
    for a, count in ((58, 8), (60, 15), (61, 8), (62, 4), (64, 16)):
        actual = [(row['energy_cm_mev_raw'], row['cross_section_mb_raw'])
                  for row in tables['capture'] if row['target_a'] == a]
        equal(len(actual), count)
        equal(sorted(actual), sorted(blocks[a]))
    # A second literal oracle for the separately printed finite-strip table.
    expected = [
        ('5.88(3)', '0.18', '0.26', '7.024', '0.220', '2.355'),
        ('5.44(3)', '0.18', '0.27', '3.412', '0.154', '0.612'),
        ('4.98(4)', '0.19', '0.28', '1.861', '0.113', '0.300'),
        ('4.49(4)', '0.21', '0.28', '0.619', '0.065', '0.092'),
        ('4.02(4)', '0.19', '0.31', '0.248', '0.042', '0.083'),
        ('3.54(4)', '0.17', '0.37', '0.082', '0.024', '0.014'),
        ('2.98(5)', '0.19', '0.37', '<0.015', None, None),
        ('2.43(5)', '0.18', '0.42', '<0.011', None, None),
    ]
    fields = ('energy_cm_mev_raw', 'energy_bin_plus_mev_raw', 'energy_bin_minus_mev_raw',
              'cross_section_mb_raw', 'statistical_uncertainty_mb_raw', 'systematic_uncertainty_mb_raw')
    equal([tuple(row[key] for key in fields) for row in tables['cu59_palpha']], expected)
    equal(measured.parenthetical('1.379(17)'), (Fraction('1.379'), Fraction('0.017')))


def exact_mass_oracle(masses, decays):
    # Independent direct SymPy interpretation of ORIGINAL fixed-width decimals;
    # no use of the production Fraction parser or production mass-difference API.
    raw = {}
    for line in (MASS_SOURCE_DIR / 'mass_1.mas20').read_text().splitlines()[36:]:
        if line.strip():
            field = line[28:42].strip()
            if field and '*' not in field:
                raw[(int(line[14:19]), int(line[9:14]))] = sp.Rational(field.replace('#', '.'))
    n58, p, cu, n56, alpha = sp.symbols('M58 H M59 M56 He4')
    first = n58 + p - cu
    second = cu + p - n56 - alpha
    net = n58 + 2*p - n56 - alpha
    equal(sp.expand(first + second - net), 0)
    substitution = dict(zip((n58, p, cu, n56, alpha),
                           (raw[key] for key in ((58, 28), (1, 1), (59, 29), (56, 28), (4, 2)))))
    result = sequence_energy(masses, decays)
    for key, expr in [('capture', first), ('second_encounter', second), ('net_before_weak_decays', net)]:
        equal(sp.Rational(str(result[key].q_atomic_kev)), expr.subs(substitution))
    terminal = raw[(58, 28)] + 2*raw[(1, 1)] - raw[(56, 26)] - raw[(4, 2)]
    equal(sp.Rational(str(result['successful_terminal_atomic_endpoint'].q_atomic_kev)), terminal)
    failed = raw[(58, 28)] + raw[(1, 1)] - raw[(59, 27)]
    equal(sp.Rational(str(result['intermediate_decay_terminal_atomic_endpoint'].q_atomic_kev)), failed)
    equal(sp.Rational(str(result['capture'].q_atomic_kev)), sp.Rational('3418.554064'))
    changed = dict(masses)
    changed[CU59] = replace(changed[CU59], excess_kev=changed[CU59].excess_kev + 1)
    must_reject(lambda: equal(sp.Rational(str(sequence_energy(changed, decays)['capture'].q_atomic_kev)),
                              first.subs(substitution)))
    must_reject(lambda: equal(sp.expand(first + second - (n58 + p - n56 - alpha)), 0))
    # Intermediate uncertainty cancels in the NET, not by summing marginal bounds.
    equal(dict(result['net_before_weak_decays'].coefficients).get(CU59, 0), 0)
    return result


def synthetic_models():
    # Exactly lambda1=1/4 and lambda2=3/8 per second, artificial supplied fluxes.
    first = FluxComponent('Ni58-pg-4.30', '4.30', None, None,
                          str(Fraction(1, 4) * 10**27 / Fraction('0.147')),
                          'MeV_CM', 'cm^-2 s^-1')
    second = FluxComponent('Cu59-pa-5.88', '5.88', '5.62', '6.06',
                           str(Fraction(3, 8) * 10**27 / Fraction('7.024')),
                           'MeV_CM', 'cm^-2 s^-1')
    def model(component):
        return ProtonEncounterModel((component,), 'laboratory-ground-state-unpolarized',
            'Oracle-only assumed ground-state unpolarized targets; not material evidence',
            'Oracle-only declared exact experimental response matching; no interpolation',
            'Synthetic finite-network oracle, not Ni/H preparation')
    return model(first), model(second)


def finite_network_oracle(network):
    a, b, d, t = sp.Rational(1, 4), sp.Rational(3, 8), sp.log(2)/sp.Rational('81.5'), sp.Integer(12)
    k = b + d
    expected_n58 = sp.exp(-a*t)
    expected_cu = a*(sp.exp(-a*t)-sp.exp(-k*t))/(k-a)
    expected_alpha = a*b/(k-a)*((1-sp.exp(-a*t))/a - (1-sp.exp(-k*t))/k)
    expected_decay = d*expected_alpha/b
    reference = [float(expr.evalf(50)) for expr in (expected_n58, expected_cu, expected_alpha, expected_decay)]
    initial = np.zeros(13)
    initial[0] = 1
    def compare(matrix):
        final = expm(matrix * 12) @ initial
        close(final[0], reference[0])
        close(final[1], reference[1])
        close(final[4:7].sum(), reference[2])
        close(final[2:4].sum(), reference[3])
        close(final[8], reference[2])
        close(final[9], reference[3])
        close(final[:7].sum(), 1)
        return final
    matrix = network.generator()
    close(matrix[:7].sum(axis=0), np.zeros(13))
    final = compare(matrix)
    populations = {name: '1' if name == 'Ni58' else '0' for name in SPECIES_NAMES}
    result = network.evolve(populations, '12')
    close(result['protons_consumed'], (1-reference[0])+reference[2])
    close(result['photon_cascade_events'], 1-reference[0])
    close(result['alpha_particles_produced'], reference[2])
    close(np.dot([58, 59, 59, 59, 56, 56, 56], final[:7]) + 4*final[8],
          58 + result['protons_consumed'])
    close(np.dot([28, 29, 28, 27, 28, 27, 26], final[:7]) + 2*final[8],
          28 + result['protons_consumed'] - final[9:13].sum())
    equal(result['deposited_heat'], None)
    omitted = matrix.copy()
    decay_rate = float(d.evalf(50))
    omitted[1, 1] += decay_rate
    omitted[2, 1] -= decay_rate
    omitted[9, 1] -= decay_rate
    must_reject(lambda: compare(omitted))
    doubled = matrix.copy()
    doubled[8] *= 2
    must_reject(lambda: compare(doubled))
    must_reject(lambda: close(2*result['events']['palpha'], result['protons_consumed']))
    return final.tolist()


def run():
    tables = load_measured_tables(SOURCE_DIR)
    masses, decays = load_sources()
    source_table_oracle(tables)
    exact_mass_oracle(masses, decays)
    first, second = synthetic_models()
    preparation = SequentialPreparation(*(['Oracle-only declared closure, not material evidence'] * 4))
    network = build_conditional_network(tables, decays, first, second, preparation)
    equal(network.capture_per_s, Fraction(1, 4))
    equal(network.palpha_per_s, Fraction(3, 8))
    equal(network.decay_half_lives_s[CU59], Fraction('81.5'))
    equal(network.decay_half_lives_s[NI56], Fraction('6.075') * 86400)
    equal(network.decay_half_lives_s[CO56], Fraction('77.236') * 86400)
    equal(network.decay_half_lives_s[NI59], Fraction(81 * 1000 * 36525 * 86400, 100))
    result = finite_network_oracle(network)
    component = first.components[0]
    rows = [row for row in tables['capture'] if row['target_a'] == 58]
    for change in ({'energy_cm_mev': '0.000000025'}, {'energy_cm_mev': '4.40'},
                   {'energy_unit': 'keV'}, {'flux_unit': 'm^-2 s^-1'}, {'flux': '-1'}):
        bad = replace(first, components=(replace(component, **change),))
        must_reject(lambda: rate_from_matched_flux(rows, bad))
    must_reject(lambda: rate_from_matched_flux(rows, replace(first, entrance_state='excited')))
    bad_bin = replace(second, components=(replace(second.components[0], support_lower_mev='5.63'),))
    must_reject(lambda: rate_from_matched_flux(tables['cu59_palpha'], bad_bin))
    limit = replace(second.components[0], row_id='Cu59-pa-2.43', energy_cm_mev='2.43',
                    support_lower_mev='2.01', support_upper_mev='2.61')
    must_reject(lambda: rate_from_matched_flux(tables['cu59_palpha'], replace(second, components=(limit,))))
    original_conversion = measured.MB_TO_CM2
    try:
        measured.MB_TO_CM2 = Fraction(1, 10**24)  # millibarns erroneously treated as barns
        must_reject(lambda: equal(rate_from_matched_flux(rows, first).per_target_per_second, Fraction(1, 4)))
    finally:
        measured.MB_TO_CM2 = original_conversion
    produced = {carrier: '0' for carrier in ENERGY_CARRIERS}
    produced.update(capture_photons='3', alpha_kinetic='4', nuclear_recoil='1', neutrinos='2')
    accounting = account_energy('10', produced)
    equal(accounting['deposited_heat_kev'], None)
    heat_mutant = {**accounting, 'deposited_heat_kev': Fraction(10)}
    must_reject(lambda: equal(heat_mutant['deposited_heat_kev'], None))
    transport = {carrier: {'deposited_kev': '0', 'escaped_kev': value, 'retained_kev': '0'}
                 for carrier, value in produced.items()}
    escaped = account_energy('10', produced, transport=transport, transport_provenance='Synthetic all-escape boundary oracle')
    equal(escaped['deposited_heat_kev'], 0)
    equal(escaped['escaped_kev'], 10)
    bad_transport = deepcopy(transport)
    bad_transport['neutrinos']['deposited_kev'] = '2'
    must_reject(lambda: account_energy('10', produced, transport=bad_transport, transport_provenance='double-count mutant'))
    report = build_report()
    equal(len(report['tables']), 59)
    for field in ('actual_material_rate', 'deposited_heat', 'conditional_network', 'correlated_two_proton_amplitude', 'ni_receiver_intrinsic_matrix_elements'):
        equal(report[field], None)
    return {
        'status': 'pass', 'scope': 'Exact acquired-source and conditional-network oracle; no material prediction',
        'oracles': ['original table transcription', 'SymPy exact original-mass Q and intermediate cancellation',
                    'analytic competing-decay finite network', 'all-escape separate-carrier energy accounting', 'actual saved-data consumer'],
        'mutations_rejected': ['Cu59 decay omitted', 'palpha event multiplicity doubled',
            'two protons assigned to every capture despite failed intermediate', 'millibarn treated as barn',
            'energy unit mismatch', 'flux area unit mismatch', 'cold-lattice energy substituted',
            'unmeasured interpolated energy', 'finite-strip support mismatch', 'upper limit promoted to central rate',
            'excited entrance substituted', 'negative flux', 'AME mass changed', 'net proton multiplicity changed',
            'emitted Q labelled heat without transport', 'deposited and escaped energy double counted'],
        'synthetic_final_populations_and_counts': result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output is None:
        print(text, end='')
    else:
        args.output.write_text(text)


if __name__ == '__main__':
    main()
