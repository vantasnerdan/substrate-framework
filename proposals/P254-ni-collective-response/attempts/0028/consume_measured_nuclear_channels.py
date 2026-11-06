#!/usr/bin/env python3
"""Publish actual acquired channel tables, evaluated Q and conditional-input gaps.

Default: no flux, no rate, no evolution, no heat. Optional --conditional-input
reads an explicit JSON object with first_encounter/second_encounter (the
ProtonEncounterModel fields and FluxComponent fields), preparation (the four
SequentialPreparation justification fields), initial_population (all seven
species, decimal strings), and duration_s (decimal string). This mode computes
only the declared constant-flux, ground-state, closed-network construction.

Run: python consume_measured_nuclear_channels.py --output PATH
Dependencies for optional evolution: numpy/scipy; default consumption is exact.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
from consume_nuclear_ledger import load_sources, energy_record, label
from substrate_framework.ni_h_measured_channels import (
    FluxComponent, ProtonEncounterModel, channel_domains, load_measured_tables,
    measured_reaction, parenthetical,
)
from substrate_framework.ni_h_nuclear_ledger import daughter_network, mass_difference, reaction_energy
from substrate_framework.ni_h_sequential_channels import (
    CU59, NI56, HE4, SequentialPreparation, build_conditional_network,
    decay_half_lives_seconds, sequence_energy,
)

SOURCE_DIR = ROOT / 'proposals/P254-ni-collective-response/sources/ni-measured-channels-original'
HERE = Path(__file__).resolve().parent


def json_value(value):
    if isinstance(value, Fraction):
        return str(value)
    if is_dataclass(value):
        return json_value(asdict(value))
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    return value


def encounter_from_json(raw: dict) -> ProtonEncounterModel:
    return ProtonEncounterModel(**{**raw, 'components': tuple(FluxComponent(**item) for item in raw['components'])})


def check_ledger_binding(ledger, a, reaction, energy):
    matches = [row for row in ledger['reaction_rows'] if row['channel'] == f'Ni{a}+p->Cu{a+1}+gamma']
    if len(matches) != 1:
        raise ValueError('Existing ledger lacks unique ordinary capture row')
    row = matches[0]
    q = row['q_atomic_ground_state']
    if q['unit'] != 'keV' or Fraction(q['rational']) != energy.q_atomic_kev:
        raise ValueError('Measured-channel / existing ledger mass-Q mismatch')
    if row['entrance'] != [f'Ni{a}', 'H1'] or row['nuclear_products'] != [f'Cu{a+1}'] or not row['radiation']:
        raise ValueError('Existing ledger capture product/radiation mismatch')
    if len(reaction.radiation) == 0:
        raise ValueError('Radiative capture cannot lose its photons')
    return row['channel']


def build_report(conditional_input: dict | None = None) -> dict:
    tables = load_measured_tables(SOURCE_DIR)
    masses, decays = load_sources()
    ledger_bytes = (HERE / 'nuclear-product-energy-ledger.json').read_bytes()
    ledger = json.loads(ledger_bytes)
    registration = json.loads((HERE / 'measured-nuclear-registration.json').read_text())
    if ledger['source_bindings'] != registration['evaluated_inputs']:
        raise ValueError('Existing nuclear ledger uses different evaluated sources')
    channels = []
    for row in tables['capture'] + tables['cu59_palpha']:
        reaction = measured_reaction(row)
        energy = reaction_energy(reaction, masses)
        binding = None
        if row.get('commissioned_even_isotope'):
            binding = check_ledger_binding(ledger, row['target_a'], reaction, energy)
        sigma_raw = row['cross_section_mb_raw']
        sigma, uncertainty = parenthetical(sigma_raw.lstrip('<'))
        channels.append({
            'source_table_row': row,
            'central_cross_section_mb': None if sigma_raw.startswith('<') else sigma,
            'upper_limit_mb': sigma if sigma_raw.startswith('<') else None,
            'parenthetical_cross_section_uncertainty_mb': uncertainty,
            'energy_cm_mev_exact': parenthetical(row['energy_cm_mev_raw'])[0],
            'energy_cm_parenthetical_uncertainty_mev': parenthetical(row['energy_cm_mev_raw'])[1],
            'reactants': [label(key, masses) for key in reaction.reactants],
            'nuclear_products': [label(key, masses) for key in reaction.products],
            'radiation': reaction.radiation, **energy_record(energy, masses),
            'existing_ledger_row': binding,
            'data_scope': 'Inclusive observable, not exclusive gamma/level feeding or complex amplitude',
        })
    seq = sequence_energy(masses, decays)
    alpha_row = next(row for row in ledger['reaction_rows'] if row['channel'] == 'Ni58+2p->Ni+alpha')
    if Fraction(alpha_row['q_atomic_ground_state']['rational']) != seq['net_before_weak_decays'].q_atomic_kev:
        raise ValueError('Sequential net Q disagrees with identical nuclear mass endpoints')
    topology = daughter_network(tuple(key for row in tables['capture'] + tables['cu59_palpha']
                                      for key in measured_reaction(row).products), decays)
    decay_rows = []
    for parent, paths in sorted(topology.items()):
        for branch in paths:
            q = mass_difference((parent,), (branch.daughter,) + branch.nuclear_ejectiles, masses)
            decay_rows.append({'parent': label(parent, masses), 'daughter': label(branch.daughter, masses),
                               'source_record': asdict(decays[parent]), 'path': asdict(branch),
                               **energy_record(q, masses), 'spectra_and_exclusive_feeding': None})
    # Exact Table I literals, NOT adopted energies. The Ni58 printed discrepancy
    # is exposed rather than silently changing Table II or the evaluated ledger.
    reported_q = {58: '3.814', 60: '4.800', 61: '5.866', 62: '6.122', 64: '7.452'}
    q_comparison = []
    for a, source_q in reported_q.items():
        row = next(row for row in tables['capture'] if row['target_a'] == a)
        q = reaction_energy(measured_reaction(row), masses).q_atomic_kev
        difference = Fraction(source_q) * 1000 - q
        q_comparison.append({'target': f'Ni{a}', 'source': 'Simon PDF p3 Table I',
                             'printed_q_mev_raw': source_q, 'ame2020_q_kev': q,
                             'printed_minus_ame2020_kev': difference,
                             'exceeds_half_last_printed_digit': abs(difference) > Fraction(1, 2),
                             'policy': 'AME2020 controls energy accounting; quoted measured Table II unchanged. No gamma spectrum inferred from printed Q.'})
    conditional = None
    if conditional_input is not None:
        network = build_conditional_network(
            tables, decays, encounter_from_json(conditional_input['first_encounter']),
            encounter_from_json(conditional_input['second_encounter']),
            SequentialPreparation(**conditional_input['preparation']))
        conditional = {'explicit_inputs': conditional_input,
                       'capture_rate_per_target_s^-1': network.capture_per_s,
                       'second_rate_per_target_s^-1': network.palpha_per_s,
                       'evolution': network.evolve(conditional_input['initial_population'], conditional_input['duration_s']),
                       'not_a_material_rate': True}
    return json_value({
        'schema': 'P254-A28-measured-NiH-channels-v1',
        'source_acquisition': json.loads((SOURCE_DIR / 'acquisition.json').read_text()),
        'verified_source_bindings': json.loads((SOURCE_DIR / 'bindings.json').read_text()),
        'existing_ledger_sha256': hashlib.sha256(ledger_bytes).hexdigest(),
        'registration': 'measured-nuclear-registration.json',
        'tables': channels, 'source_domains': channel_domains(tables),
        'source_printed_q_comparison': q_comparison,
        'printed_q_original_page_and_ame_inspection': json.loads((SOURCE_DIR / 'capture-q-source-inspection.json').read_text()),
        'original_page_bindings': json.loads((SOURCE_DIR / 'page-image-bindings.json').read_text()),
        'sequential_energy': seq,
        'decay_topology': decay_rows,
        'sequential_half_lives_seconds': {label(key, masses): value for key, value in decay_half_lives_seconds(decays).items()},
        'decay_scope': 'NUBASE inclusive isotope topology; beta+/EC split, nuclear feeding and photon/neutrino spectra remain unspecified',
        'isotope_scope': 'All five Simon stable-Ni tables retained; 61Ni is extra source coverage, not one of the commissioned even-isotope quartet',
        'conditional_network': conditional,
        'actual_material_entrance_distribution': None, 'actual_material_rate': None,
        'correlated_two_proton_amplitude': None, 'ni_receiver_intrinsic_matrix_elements': None,
        'deposited_heat': None, 'escaping_energy': None,
        'missing_inputs': registration['missing_inputs'],
        'scope': 'Concrete measured ordinary binary input and explicitly distinct conditional sequential construction; not quantitative LENR, Kim correlated-2p amplitude or complete N-node closure',
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--conditional-input', type=Path,
                        help='Explicit encounter/preparation/population JSON; omitted means no rate calculation')
    args = parser.parse_args()
    inputs = None if args.conditional_input is None else json.loads(args.conditional_input.read_text())
    text = json.dumps(build_report(inputs), indent=2, allow_nan=False) + '\n'
    if args.output is None:
        print(text, end='')
    else:
        args.output.write_text(text)


if __name__ == '__main__':
    main()
