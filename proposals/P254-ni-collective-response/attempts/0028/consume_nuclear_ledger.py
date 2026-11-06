#!/usr/bin/env python3
"""Emit the source-backed product/Q ledger; never a rate or heat prediction.

Run from any directory: python path/to/consume_nuclear_ledger.py [--output FILE]
The output contains exact rational keV and decimal renderings. Raw mass/decay
records travel with each isotope; missing transition/transport inputs stay open.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
from substrate_framework.ni_h_nuclear_ledger import (
    Nuclide, daughter_network, mass_difference, nickel_channels,
    reaction_energy, read_ame2020, read_nubase2020,
)

SOURCE_DIR = ROOT / 'proposals/P254-ni-collective-response/sources/ame2020-nubase2020-original'
BINDINGS = {
    'mass_1.mas20': 'e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307',
    'nubase_4.mas20': '1585a5eea86c5e17e90307c7e6e786d060049c4039e392a261ff6db977df9859',
}


def load_sources():
    for name, expected in BINDINGS.items():
        actual = hashlib.sha256((SOURCE_DIR / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Original source hash mismatch: {name}')
    return read_ame2020(SOURCE_DIR / 'mass_1.mas20'), read_nubase2020(SOURCE_DIR / 'nubase_4.mas20')


def label(key: Nuclide, masses) -> str:
    return 'n' if key.z == 0 else f'{masses[key].element}{key.a}'


def energy_json(value: Fraction | None):
    if value is None:
        return None
    with localcontext() as context:
        context.prec = 40
        decimal = str(Decimal(value.numerator) / Decimal(value.denominator))
    return {'rational': str(value), 'decimal': decimal, 'unit': 'keV'}


def energy_record(energy, masses):
    return {
        'q_atomic_ground_state': energy_json(energy.q_atomic_kev),
        'central_value_sign': 'exothermic' if energy.q_atomic_kev > 0 else 'endothermic' if energy.q_atomic_kev < 0 else 'zero',
        'uncertainty_l1_upper_bound': energy_json(energy.uncertainty_l1_kev),
        'uncertainty_meaning': 'Upper bound on standard uncertainty with arbitrary correlations; NOT a confidence interval. AME covariance not acquired.',
        'mass_excess_coefficients': {label(key, masses): coefficient for key, coefficient in energy.coefficients},
        'extrapolated_mass_inputs': [label(key, masses) for key in energy.extrapolated_inputs],
    }


def build_ledger():
    masses, decays = load_sources()
    channels = nickel_channels(include_binary_capture=True)
    network = daughter_network(tuple(key for reaction in channels for key in reaction.products), decays)
    species = set(network)
    species.update(key for reaction in channels for key in reaction.reactants)
    rows = []
    for reaction in channels:
        energy = reaction_energy(reaction, masses)
        row = {
            'channel': reaction.name, 'mechanism': reaction.mechanism,
            'entrance': [label(key, masses) for key in reaction.reactants],
            'nuclear_products': [label(key, masses) for key in reaction.products],
            'radiation': list(reaction.radiation),
            'hydrogen_helium_notation': 'H1 and He4 in mass bookkeeping represent proton and alpha nuclear products, with neutral-atom mass endpoints.',
            'endpoint': 'Ground-state neutral-atom energy reference; no final-state excitation assigned.',
            **energy_record(energy, masses),
            'energy_disposal': {
                'available_cm_energy_expression': 'Q_atomic + entrance_CM_kinetic - final_nuclear_excitation + endpoint_electronic_correction',
                'threshold': 'Positive available CM energy is necessary, not a transmission or reaction probability. No binary-lab threshold is assigned to the three-body entrance.',
                'particle_partition': 'Two-body CM kinematics available via two_body_kinetic_energies with supplied consistent product rest energies; unspecified excitation/cascades and entrance momentum prevent assigning actual lab spectra.',
                'transport_needed': ['charged-product stopping and recoil', 'neutron scattering/capture/escape where present', 'gamma conversion/absorption/escape', 'beta/positron stopping and annihilation', 'neutrino energy escape'],
                'deposited_energy': None,
            },
        }
        rows.append(row)
    decay_rows = []
    for parent, paths in sorted(network.items()):
        for path in paths:
            energy = mass_difference((parent,), (path.daughter,) + path.nuclear_ejectiles, masses)
            decay_rows.append({
                'parent': label(parent, masses), 'daughter': label(path.daughter, masses),
                'mode': path.mode, 'source_branch_raw': path.source_branch_raw,
                'inclusive_parent_mode': path.inclusive_parent_mode,
                'qualification': 'Unobserved/unknown qualifier retained; not assigned a probability' if '?' in path.source_branch_raw else 'Evaluated aggregate mode; feeding/spectra unresolved',
                'nuclear_ejectiles': [label(key, masses) for key in path.nuclear_ejectiles],
                'leptons': list(path.leptons), 'energy_convention': path.energy_convention,
                'atomic_mass_difference': energy_json(energy.q_atomic_kev),
                'uncertainty_l1_upper_bound': energy_json(energy.uncertainty_l1_kev),
                'extrapolated_mass_inputs': [label(key, masses) for key in energy.extrapolated_inputs],
                'excited_daughter_cascade': 'NUBASE aggregate branches do not resolve gamma/internal-conversion feeding; ENSDF level and radiation data needed.',
                'energy_is_not_heat': True,
            })
    inputs = {}
    for key in sorted(species):
        mass, decay = masses[key], decays[key]
        inputs[label(key, masses)] = {
            'A': key.a, 'Z': key.z,
            'mass_excess': energy_json(mass.excess_kev),
            'mass_uncertainty': energy_json(mass.uncertainty_kev),
            'mass_extrapolated': mass.extrapolated,
            'ame_source_line': mass.line_number, 'ame_raw_line': mass.raw_line,
            'nubase': {k: v for k, v in asdict(decay).items() if k != 'nuclide'},
            'isotope_abundance_scope': 'IS in raw NUBASE data is evaluated natural abundance, NOT specimen composition.',
        }
    return {
        'schema': 'P254-NiH-source-backed-product-energy-ledger-v1',
        'source_bindings': BINDINGS,
        'source_acquisition': 'nuclear-source-acquisition.json',
        'documentation': 'nuclear-source-registration.json',
        'units': 'keV; exact decimal/rational input, no binary floating arithmetic',
        'reaction_rows': rows, 'decay_topology': decay_rows, 'isotope_inputs': inputs,
        'decay_topology_semantics': 'Inclusive B+ means EC plus positron decay; B+p is its subset, never an extra independent probability. Unobserved ? modes remain qualified. No branching network time/yield calculation.',
        'atomic_reference_correction': 'For bare nuclei Q_nuclear = Q_atomic + sum(B_e_initial) - sum(B_e_final), where B_e is positive total electronic binding. Material ionization/bond/charge changes require corresponding endpoint energies. Electron rest masses cancel for A,Z-conserving reactions.',
        'missing_physics': [
            'Specimen isotope fractions and nuclear entrance energy/spin/position state.',
            'A supported two-proton spatial/spin correlation, lifetime, population, and three-body transition amplitude OR a justified sequential channel model. Opposite spins alone do not bind protons.',
            'Channel-resolved cross sections/transmissions, final excited-state populations, angular distributions and branching; ordinary binary capture is not the proposed Ni+2p amplitude.',
            'Ni receiver level energies/J^pi/widths plus spin-momentum transition density for the specified nuclear operator; a gamma strength is not automatically the Hagelstein a-operator matrix.',
            'Nuclear-to-medium coupling and material spectral support at the nuclear transition frequency/momentum, not just valence-response coherence.',
            'ENSDF daughter level feeding/radiation and EC/positron split; electron binding corrections, stopping/capture, geometry and escape for deposited energy.',
        ],
        'calculated_fusion_rate': None, 'calculated_heat_rate': None,
        'scope': 'Energetics and product topology only. Positive Q, nuclear absorption, and norm loss do not establish a rate or deposited heat.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    text = json.dumps(build_ledger(), indent=2) + '\n'
    if args.output:
        args.output.write_text(text, encoding='utf-8')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
