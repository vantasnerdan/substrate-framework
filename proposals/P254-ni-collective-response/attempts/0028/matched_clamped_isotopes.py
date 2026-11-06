"""Saved NIST atomic-isotope masses as explicit inertial-only H1/H2 variants.

The original case is never silently relabeled protium. These atomic masses
are an effective ionic mass substitution, not a nuclear Hamiltonian change.
"""
from pathlib import Path
import re

import numpy as np
from material_run_contract import binding

SOURCE = Path(__file__).resolve().parents[2]/'sources/nist-hydrogen-isotopes.txt'


def nist_hydrogen_masses(path=SOURCE):
    text = Path(path).read_text()
    values = {}
    for block in text.split('Atomic Number = ')[1:]:
        number = re.search(r'^1\s', block)
        isotope = re.search(r'Mass Number = ([12])\s', block)
        mass = re.search(r'Relative Atomic Mass = ([0-9.]+)\(([0-9]+)\)', block)
        if number and isotope and mass:
            key = 'H'+isotope.group(1)
            if key in values:
                raise ValueError('duplicate saved NIST isotope entry')
            values[key] = float(mass.group(1))
    if set(values) != {'H1', 'H2'} or not 0 < values['H1'] < values['H2']:
        raise ValueError('saved original NIST H1/H2 mass entries required')
    return values


def isotope_masses(arrays, variant):
    if variant not in ('original', 'H1', 'H2'):
        raise ValueError('explicit original, H1 or H2 mass variant required')
    original = np.asarray(arrays['atomic_masses_amu'], dtype=float)
    names = np.asarray(arrays['atom_names']).astype(str)
    hydrogen = names == 'H'
    masses = original.copy()
    values = nist_hydrogen_masses()
    if variant != 'original':
        if not np.any(hydrogen):
            raise ValueError('H isotope variant requires original H atoms')
        masses[hydrogen] = values[variant]
    return masses, dict(variant=variant, original_masses_amu=original,
        effective_masses_amu=masses, hydrogen_atom_indices=np.asarray(arrays['atom_indices'])[hydrogen],
        nist_source=binding(SOURCE), nist_relative_atomic_masses_amu=values,
        original_H_exactly_NIST_H1=bool(np.any(hydrogen) and np.all(original[hydrogen] == values['H1'])),
        scope='Mass-only effective ionic substitution; fixed physical potential, geometry, charge and electronic preparation. Not nuclear physics or convergence.')
