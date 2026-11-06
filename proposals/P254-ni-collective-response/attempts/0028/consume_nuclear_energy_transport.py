#!/usr/bin/env python3
"""Consume acquired primary NIST tables, optionally an explicit energy preparation.

No input means an actual original-table inventory, NOT a material heat run.
Input JSON is deliberately explicit. Photon nodes represent supplied discrete
lines or supplied quadrature nodes; count already includes any quadrature weight.
No interpolated spectrum or event rate is manufactured. Example input:

{
  "normalization": "one hypothetical emitted photon; not a measured reaction",
  "spectrum_evidence": "explicit conditional 1 MeV line, not an observed cascade",
  "expected_emitted_energy_mev": 1.0,
  "carriers": [{
    "kind": "photon", "energy_mev": 1.0, "count": 1.0,
    "paths": [{"probability": 1.0,
      "geometry_evidence": "hypothetical ray to external boundary, 1 g/cm2 Ni",
      "layers": [{"areal_density_g_cm2": 1.0, "mass_fractions": {"Ni": 1.0}}]
    }]
  }]
}

Photon: optional edge_side="below"/"above" (required at an elemental edge).
Every path requires probability, geometry_evidence and ordered layers; [] means
explicit vacuum. Each layer takes mass_fractions {"Ni":wNi,"H":wH}, sum=1,
NOT atomic fractions, and areal_density_g_cm2 (rho*physical path length).
Mixed Ni/H attenuation is supported but mixed mu_en is intentionally unavailable.

Charged carrier: kind electron/proton/alpha, energy_mev, count, and transport:
{"medium_element":"Ni", "track_areal_g_cm2":0.1,
 "track_evidence":"explicit track, medium and state applicability assumption"}.
Only electron Ni/H and proton/alpha H have acquired stopping tables. Track mass
length is NOT geometric chord, projected range or the periodic 4^3 cell size.
Alternatively supply unresolved_reason instead of transport; this also supports
kind recoil/neutron/positron, whose kinetic transport is not available here.
Neutrino: kind="neutrino", energy_mev, count, free_streaming=true (explicit).
All energies are positive MeV, counts nonnegative per normalization. No supplied
kinetic energy is silently supplemented by Q, rest mass or annihilation photons.
Neither mode reports actual material heat or power.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
from substrate_framework.nuclear_energy_transport import load_primary_tables, conditional_energy_ledger

HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE / 'sources/ni-energy-transport-original'


def binding(path):
    raw = Path(path).read_bytes()
    return {'file': str(Path(path).relative_to(ROOT)), 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}


def build_report(request=None):
    registration_path = HERE / 'nuclear-transport-registration.json'
    registration = json.loads(registration_path.read_text())
    acquisition = SOURCE_DIR / 'acquisition.json'
    source_binding = binding(acquisition)
    if source_binding['sha256'] != registration['acquisition_sha256'] or source_binding['bytes'] != registration['acquisition_bytes']:
        raise ValueError('Acquisition manifest no longer matches preregistration')
    photons, stopping, originals = load_primary_tables(SOURCE_DIR)
    gap_path = HERE / registration['supplemental_gap_receipt']['file']
    gap_binding = binding(gap_path)
    for key in ('sha256', 'bytes'):
        if gap_binding[key] != registration['supplemental_gap_receipt'][key]:
            raise ValueError('Supplemental availability receipt binding mismatch')
    gap_receipt = json.loads(gap_path.read_text())
    for source in gap_receipt['sources']:
        observed = binding(SOURCE_DIR / source['file'])
        if any(observed[key] != source[key] for key in ('sha256', 'bytes')):
            raise ValueError('Original empty-table response binding mismatch')
    report = {
        'schema': 'P254.A28.nuclear-energy-transport.v1',
        'mode': 'actual-primary-table-inventory' if request is None else 'explicit-conditional-energy-transfer',
        'registration': binding(registration_path), 'acquisition': source_binding,
        'implementation': binding(ROOT / 'src/substrate_framework/nuclear_energy_transport.py'),
        'consumer': binding(Path(__file__)), 'original_sources': originals,
        'photon_tables': {}, 'stopping_tables': [],
        'source_gaps': json.loads(acquisition.read_text())['gaps'],
        'supplemental_gap': {'binding': gap_binding, 'finding': gap_receipt['finding'], 'sources': gap_receipt['sources']},
        'primary_table_limitations': registration['primary_table_limitations'],
        'actual_material_heat_mev': None, 'actual_material_power_w': None,
        'missing_material_inputs': [
            'Exclusive emitted spectra/branch populations and matching actual reaction rates, not inclusive cross sections or AME Q alone.',
            'Measured finite sample geometry, source position/direction distributions and Ni/H composition/phase; periodic 4^3 lattice is not an absorber.',
            'Secondary photon and charged transport, containment and thermalization; compound-specific coefficients where needed.',
        ],
    }
    for element, table in photons.items():
        report['photon_tables'][element] = {
            'rows': len(table.rows), 'energy_domain_mev': [table.rows[0].energy_mev, table.rows[-1].energy_mev],
            'edges': [{'label': row.edge_label, 'energy_mev': row.energy_mev,
                       'below_raw': table.rows[index - 1].raw, 'above_raw': row.raw}
                      for index, row in enumerate(table.rows) if row.edge_label],
            'all_original_rows': [{'raw': row.raw, 'edge_label': row.edge_label} for row in table.rows],
        }
    for (particle, element), table in stopping.items():
        report['stopping_tables'].append({
            'particle': particle, 'element': element, 'rows': len(table.rows),
            'energy_domain_mev': [table.rows[0].energy_mev, table.rows[-1].energy_mev],
            'csda_range_domain_g_cm2': [table.rows[0].csda_g_cm2, table.rows[-1].csda_g_cm2],
            'columns': ['energy_MeV', 'collision_MeV_cm2_g',
                        'radiative_MeV_cm2_g' if particle == 'electron' else 'nuclear_MeV_cm2_g',
                        'total_MeV_cm2_g', 'CSDA_g_cm2',
                        'radiation_yield' if particle == 'electron' else 'projected_range_g_cm2',
                        'density_effect' if particle == 'electron' else 'detour_factor'],
            'all_original_rows': [row.raw for row in table.rows],
        })
    if request is not None:
        report['preparation'] = request
        report['conditional_ledger'] = conditional_energy_ledger(request, photons, stopping)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=Path, help='Explicit conditional carrier/path JSON; omitted means inventory only')
    parser.add_argument('--output', type=Path, help='Output JSON (otherwise stdout)')
    args = parser.parse_args()
    request = json.loads(args.input.read_text()) if args.input else None
    result = build_report(request)
    if args.input:
        raw = args.input.read_bytes()
        result['input_binding'] = {'file': str(args.input), 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
