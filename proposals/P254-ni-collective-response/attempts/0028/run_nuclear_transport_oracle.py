#!/usr/bin/env python3
"""Independent original-table, Decimal-transfer and CSDA endpoint oracle.

For Root to execute once after integration. Conditional fixtures below are
explicit mathematical preparations, NOT observed nuclear spectra, sample sizes,
material rates, or heat measurements. No cloud or material trajectory is used.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from decimal import Decimal, localcontext
import hashlib
import html
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
from substrate_framework.nuclear_energy_transport import (
    load_primary_tables, photon_transfer, charged_slowing, conditional_energy_ledger,
)
from consume_nuclear_energy_transport import SOURCE_DIR, build_report, binding


def close(actual, expected, rtol=3e-12, atol=1e-15):
    assert math.isclose(actual, float(expected), rel_tol=rtol, abs_tol=atol), (actual, expected)


def reject(call):
    try:
        call()
    except (ValueError, TypeError, AssertionError):
        return
    raise AssertionError('Sensitive mutation survived')


def independent_photon_rows(path):
    # Use HTML table cells, NOT the production ASCII PRE text.
    source = path.read_text(encoding='windows-1252')
    rows = []
    for block in re.findall(r'<TR\b[^>]*>(.*?)</TR>', source, flags=re.I | re.S):
        fields = [html.unescape(re.sub(r'<[^>]*>', '', value)).strip()
                  for value in re.findall(r'<TD\b[^>]*>(.*?)</TD>', block, flags=re.I | re.S)]
        fields = [value for value in fields if value]
        if len(fields) == 4 and fields[0] in ('K', 'L1'):
            edge, fields = fields[0], fields[1:]
        else:
            edge = ''
        if len(fields) == 3 and all(re.fullmatch(r'\d+\.\d+E[+-]\d+', value) for value in fields):
            rows.append((tuple(fields), edge))
    return rows


def independent_star_rows(path):
    # Direct original BR-delimited numeric rows; independent of HTMLParser,
    # production table objects, column decoding, or any derived table JSON.
    source = path.read_text(encoding='windows-1252')
    values = []
    for fragment in re.split(r'<br\s*/?>', source, flags=re.I):
        match = re.fullmatch(r'\s*((?:[0-9.E+\-]+\s+){6}[0-9.E+\-]+)\s*', fragment)
        if match:
            fields = tuple(match.group(1).split())
            tuple(Decimal(value) for value in fields)
            values.append(fields)
    return values


def source_oracle(photons, stopping):
    for element, filename in [('Ni', 'nist-nickel-z28.html'), ('H', 'nist-hydrogen-z01.html')]:
        independent = independent_photon_rows(SOURCE_DIR / filename)
        assert independent == [(row.raw, row.edge_label) for row in photons[element].rows]
        assert len(independent) == (41 if element == 'Ni' else 36)
    ni = photons['Ni']
    assert ni.coefficients(0.0083328, 'below') == (44.28, 42.42)
    assert ni.coefficients(0.0083328, 'above') == (329.4, 224.0)
    assert ni.coefficients(0.0010081, 'below') == (9654.0, 9598.0)
    assert ni.coefficients(0.0010081, 'above') == (10990.0, 10930.0)
    assert photons['H'].coefficients(1.0) == (0.1263, 0.05556)
    assert ni.coefficients(1.0) == (0.06160, 0.02674)
    for particle, element, filename in [
        ('proton', 'H', 'nist-pstar-hydrogen-response.html'),
        ('alpha', 'H', 'nist-astar-hydrogen-response.html'),
        ('electron', 'Ni', 'nist-estar-nickel-response.html'),
        ('electron', 'H', 'nist-estar-hydrogen-response.html'),
    ]:
        independent = independent_star_rows(SOURCE_DIR / filename)
        assert independent == [row.raw for row in stopping[particle, element].rows]
        assert len(independent) > 50
    assert stopping['proton', 'H'].value(1.0, 'csda_g_cm2') == 8.476e-4
    assert stopping['alpha', 'H'].value(1.0, 'csda_g_cm2') == 1.671e-4
    assert stopping['electron', 'Ni'].value(1.0, 'csda_g_cm2') == 0.6087
    assert stopping['electron', 'Ni'].value(1.0, 'extra') == 0.01864
    assert stopping['electron', 'H'].value(1.0, 'csda_g_cm2') == 0.2117


def path(layers, probability=1.0):
    return {'probability': probability, 'geometry_evidence': 'oracle-only supplied external ray, not a periodic cell', 'layers': layers}


def layer(element, thickness):
    return {'areal_density_g_cm2': thickness, 'mass_fractions': {element: 1.0}}


def decimal_ray(energy, layers):
    # Literal original 1 MeV coefficients, not production interpolation/data.
    coefficients = {'Ni': (Decimal('0.06160'), Decimal('0.02674')),
                    'H': (Decimal('0.1263'), Decimal('0.05556'))}
    with localcontext() as context:
        context.prec = 60
        remaining = Decimal(str(energy))
        kerma = Decimal(0)
        unresolved = Decimal(0)
        for element, thickness in layers:
            mu, mu_en = coefficients[element]
            survival = (-mu * Decimal(str(thickness))).exp()
            removed = remaining * (1 - survival)
            kerma += removed * mu_en / mu
            unresolved += removed * (1 - mu_en / mu)
            remaining *= survival
        return tuple(float(value) for value in (remaining, kerma, unresolved))


def photon_oracle(photons):
    paths = [path([layer('Ni', 2.0), layer('H', 0.5)], 0.3), path([], 0.7)]
    result, detail = photon_transfer(photons, 1.0, paths)
    expected = decimal_ray(1.0, [('Ni', 2.0), ('H', 0.5)])
    for key, target in zip(('photon_uncollided_escape', 'photon_collision_kerma_budget', 'photon_removed_unresolved'),
                           (0.7 + 0.3 * expected[0], 0.3 * expected[1], 0.3 * expected[2])):
        close(result[key], target)
    close(sum(result.values()), 1.0)
    reverse, _ = photon_transfer(photons, 1.0, [path([layer('H', 0.5), layer('Ni', 2.0)])])
    close(reverse['photon_uncollided_escape'], expected[0])
    reject(lambda: close(reverse['photon_collision_kerma_budget'], expected[1]))
    # Distinguish attenuation from absorption, and reject exp(-mu_en*X).
    reject(lambda: close(result['photon_collision_kerma_budget'], 1-result['photon_uncollided_escape']))
    pure, _ = photon_transfer(photons, 1.0, [path([layer('Ni', 2.0)])])
    reject(lambda: close(pure['photon_collision_kerma_budget'], -math.expm1(-0.02674*2)))
    mixed = {'areal_density_g_cm2': 2, 'mass_fractions': {'Ni': 0.75, 'H': 0.25}}
    mixture, details = photon_transfer(photons, 1.0, [path([mixed])])
    close(mixture['photon_uncollided_escape'], math.exp(-2*(0.75*0.06160+0.25*0.1263)))
    assert mixture['photon_collision_kerma_budget'] == 0
    assert details[0]['layers'][0]['energy_absorption_cm2_g'] is None
    close(mixture['photon_removed_unresolved'], 1-mixture['photon_uncollided_escape'])
    for thickness in (0.0, 1e-16, 1e6):
        values, _ = photon_transfer(photons, 1.0, [path([layer('Ni', thickness)])])
        close(sum(values.values()), 1.0)
        assert min(values.values()) >= 0
    tiny, _ = photon_transfer(photons, 1.0, [path([layer('Ni', 1e-16)])])
    close(tiny['photon_collision_kerma_budget'], 0.02674e-16, atol=1e-32)
    ni = photons['Ni']
    close(ni.coefficients(math.sqrt(0.008*0.0083328))[0], math.sqrt(49.52*44.28))
    close(ni.coefficients(math.sqrt(0.0083328*0.01))[0], math.sqrt(329.4*209.0))
    reject(lambda: ni.coefficients(0.0083328))
    reject(lambda: ni.coefficients(0.0010081))
    reject(lambda: ni.coefficients(1e-6))
    reject(lambda: ni.coefficients(21))
    reject(lambda: ni.coefficients(1000))  # keV supplied as MeV.
    reject(lambda: ni.coefficients(float('nan')))
    reject(lambda: photon_transfer(photons, 1.0, []))
    reject(lambda: photon_transfer(photons, 1.0, [path([], 0.8)]))
    reject(lambda: photon_transfer(photons, 1.0, [path([layer('Ni', -1)])]))
    bad = deepcopy(mixed)
    bad['mass_fractions'] = {'Ni': 1, 'H': 1}
    reject(lambda: photon_transfer(photons, 1.0, [path([bad])]))
    missing = path([])
    del missing['geometry_evidence']
    reject(lambda: photon_transfer(photons, 1.0, [missing]))
    changed_rows = tuple(replace(row, energy_absorption_cm2_g=row.attenuation_cm2_g) for row in ni.rows)
    changed = dict(photons, Ni=replace(ni, rows=changed_rows))
    mutated, _ = photon_transfer(changed, 1.0, [path([layer('Ni', 2)])])
    reject(lambda: close(mutated['photon_collision_kerma_budget'], pure['photon_collision_kerma_budget']))


def csda_oracle(stopping):
    electron = stopping['electron', 'Ni']
    budget, diagnostic = charged_slowing(electron, 2.0, 1.329-0.6087)
    close(budget['charged_residual_kinetic'], 1.0)
    expected_brems = 2*0.03384-0.01864
    close(budget['electron_bremsstrahlung_untransported'], expected_brems)
    close(budget['electron_collision_slowing'], 1-expected_brems)
    close(sum(budget.values()), 2.0)
    reject(lambda: close(budget['electron_collision_slowing'], 1.0))
    full, _ = charged_slowing(electron, 2.0, 1.329)
    assert full['charged_residual_kinetic'] == 0
    close(full['electron_bremsstrahlung_untransported'], 2*0.03384)
    beyond, _ = charged_slowing(electron, 2.0, 10)
    assert beyond == full
    vacuum, _ = charged_slowing(electron, 2.0, 0)
    assert vacuum == {'charged_residual_kinetic': 2.0, 'electron_bremsstrahlung_untransported': 0.0, 'electron_collision_slowing': 0.0}
    # Every adjacent table interval independently checks the declared log rule.
    for table in stopping.values():
        for a, b in zip(table.rows, table.rows[1:]):
            middle_energy = math.sqrt(a.energy_mev*b.energy_mev)
            middle_range = math.sqrt(a.csda_g_cm2*b.csda_g_cm2)
            close(table.value(middle_energy, 'csda_g_cm2'), middle_range)
            close(table.energy_at_range(middle_range), middle_energy)
        lowest = table.rows[0]
        reject(lambda: charged_slowing(table, lowest.energy_mev, lowest.csda_g_cm2/2))
        reject(lambda: charged_slowing(table, lowest.energy_mev/2, 0))
        reject(lambda: charged_slowing(table, lowest.energy_mev, -1))
    proton = stopping['proton', 'H']
    p, diagnostic = charged_slowing(proton, 2.0, 0.002883-0.0008476)
    close(p['charged_residual_kinetic'], 1.0)
    close(p['ion_electronic_plus_recoil_slowing'], 1.0)
    assert diagnostic['initial_nuclear_recoil_stopping_mev_cm2_g'] == 0.3772
    projected_as_csda = replace(proton, rows=tuple(replace(row, csda_g_cm2=row.extra) for row in proton.rows))
    wrong, _ = charged_slowing(projected_as_csda, 2.0, 0.002883-0.0008476)
    reject(lambda: close(wrong['charged_residual_kinetic'], 1.0))


def ledger_oracle(photons, stopping):
    request = {'normalization': 'one oracle-only conditional energy preparation',
               'spectrum_evidence': 'explicit oracle nodes, not an observed reaction branch',
               'expected_emitted_energy_mev': 5.2,
               'carriers': [
                   {'kind': 'photon', 'energy_mev': 1.0, 'count': 2, 'paths': [path([layer('Ni', 2)])]},
                   {'kind': 'electron', 'energy_mev': 2.0, 'count': 1,
                    'transport': {'medium_element': 'Ni', 'track_areal_g_cm2': 1.329,
                                  'track_evidence': 'oracle-only complete CSDA track in pure Ni'}},
                   {'kind': 'neutrino', 'energy_mev': 1.0, 'count': 1, 'free_streaming': True},
                   {'kind': 'recoil', 'energy_mev': 0.1, 'count': 2,
                    'unresolved_reason': 'No acquired heavy-ion table; no invented recoil stopping'},
               ]}
    report = conditional_energy_ledger(request, photons, stopping)
    close(report['emitted_energy_mev'], 5.2)
    close(report['accounted_energy_mev'], 5.2)
    close(report['budget_mev']['neutrino_escape'], 1)
    close(report['budget_mev']['recoil_kinetic_unresolved'], 0.2)
    close(report['budget_mev']['photon_uncollided_escape'], 2*decimal_ray(1, [('Ni', 2)])[0])
    assert report['actual_material_heat_mev'] is None and report['actual_material_power_w'] is None
    for mutate in (
        lambda r: r.update(expected_emitted_energy_mev=10),
        lambda r: r['carriers'][0].update(count=-1),
        lambda r: r['carriers'][0].update(energy_mev=float('inf')),
        lambda r: r['carriers'][1].update(kind='positron'),
        lambda r: r['carriers'][1].update(kind='proton'),
        lambda r: r['carriers'][1]['transport'].update(medium_element='NiH'),
        lambda r: r['carriers'][2].update(free_streaming=False),
        lambda r: r['carriers'][3].pop('unresolved_reason'),
        lambda r: r.update(periodic_cell_as_sample=True),
    ):
        changed = deepcopy(request)
        mutate(changed)
        reject(lambda: conditional_energy_ledger(changed, photons, stopping))
    return request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    photons, stopping, bindings = load_primary_tables(SOURCE_DIR)
    source_oracle(photons, stopping)
    photon_oracle(photons)
    csda_oracle(stopping)
    preparation = ledger_oracle(photons, stopping)
    actual = build_report()
    assert actual['mode'] == 'actual-primary-table-inventory'
    assert actual['actual_material_heat_mev'] is None
    conditional = build_report(preparation)
    close(conditional['conditional_ledger']['accounted_energy_mev'], 5.2)
    result = {'schema': 'P254.A28.nuclear-transport-oracle.v1', 'status': 'passed',
              'oracle': binding(Path(__file__)), 'implementation': actual['implementation'],
              'registration': actual['registration'], 'acquisition': actual['acquisition'],
              'sources': bindings,
              'checks': ['independent-original-HTML-versus-PRE-rows', 'literal-edges-and-source-domains',
                         'Decimal-layered-path-primary-transfer', 'mixture-absorption-not-fabricated',
                         'stable-thin-thick-energy-conservation', 'edge-sided-interpolation',
                         'CSDA-range-endpoints-and-inversion', 'electron-radiative-yield-partition',
                         'projected-range-mutation', 'explicit-carrier-conservation-and-neutrino-escape',
                         'missing-geometry-spectrum-and-transport-rejections', 'actual-primary-table-consumer'],
              'scientific_scope': 'conditional mathematical oracle only; no actual material spectrum, heat, rate or sample geometry asserted'}
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
