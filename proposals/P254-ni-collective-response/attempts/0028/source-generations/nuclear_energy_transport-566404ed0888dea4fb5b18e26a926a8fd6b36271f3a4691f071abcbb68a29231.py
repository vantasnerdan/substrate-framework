"""Conditional energy transport from original NIST photon and STAR tables.

No reaction model, measured sample geometry, cascade spectrum or heat rate is
inferred. See A28 nuclear-transport-registration.json for the mathematical
construction, source receipts and the distinction between kerma and finite-body
heat. All energies are MeV; mass thicknesses and CSDA tracks are g/cm2.
"""
from __future__ import annotations

from dataclasses import dataclass
from html.parser import HTMLParser
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Mapping


class _PreText(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.active = False
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag == 'pre':
            self.active = True
        elif tag == 'br' and self.active:
            self.parts.append('\n')

    def handle_endtag(self, tag):
        if tag == 'pre':
            self.active = False
            self.parts.append('\n')

    def handle_data(self, data):
        if self.active:
            self.parts.append(data)


def _pre(path):
    parser = _PreText()
    parser.feed(Path(path).read_text(encoding='windows-1252'))
    return ''.join(parser.parts)


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} must be a numeric scalar')
    value = float(value)
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f'{name} must be finite and {"positive" if positive else "nonnegative"}')
    return value


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'Explicit {name} is required')
    return value


def _keys(value, required, optional=()):
    if not isinstance(value, dict) or set(value) - set(required) - set(optional) or set(required) - set(value):
        raise ValueError(f'Expected keys {sorted(required)}, optional {sorted(optional)}')


@dataclass(frozen=True)
class PhotonRow:
    energy_mev: float
    attenuation_cm2_g: float
    energy_absorption_cm2_g: float
    edge_label: str
    raw: tuple[str, str, str]


@dataclass(frozen=True)
class PhotonTable:
    element: str
    rows: tuple[PhotonRow, ...]

    def coefficients(self, energy_mev, side=None):
        energy = _number(energy_mev, 'photon energy MeV', positive=True)
        if side not in (None, 'below', 'above'):
            raise ValueError('edge side must be below or above')
        exact = [row for row in self.rows if row.energy_mev == energy]
        if exact:
            if len(exact) == 2 and side is None:
                raise ValueError(f'{self.element} absorption edge requires below/above side')
            row = exact[-1] if side == 'above' else exact[0]
            return row.attenuation_cm2_g, row.energy_absorption_cm2_g
        if not self.rows[0].energy_mev < energy < self.rows[-1].energy_mev:
            raise ValueError(f'{self.element} photon energy outside original source domain')
        # Last lower row is above any lower edge; first upper row is below it.
        lower = next(row for row in reversed(self.rows) if row.energy_mev < energy)
        upper = next(row for row in self.rows if row.energy_mev > energy)
        t = math.log(energy / lower.energy_mev) / math.log(upper.energy_mev / lower.energy_mev)
        return tuple(_log_value(a, b, t) for a, b in (
            (lower.attenuation_cm2_g, upper.attenuation_cm2_g),
            (lower.energy_absorption_cm2_g, upper.energy_absorption_cm2_g),
        ))


def _log_value(a, b, t):
    return math.exp(math.log(a) + t * math.log(b / a))


def read_photon_table(path, element):
    if element not in ('Ni', 'H'):
        raise ValueError('Only acquired Ni/H photon tables are supported')
    raw_html = Path(path).read_text(encoding='windows-1252')
    name = 'Nickel' if element == 'Ni' else 'Hydrogen'
    if name not in raw_html:
        raise ValueError('Wrong photon source material')
    rows = []
    for line in _pre(path).splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[0] in ('L1', 'K'):
            edge, fields = fields[0], fields[1:]
        else:
            edge = ''
        if len(fields) != 3 or not all(re.fullmatch(r'\d+\.\d+E[+-]\d+', item) for item in fields):
            continue
        energy, mu, mu_en = map(float, fields)
        if not 0 < mu_en <= mu or energy <= 0:
            raise ValueError('Invalid photon source row')
        if rows and (energy < rows[-1].energy_mev or (energy == rows[-1].energy_mev) != bool(edge)):
            raise ValueError('Invalid source edge order')
        rows.append(PhotonRow(energy, mu, mu_en, edge, tuple(fields)))
    if len(rows) < 2:
        raise ValueError('No original photon table found')
    return PhotonTable(element, tuple(rows))


@dataclass(frozen=True)
class StoppingRow:
    energy_mev: float
    collision_mev_cm2_g: float
    other_mev_cm2_g: float
    total_mev_cm2_g: float
    csda_g_cm2: float
    extra: float  # Electron radiation yield OR ion projected range.
    final: float  # Electron density effect OR ion detour factor.
    raw: tuple[str, ...]


@dataclass(frozen=True)
class StoppingTable:
    particle: str
    element: str
    rows: tuple[StoppingRow, ...]

    def value(self, energy_mev, field):
        energy = _number(energy_mev, 'charged energy MeV', positive=True)
        for row in self.rows:
            if energy == row.energy_mev:
                return getattr(row, field)
        if not self.rows[0].energy_mev < energy < self.rows[-1].energy_mev:
            raise ValueError('Charged energy outside original source domain')
        lo = next(row for row in reversed(self.rows) if row.energy_mev < energy)
        hi = next(row for row in self.rows if row.energy_mev > energy)
        t = math.log(energy / lo.energy_mev) / math.log(hi.energy_mev / lo.energy_mev)
        return _log_value(getattr(lo, field), getattr(hi, field), t)

    def energy_at_range(self, residual_range):
        if residual_range == 0:
            return 0.0  # Complete stop, not a low-energy extrapolated point.
        for row in self.rows:
            if residual_range == row.csda_g_cm2:
                return row.energy_mev
        if not self.rows[0].csda_g_cm2 < residual_range < self.rows[-1].csda_g_cm2:
            raise ValueError('Residual CSDA range falls outside table; no sub-domain extrapolation')
        lo = next(row for row in reversed(self.rows) if row.csda_g_cm2 < residual_range)
        hi = next(row for row in self.rows if row.csda_g_cm2 > residual_range)
        t = math.log(residual_range / lo.csda_g_cm2) / math.log(hi.csda_g_cm2 / lo.csda_g_cm2)
        return _log_value(lo.energy_mev, hi.energy_mev, t)


def read_stopping_table(path, particle, element):
    text = _pre(path)
    title = {'electron': 'Electrons in', 'proton': 'PROTON STOPPING', 'alpha': 'HELIUM ION STOPPING'}
    if particle not in title or element not in ('Ni', 'H'):
        raise ValueError('Unsupported STAR material or particle')
    # ASTAR header uses HELIUM IONS; require its distinctive range terminology.
    if particle == 'alpha':
        valid = 'HELIUM' in text and 'STOPPING' in text
    else:
        valid = title[particle] in text
    if not valid or ('NICKEL' if element == 'Ni' else 'HYDROGEN') not in text:
        raise ValueError('Wrong STAR particle or material')
    rows = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 7:
            continue
        try:
            numbers = tuple(float(item) for item in fields)
        except ValueError:
            continue
        if not all(math.isfinite(x) and x >= 0 for x in numbers) or min(numbers[:5]) <= 0:
            raise ValueError('Invalid STAR row')
        if abs(numbers[1] + numbers[2] - numbers[3]) > 0.002 * numbers[3]:
            raise ValueError('STAR stopping columns do not sum within printed rounding')
        if rows and (numbers[0] <= rows[-1].energy_mev or numbers[4] <= rows[-1].csda_g_cm2):
            raise ValueError('Nonmonotone STAR energy/range')
        rows.append(StoppingRow(*numbers, tuple(fields)))
    if len(rows) < 2:
        raise ValueError('No original STAR table found')
    return StoppingTable(particle, element, tuple(rows))


def load_primary_tables(source_dir):
    """Bind every original byte, then parse originals directly (no derived fit)."""
    directory = Path(source_dir)
    receipt = json.loads((directory / 'acquisition.json').read_text())
    bindings = []
    for source in receipt['sources']:
        if source.get('file') is None:
            continue
        data = (directory / source['file']).read_bytes()
        if len(data) != source['bytes'] or hashlib.sha256(data).hexdigest() != source['sha256']:
            raise ValueError(f'Original source binding mismatch: {source["file"]}')
        bindings.append({key: source[key] for key in ('file', 'sha256', 'bytes', 'status')})
    photons = {element: read_photon_table(directory / filename, element) for element, filename in (
        ('Ni', 'nist-nickel-z28.html'), ('H', 'nist-hydrogen-z01.html'))}
    stopping = {}
    for particle, element, name in (
        ('electron', 'Ni', 'nist-estar-nickel-response.html'),
        ('electron', 'H', 'nist-estar-hydrogen-response.html'),
        ('proton', 'H', 'nist-pstar-hydrogen-response.html'),
        ('alpha', 'H', 'nist-astar-hydrogen-response.html'),
    ):
        stopping[particle, element] = read_stopping_table(directory / name, particle, element)
    return photons, stopping, bindings


def _composition(fractions):
    if not isinstance(fractions, dict) or not fractions or set(fractions) - {'Ni', 'H'}:
        raise ValueError('Explicit Ni/H mass fractions required')
    fractions = {key: _number(value, 'mass fraction') for key, value in fractions.items()}
    if not math.isclose(math.fsum(fractions.values()), 1.0, rel_tol=0, abs_tol=1e-12):
        raise ValueError('Mass fractions must sum to one, not atom fractions')
    return {key: value for key, value in fractions.items() if value > 0}


def photon_transfer(tables: Mapping[str, PhotonTable], energy_mev, paths, edge_side=None):
    """Primary uncollided escape plus kerma/unresolved partition of removed energy.

    This is NOT exp(-mu_en*X) absorption and NOT a scattered-photon solver.
    Homogeneous mixtures intentionally have no fabricated mu_en coefficient.
    """
    energy = _number(energy_mev, 'photon energy MeV', positive=True)
    if not isinstance(paths, list) or not paths:
        raise ValueError('Explicit nonempty photon path ensemble required')
    escape = kerma = unresolved = 0.0
    probability_sum = 0.0
    details = []
    for path in paths:
        _keys(path, ('probability', 'geometry_evidence', 'layers'))
        _text(path['geometry_evidence'], 'physical path geometry evidence')
        probability = _number(path['probability'], 'path probability')
        probability_sum += probability
        if not isinstance(path['layers'], list):
            raise ValueError('layers must be an ordered list; [] explicitly means vacuum')
        remaining = energy
        path_kerma = path_unresolved = 0.0
        layers = []
        for layer in path['layers']:
            _keys(layer, ('areal_density_g_cm2', 'mass_fractions'))
            thickness = _number(layer['areal_density_g_cm2'], 'areal density g/cm2')
            fractions = _composition(layer['mass_fractions'])
            coefficients = {element: tables[element].coefficients(energy, edge_side) for element in fractions}
            mu = math.fsum(fractions[element] * coefficients[element][0] for element in fractions)
            mu_en = next(iter(coefficients.values()))[1] if len(fractions) == 1 else None
            removed = remaining * -math.expm1(-mu * thickness)
            allocated = removed * mu_en / mu if mu_en is not None else 0.0
            path_kerma += allocated
            path_unresolved += removed - allocated
            remaining *= math.exp(-mu * thickness)
            layers.append({'attenuation_cm2_g': mu, 'energy_absorption_cm2_g': mu_en,
                           'primary_removed_mev': removed, 'collision_kerma_budget_mev': allocated,
                           'absorption_status': 'pure-element source' if mu_en is not None else 'mixture mu_en unavailable'})
        escape += probability * remaining
        kerma += probability * path_kerma
        unresolved += probability * path_unresolved
        details.append({'probability': probability, 'layers': layers})
    if not math.isclose(probability_sum, 1, rel_tol=0, abs_tol=1e-12):
        raise ValueError('Path probabilities must sum to one')
    return {'photon_uncollided_escape': escape, 'photon_collision_kerma_budget': kerma,
            'photon_removed_unresolved': unresolved}, details


def charged_slowing(table, energy_mev, track_areal_g_cm2):
    """CSDA kinetic budget on a specified TRACK, not geometric chord/depth.

    The loss is not a local heat measurement. Electronic/nuclear stopping ratios
    at the initial energy are diagnostics only. Full-stop ranges/yields include
    NIST's own integration below the lowest output energy; we never interpolate
    a residual energy there. Electrons mean electrons, not positrons.
    """
    energy = _number(energy_mev, 'charged energy MeV', positive=True)
    track = _number(track_areal_g_cm2, 'charged track g/cm2')
    initial_range = table.value(energy, 'csda_g_cm2')
    if track == 0:
        residual = energy
    else:
        residual = table.energy_at_range(max(initial_range - track, 0.0))
    loss = energy - residual
    radiation = 0.0
    if table.particle == 'electron':
        radiation = energy * table.value(energy, 'extra')
        if residual:
            radiation -= residual * table.value(residual, 'extra')
        if radiation < -1e-12 or radiation > loss + 1e-12:
            raise ValueError('Interpolated electron radiation yield violates kinetic budget')
        radiation = min(loss, max(0.0, radiation))
    return {
        'charged_residual_kinetic': residual,
        'electron_bremsstrahlung_untransported': radiation,
        'electron_collision_slowing' if table.particle == 'electron' else 'ion_electronic_plus_recoil_slowing': loss - radiation,
    }, {'initial_csda_g_cm2': initial_range, 'used_track_g_cm2': min(track, initial_range),
        'initial_collision_stopping_mev_cm2_g': table.value(energy, 'collision_mev_cm2_g'),
        'initial_radiative_stopping_mev_cm2_g' if table.particle == 'electron' else 'initial_nuclear_recoil_stopping_mev_cm2_g': table.value(energy, 'other_mev_cm2_g'),
        'residual_status': 'at supplied track endpoint; not proven escaping' if residual else 'CSDA stop; not proof of finite-sample containment'}


def conditional_energy_ledger(request, photons, stopping):
    """Account supplied line/quadrature nodes; no implicit spectral integration.

    Each count is the explicit expected number in the stated normalization.
    Neutrino, recoil, neutron, positron and unresolved charged energies must be
    supplied separately; an optional expected energy checks only that supplied
    budget, never infers a spectrum from Q.
    """
    _keys(request, ('normalization', 'spectrum_evidence', 'carriers'), ('expected_emitted_energy_mev',))
    _text(request['normalization'], 'normalization (not an inferred event rate)')
    _text(request['spectrum_evidence'], 'spectrum provenance or explicit conditional assumption')
    if not isinstance(request['carriers'], list) or not request['carriers']:
        raise ValueError('Explicit emitted carriers required')
    budgets = {}
    details = []
    input_terms = []
    for carrier in request['carriers']:
        _keys(carrier, ('kind', 'energy_mev', 'count'), ('paths', 'edge_side', 'transport', 'free_streaming', 'unresolved_reason'))
        kind = carrier['kind']
        energy = _number(carrier['energy_mev'], 'carrier energy MeV', positive=True)
        count = _number(carrier['count'], 'carrier count')
        input_terms.append(energy * count)
        if kind == 'photon':
            if set(carrier) - {'kind', 'energy_mev', 'count', 'paths', 'edge_side'}:
                raise ValueError('Photon carrier has incompatible fields')
            terms, detail = photon_transfer(photons, energy, carrier.get('paths'), carrier.get('edge_side'))
        elif kind == 'neutrino':
            if set(carrier) != {'kind', 'energy_mev', 'count', 'free_streaming'} or carrier['free_streaming'] is not True:
                raise ValueError('Neutrino escape requires explicit free_streaming=true')
            terms, detail = {'neutrino_escape': energy}, {}
        elif kind in ('electron', 'proton', 'alpha', 'recoil', 'neutron', 'positron'):
            if 'transport' in carrier:
                if set(carrier) != {'kind', 'energy_mev', 'count', 'transport'}:
                    raise ValueError('Charged transport has incompatible fields')
                transport = carrier['transport']
                _keys(transport, ('medium_element', 'track_areal_g_cm2', 'track_evidence'))
                _text(transport['track_evidence'], 'charged track and medium applicability evidence')
                key = kind, transport['medium_element']
                if key not in stopping:
                    raise ValueError(f'No acquired primary stopping table for {key}; explicitly keep unresolved')
                terms, detail = charged_slowing(stopping[key], energy, transport['track_areal_g_cm2'])
            else:
                if set(carrier) != {'kind', 'energy_mev', 'count', 'unresolved_reason'}:
                    raise ValueError('Unavailable carrier transport needs explicit unresolved_reason')
                _text(carrier['unresolved_reason'], 'unresolved reason')
                terms, detail = {kind + '_kinetic_unresolved': energy}, {'reason': carrier['unresolved_reason']}
        else:
            raise ValueError('Unsupported energy carrier kind')
        for name, value in terms.items():
            budgets[name] = budgets.get(name, 0.0) + count * value
        details.append({'kind': kind, 'energy_mev': energy, 'count': count, 'per_particle_budget_mev': terms, 'transport': detail})
    emitted = math.fsum(input_terms)
    accounted = math.fsum(budgets.values())
    if not math.isfinite(emitted) or not math.isfinite(accounted):
        raise ValueError('Energy budget overflow')
    tolerance = 2e-12 * max(emitted, 1e-300)
    if abs(emitted - accounted) > tolerance:
        raise ValueError('Conditional transfer failed energy conservation')
    if 'expected_emitted_energy_mev' in request:
        expected = _number(request['expected_emitted_energy_mev'], 'expected emitted energy')
        if not math.isclose(emitted, expected, rel_tol=2e-12, abs_tol=0):
            raise ValueError('Explicit carrier budget does not match supplied expected energy')
    return {'normalization': request['normalization'], 'spectrum_evidence': request['spectrum_evidence'],
            'emitted_energy_mev': emitted, 'accounted_energy_mev': accounted,
            'conservation_residual_mev': emitted - accounted, 'budget_mev': budgets,
            'carriers': details, 'actual_material_heat_mev': None, 'actual_material_power_w': None,
            'limitations': [
                'Collision kerma is a primary-field budget, not finite-sample absorbed dose or heat; secondary charged escape and secondary photon transport are unresolved.',
                'Only uncollided photons are proven escaping along the supplied open paths; removed radiation is not synonymous with absorbed or escaped energy.',
                'CSDA track length is not projected range, escape probability, microscopic recoil spectrum or thermalization time. Residual charged kinetic energy is not automatically escape.',
                'Elemental Ni/H attenuation is supported; mixed-medium mu_en and mixed-medium stopping are not fabricated.',
                'Source photon tables omit photonuclear interactions and edge chemical/phase fine structure. Table interpolation is a declared approximation, not the native NIST interpolation.',
                'Hydrogen stopping tables describe the source hydrogen medium, not interstitial H in Ni; applying any table requires matching medium/state evidence.',
                'Input spectrum/counts/geometry are external conditions; no cascade or rate comes from AME Q, inclusive cross sections, or the periodic 4^3 cell.',
                'Energy ledger counts supplied kinetic/radiative energy only. Positron annihilation/rest-energy accounting requires separately supplied radiation and endpoint convention.',
            ]}
