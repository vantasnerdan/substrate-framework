"""Measured ordinary p-Ni channels, with no extrapolation or material rate default.

Data are an exact transcription of original Simon Table II and Lopez-Saavedra
v2 Table I. A response-matched flux calculation is conditional on caller-supplied
state/response evidence; it is NOT a reconstructed continuous excitation function.
Inclusive measurements do not specify amplitudes, exclusive branches or spectra.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from typing import Mapping

from .ni_h_nuclear_ledger import H1, HE4, Nuclide, Reaction

BINDINGS_SHA256 = '920f5f6f7e141e34644a3e5c0780447457bb8cf0569179b5c5c41fa9291b747b'
MB_TO_CM2 = Fraction(1, 10**27)


def exact_nonnegative(value: str | int | Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (str, int, Fraction)):
        raise TypeError('Supply exact decimal strings, integers or Fractions')
    result = Fraction(value)
    if result < 0:
        raise ValueError('Negative physical input')
    return result


def parenthetical(raw: str) -> tuple[Fraction, Fraction | None]:
    """Keep raw text elsewhere; interpret parenthetical last-digit uncertainty."""
    match = re.fullmatch(r'(\d+(?:\.\d+)?)(?:\((\d+)\))?', raw)
    if match is None:
        raise ValueError(f'Not a measured decimal: {raw!r}')
    central, digits = match.groups()
    places = len(central.split('.')[1]) if '.' in central else 0
    return Fraction(central), None if digits is None else Fraction(int(digits), 10**places)


def load_measured_tables(source_dir: str | Path) -> dict:
    """Verify original PDFs and derived text/transcription before use."""
    source_dir = Path(source_dir)
    raw = (source_dir / 'bindings.json').read_bytes()
    if hashlib.sha256(raw).hexdigest() != BINDINGS_SHA256:
        raise ValueError('Measured source binding manifest changed')
    bindings = json.loads(raw)
    for name, digest in {**bindings['originals'], **bindings['derived_files']}.items():
        if hashlib.sha256((source_dir / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Measured source hash mismatch: {name}')
    return json.loads((source_dir / 'measured-tables.json').read_text())


def measured_reaction(row: Mapping) -> Reaction:
    target = Nuclide(row['target_a'], row['target_z'])
    product = Nuclide(row['product_a'], row['product_z'])
    if target.z == 28:
        reaction = Reaction(row['id'], 'ordinary binary radiative capture',
                            (target, H1), (product,), ('gamma cascade',))
    elif target == Nuclide(59, 29):
        reaction = Reaction(row['id'], 'ordinary binary proton alpha reaction',
                            (target, H1), (product, HE4))
    else:
        raise ValueError('Unsupported measured channel')
    reaction.check_conservation()
    return reaction


@dataclass(frozen=True)
class FluxComponent:
    """Integrated flux assigned to one experimental response, not per MeV.

    Cu strip bounds must equal the published effective energy +/- bin offsets.
    Simon does not tabulate a response interval: both bounds must be None and
    the encounter model must independently justify matching that target response.
    """
    row_id: str
    energy_cm_mev: str
    support_lower_mev: str | None
    support_upper_mev: str | None
    flux: str
    energy_unit: str
    flux_unit: str


@dataclass(frozen=True)
class ProtonEncounterModel:
    components: tuple[FluxComponent, ...]
    entrance_state: str
    state_matching_justification: str
    response_matching_justification: str
    provenance: str


@dataclass(frozen=True)
class ConditionalRate:
    per_target_per_second: Fraction
    row_ids: tuple[str, ...]
    encounter_model: ProtonEncounterModel
    interpretation: str = ('Conditional central-value rate only; measurement uncertainties '
                           'and missing cross-row covariance are not a rate confidence interval')


def rate_from_matched_flux(rows: list[dict], model: ProtonEncounterModel) -> ConditionalRate:
    """Sum sigma_i Phi_i only for explicitly matched measured responses.

    No interpolation even inside the quoted energy range. Upper-limit rows may
    be published but cannot be promoted to central-value reaction hazards.
    Energies are relative CM energies; a cold-lattice energy is not a beam input.
    This checks declarations/support, not the truth of their physical evidence.
    """
    if model.entrance_state != 'laboratory-ground-state-unpolarized':
        raise ValueError('No polarized/excited-state cross sections acquired')
    if not all(x.strip() for x in (model.state_matching_justification,
                                    model.response_matching_justification, model.provenance)):
        raise ValueError('Explicit state, response and encounter provenance required')
    if not model.components:
        raise ValueError('Supply an explicit flux model (zero flux may be explicit)')
    index = {row['id']: row for row in rows}
    if len(index) != len(rows):
        raise ValueError('Duplicate measured row')
    ids = [component.row_id for component in model.components]
    if len(set(ids)) != len(ids):
        raise ValueError('Repeated response would double-count flux')
    targets = set()
    total = Fraction(0)
    for component in model.components:
        if component.energy_unit != 'MeV_CM' or component.flux_unit != 'cm^-2 s^-1':
            raise ValueError('Require MeV_CM and integrated cm^-2 s^-1 flux')
        if component.row_id not in index:
            raise ValueError('No acquired measurement for this response')
        row = index[component.row_id]
        targets.add((row['target_a'], row['target_z']))
        energy, _ = parenthetical(row['energy_cm_mev_raw'])
        if exact_nonnegative(component.energy_cm_mev) != energy:
            raise ValueError('Energy differs from source effective energy; no interpolation')
        if 'energy_bin_plus_mev_raw' in row:
            if component.support_lower_mev is None or component.support_upper_mev is None:
                raise ValueError('Finite measured-strip support must be supplied')
            expected = (energy - Fraction(row['energy_bin_minus_mev_raw']),
                        energy + Fraction(row['energy_bin_plus_mev_raw']))
            supplied = (exact_nonnegative(component.support_lower_mev),
                        exact_nonnegative(component.support_upper_mev))
            if supplied != expected:
                raise ValueError('Flux response does not match published finite strip domain')
        elif component.support_lower_mev is not None or component.support_upper_mev is not None:
            raise ValueError('Simon Table II supplies no bin bounds; do not invent them')
        if row['cross_section_mb_raw'].startswith('<'):
            raise ValueError('Upper limit is not a detected central-value cross section')
        sigma, _ = parenthetical(row['cross_section_mb_raw'])
        total += sigma * MB_TO_CM2 * exact_nonnegative(component.flux)
    if len(targets) != 1:
        raise ValueError('Per-target rate cannot combine different target isotopes')
    return ConditionalRate(total, tuple(ids), model)


def channel_domains(tables: dict) -> list[dict]:
    """Ranges summarize support, NOT a license to interpolate between rows."""
    groups: dict[tuple[int, int], list[dict]] = {}
    for row in tables['capture'] + tables['cu59_palpha']:
        groups.setdefault((row['target_a'], row['target_z']), []).append(row)
    result = []
    for (a, z), rows in sorted(groups.items()):
        energies = [parenthetical(row['energy_cm_mev_raw'])[0] for row in rows]
        result.append({'target': {'a': a, 'z': z}, 'number_of_rows': len(rows),
                       'minimum_effective_energy_mev': str(min(energies)),
                       'maximum_effective_energy_mev': str(max(energies)),
                       'support': 'Only the individual acquired experimental responses; no continuous kernel',
                       'upper_limit_rows': [r['id'] for r in rows if r['cross_section_mb_raw'].startswith('<')],
                       'entrance_state': rows[0]['entrance_state'],
                       'observable': rows[0]['observable']})
    return result
