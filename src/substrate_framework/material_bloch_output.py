"""Authenticated complete20-source/128-slot lazy Bloch material coverage.

Attempt0024/0025 control modules must be on PYTHONPATH. This reader consumes
completed archives and reconstructs metadata, never H0, ionic actions or WFC
coefficients. Post-exit manager success is separate mandatory Root authority.
"""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from types import MappingProxyType
import xml.etree.ElementTree as ET
import numpy as np
from material_run_contract import binding,bound_json,canonical,recover_block,verify_consumers
from material_reference_output import load_arrays
from .material_bloch_admission import (RECEIPT_FIELDS,validate_recipe,wave_bindings,
                                       source_metadata,final_metadata,background_inputs)
from .material_bloch_records import BlochSlotRecord,decode_records,validate_source_arrays
from .material_bloch_probe import validate_probe_arrays
from .material_dyad_output import REFERENCE_SHA256
from .qe_saved_reference import load_saved_reference
from .qe_bloch_grid import build_bloch_grid


@dataclass(frozen=True)
class SavedBlochCoverage:
    grid: object
    slot_records: tuple[BlochSlotRecord,...]
    source_receipts: object
    coverage_identity: str
    recipe: object
    probe_receipt: object
    probe_comparisons: object


def preload_parsers():
    """Load XML dispatch before inventory; inherited IO preloads NPZ dispatch."""
    ET.fromstring('<source-inventory-preload/>')


def _receipt(entry, identity, metadata):
    expected = {name:entry[name] for name in RECEIPT_FIELDS}
    path = Path(expected['archive']['path'])
    if expected['metadata'] != metadata or path.stem != entry['key']:
        raise ValueError('archive metadata/key differs from complete recipe')
    actual = recover_block(path.parent,path.stem,identity,metadata)
    if actual != expected:
        raise ValueError('archive/receipt exact identity differs')
    return actual


def _authenticate_sources(terminal):
    recipe = terminal['recipe']
    sources = recipe['implementation_sources']
    if len({item['path'] for item in sources}) != len(sources):
        raise ValueError('duplicate producer source binding')
    final = terminal['final_implementation_sources']
    verify_consumers(recipe,final)
    all_sources = {item['path']:item for item in sources+final}
    for item in all_sources.values():
        if binding(item['path']) != item:
            raise ValueError('bound producer/consumer source changed: '+item['path'])
    if all_sources.get(str(Path(__file__).resolve())) != binding(__file__):
        raise ValueError('saved coverage reader lacks source authority')


def _background_receipt(inputs):
    record = inputs['reference_terminal']
    if record['sha256'] != REFERENCE_SHA256 or str(record['bytes']) != '1000830':
        raise ValueError('accepted original A25 background authority required')
    terminal = bound_json(record)
    recipe, identity = terminal['recipe'], terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or hashlib.sha256(canonical(recipe)).hexdigest() != identity):
        raise ValueError('accepted background recipe/completion differs')
    reduced = background_inputs(inputs)
    if any(reduced[key] != recipe['inputs'][key] for key in ('xml','density','upfs','wavefunctions')):
        raise ValueError('actual inputs do not share accepted original background')
    entries = [item for item in terminal['outputs'] if item['key'] == 'background-coordinates']
    if len(entries) != 1:
        raise ValueError('one complete accepted background archive required')
    receipt = {key:entries[0][key] for key in RECEIPT_FIELDS}
    if binding(receipt['archive']['path']) != receipt['archive']:
        raise ValueError('accepted background bytes changed')
    return receipt


def load_saved_bloch_coverage(terminal_binding):
    """Authenticate complete archives, actual source bytes and derived slot metadata."""
    terminal = bound_json(terminal_binding)
    recipe, identity = terminal['recipe'],terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or recipe['kind'] != 'material-bloch-grid'
            or hashlib.sha256(canonical(recipe)).hexdigest() != identity):
        raise ValueError('complete immutable Bloch terminal/recipe required')
    validate_recipe(recipe)
    _authenticate_sources(terminal)
    inputs = recipe['inputs']
    for record in (inputs['xml'],inputs['density'],*inputs['upfs'].values(),*wave_bindings(inputs).values()):
        if binding(record['path']) != record:
            raise ValueError('actual original input bytes changed')
    ref = load_saved_reference(inputs['xml'])
    grid = build_bloch_grid(ref,wave_bindings(inputs),symmetry_policy=recipe['settings']['symmetry_policy'])
    if len(grid.slots) != 128 or ref.nat != 37 or ref.nbnd_by_spin != (346,346):
        raise ValueError('complete actual37atom/128slot coverage required')
    entries = {entry['key']:entry for entry in terminal['outputs']}
    expected_keys = {f'source-{spin}-{ik}' for spin in (1,2) for ik in range(1,11)} | {'coverage-probe'}
    if len(entries) != len(terminal['outputs']) or set(entries) != expected_keys:
        raise ValueError('all20 unique source boundaries and final probe boundary required')
    records, receipts, receipt_bindings = [], {}, {}
    for spin in (1,2):
        for ik in range(1,11):
            key = f'source-{spin}-{ik}'
            receipt = _receipt(entries[key],identity,source_metadata(recipe,spin,ik))
            arrays = load_arrays(receipt)
            slots = [slot for slot in grid.slots if (slot.ispin,slot.source_ik)==(spin,ik)]
            records.extend(validate_source_arrays(ref,slots,arrays))
            receipts[spin,ik] = receipt
            receipt_bindings[key] = binding(Path(receipt['archive']['path']).with_suffix('.json'))
            del arrays
    records.sort(key=lambda r:(r.mesh_slot,r.ispin))
    if len(records) != 128 or len({(r.ispin,r.mesh_index) for r in records}) != 128:
        raise ValueError('complete unique typed spin/grid descriptors required')
    background = _background_receipt(inputs)
    final = _receipt(entries['coverage-probe'],identity,final_metadata(recipe,background))
    arrays = load_arrays(final)
    if decode_records(arrays.pop('slot_records_json')) != tuple(records):
        raise ValueError('final coverage list differs from complete source boundaries')
    if json.loads(str(arrays.pop('source_receipt_bindings_json'))) != receipt_bindings:
        raise ValueError('final coverage does not bind every original completed source receipt')
    comparisons = validate_probe_arrays(arrays,grid,records)
    return SavedBlochCoverage(grid,tuple(records),MappingProxyType(receipts),identity,
                              MappingProxyType(recipe),final,MappingProxyType(comparisons))
