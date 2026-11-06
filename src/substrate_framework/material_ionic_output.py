"""Authenticated completed actual ionic baseline; no Ewald rerun on readback.

Uses existing A24 immutable archive/readback and A25 source/recipe authority.
The separately authenticated post-exit manager outcome remains Root-owned;
this reader authenticates the terminal's pre-publication snapshot only.
"""
from dataclasses import dataclass
import hashlib
from pathlib import Path
from types import MappingProxyType

from material_run_contract import (HERE, binding, bound_json, canonical,
                                   implementation, recover_block, verify_consumers)
from material_reference_output import load_arrays
from .material_ionic_recipe import (KEY, KIND, SCHEMA, block_metadata, load_ionic_inputs,
                                    preload_parsers, validate_recipe)
from .material_ionic_arrays import restore_term, validate_arrays

RECEIPT_FIELDS = ('recipe_sha256', 'archive', 'arrays', 'metadata', 'readback')


@dataclass(frozen=True)
class SavedIonicLinearTerm:
    term: object
    arrays: object
    comparison: object
    receipt: object
    recipe: object
    recipe_sha256: str


def authenticate_sources(recipe, producer_final_sources=()):
    """Bind original compound runtime authority and actual current consumers."""
    registration = bound_json(recipe['registration'])
    if (Path(recipe['registration']['path']).resolve() != HERE/'implementation-registration.json'
            or recipe['source_access'] != registration['source_access']):
        raise ValueError('wrong compound A25 runtime registration/source authority')
    bound_json(recipe['source_access'])
    declared = recipe['implementation_sources']
    if len({item['path'] for item in declared}) != len(declared):
        raise ValueError('duplicate implementation source')
    verify_consumers(recipe, producer_final_sources)
    records = {item['path']: item for item in declared}
    records.update({item['path']: item for item in producer_final_sources})
    for item in records.values():
        if binding(item['path']) != item:
            raise ValueError('original implementation source changed: '+item['path'])
    if records.get(str(Path(__file__).resolve())) != binding(__file__):
        raise ValueError('current saved ionic reader is not explicitly authenticated')
    actual = implementation()
    verify_consumers(recipe, actual)
    return actual


def load_saved_ionic_block(target, identity, expected_receipt, recipe, *, data=None):
    """Recover ONE fully committed baseline against a bound expected receipt.

    Used by producer recovery and the terminal reader. Caller authenticates
    recipe and receipt bytes first. No incomplete-array/sum continuation exists.
    Optional data reuses already admitted decoded original inputs in a producer.
    """
    target = Path(target).resolve()
    expected = {name: expected_receipt[name] for name in RECEIPT_FIELDS}
    if expected['recipe_sha256'] != identity or hashlib.sha256(canonical(recipe)).hexdigest() != identity:
        raise ValueError('saved ionic recipe identity differs')
    actual = recover_block(target, KEY, identity, expected['metadata'])
    if actual != expected or actual['archive']['path'] != str(target/(KEY+'.npz')):
        raise ValueError('completed ionic archive differs from authenticated receipt')
    if actual['readback'] != 'exact dtype/shape/value; allow_pickle=False':
        raise ValueError('unsupported original readback contract')
    data = load_ionic_inputs(recipe['inputs'], recipe['settings']) if data is None else data
    base = block_metadata(recipe, data)
    metadata = actual['metadata']
    if set(metadata) != set(base) | {'comparison'} or any(metadata[name] != value for name, value in base.items()):
        raise ValueError('saved ionic physical/source/Native-stage metadata differs')
    arrays = load_arrays(actual)
    comparison = validate_arrays(arrays, data)
    if comparison != metadata['comparison'] or metadata['schema'] != SCHEMA:
        raise ValueError('saved comparison does not match all original arrays')
    term = restore_term(arrays, data, recipe['settings'])
    return SavedIonicLinearTerm(term, MappingProxyType(arrays), MappingProxyType(comparison),
                                MappingProxyType(actual), MappingProxyType(recipe), identity)


def load_saved_ionic_linear_term(terminal_binding):
    """Authenticate original terminal, recipe, source, input and whole archive.

    Returns SavedIonicLinearTerm with all111 ionic-force/gradient coordinates,
    energy/components, original geometry/valences, Native decimal identities
    and raw differences. Decodes small XML/UPF/stdout again to verify identity;
    never recomputes Ewald, density, WFCs, PBE or a prepared state. Actual reader
    consumers are inventoried before/after readback against the bound universe.
    """
    preload_parsers()
    terminal = bound_json(terminal_binding)
    recipe, identity = terminal['recipe'], terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or recipe['kind'] != KIND or hashlib.sha256(canonical(recipe)).hexdigest() != identity):
        raise ValueError('complete authenticated actual ionic terminal required')
    validate_recipe(recipe)
    authenticate_sources(recipe, terminal['final_implementation_sources'])
    if len(terminal['outputs']) != 1 or terminal['outputs'][0]['key'] != KEY:
        raise ValueError('exactly one complete37-ion baseline required')
    entry = terminal['outputs'][0]
    directory = Path(terminal_binding['path']).resolve().parent
    result = load_saved_ionic_block(directory, identity, entry, recipe)
    authenticate_sources(recipe, terminal['final_implementation_sources'])
    return result
