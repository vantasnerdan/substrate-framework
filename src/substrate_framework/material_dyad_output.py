"""Authenticated completed Gamma dyad maps; no map recomputation on readback.

The attempt0024/0025 control modules must be on PYTHONPATH (the pilot adds
these paths). This module reuses their archive, restore and publication
protocols, rather than defining a second runtime. Sources authenticate the
current consumer and the producer's bound recipe; historical background
producer bytes are provenance, not permission to rerun that producer.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from types import MappingProxyType
import xml.etree.ElementTree as ET

import numpy as np
from material_run_contract import binding, bound_json, canonical, recover_block
from material_reference_output import load_arrays, restore_background
from .qe_material_reference import load_material_inputs
from .qe_saved_native import load_wavefunctions  # preload the real cold binary reader

SCHEMA = 'material-dyad-maps-v1'
SCOPE = 'selected k=q=Gamma; both spins separately; all original PW rows and 111 coordinates; not complete BvK'
STATE_SCOPE = 'Raw non-Hermitian operator L C Rdagger; C is NOT state D or a physical preparation'
REFERENCE_SHA256 = '52fa2c8a2ed5ad7e216b4f3e041ab80dccfd8daf119c1f5cd2b247a381ffcccb'
CONTINUATION_SHA256 = 'd7af1004527f6aa3d36e0b0db69c0c87a2bf12e587a95dc8f973ef8a76e31450'
RECEIPT_FIELDS = ('recipe_sha256', 'archive', 'arrays', 'metadata', 'readback')
COMPARISONS = ('integrated_trace', 'density_potential_bilinear', 'gradient_ionic_bilinear')
UNITS = {'density': 'bohr^-3', 'positive_gradient': 'Ry/bohr', 'potential': 'Ry',
         'potential_action': 'Ry', 'ionic_action': 'Ry/bohr', 'direction': 'dimensionless',
         'left': 'dimensionless PW coefficients', 'right': 'dimensionless PW coefficients',
         'coefficient': 'dimensionless raw operator coefficient; not state D',
         'comparison_pairs': ['dimensionless', 'Ry', 'Ry/bohr']}


def preload_parsers():
    """Load actual XML and NPZ dispatch before inventory, without material work.

    material_pilot_io already exercises np.savez/np.load with an empty array.
    Expat is lazy in ElementTree: instantiate its actual parser here too.
    All material/output modules above are loaded, not just whatever happened
    to be resident when a caller asked for implementation().
    """
    ET.fromstring('<source-inventory-preload/>')


def validate_grant(grant):
    inputs, settings = grant['inputs'], grant['settings']
    if set(inputs) != {'xml', 'density', 'upfs', 'wavefunctions', 'reference_terminal'}:
        raise ValueError('exact XML/density/UPF/both-spin/reference-terminal input bindings required')
    if set(inputs['upfs']) != {'Ni', 'H'} or set(inputs['wavefunctions']) != {'1', '2'}:
        raise ValueError('both original species and spins required')
    if set(settings) != {'block_size', 'probe_recipe', 'continuation_registration'}:
        raise ValueError('explicit block_size, probe_recipe and continuation_registration required')
    if type(settings['block_size']) is not int or not 1 <= settings['block_size'] <= 16:
        raise ValueError('block_size must be an integer1..16')
    if settings['probe_recipe'] != 'saved-bands1-2-complex-all-coordinate-v1':
        raise ValueError('only the registered raw operator probe recipe is supported')
    continuation = bound_json(settings['continuation_registration'])
    if (settings['continuation_registration']['sha256'] != CONTINUATION_SHA256
            or 'A_actual_material_consumer' not in continuation):
        raise ValueError('exact0026 continuation authority required')
    if grant.get('library_binding') is not None:
        raise ValueError('no PBE backend is loaded or admitted by this consumer')
    receipts = grant['resume_receipts']
    if (not isinstance(receipts, dict) or not set(receipts) <= {'spin-1-dyad', 'spin-2-dyad'}
            or (not grant['resume'] and receipts)):
        raise ValueError('explicit completed-block receipt bindings required; fresh uses {}')


def load_accepted_background(inputs):
    """Restore accepted0025 verbatim; radial tables/FFT adapter are charged setup.

    No backend, reconstruct_reference, SCF or PH producer is invoked. The
    original immutable terminal authenticates the complete background receipt.
    """
    source = inputs['reference_terminal']
    if source['sha256'] != REFERENCE_SHA256 or str(source['bytes']) != '1000830':
        raise ValueError('requires accepted original0025 material-reference terminal')
    terminal = bound_json(source)
    recipe, identity = terminal['recipe'], terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or hashlib.sha256(canonical(recipe)).hexdigest() != identity):
        raise ValueError('accepted reference recipe or completion differs')
    for key in ('xml', 'density', 'upfs', 'wavefunctions'):
        if inputs[key] != recipe['inputs'][key]:
            raise ValueError('input differs from accepted original0025: '+key)
    entries = [entry for entry in terminal['outputs'] if entry['key'] == 'background-coordinates']
    if len(entries) != 1:
        raise ValueError('exactly one accepted background required')
    expected = {name: entries[0][name] for name in RECEIPT_FIELDS}
    archive = Path(expected['archive']['path'])
    actual = recover_block(archive.parent, archive.stem, identity, expected['metadata'])
    if actual != expected:
        raise ValueError('accepted background archive/receipt mismatch')
    material = load_material_inputs(inputs['xml'], inputs['density'], inputs['upfs'])
    ref = material.reference
    if (ref.nat != 37 or ref.atom_names != ('Ni',)*31+('H',)*6
            or ref.atom_indices != tuple(range(1, 38)) or ref.nbnd_by_spin != (346, 346)
            or ref.fft_grid != (72, 72, 72)):
        raise ValueError('original37atom/111coordinate/72cubed material required')
    background = restore_background(material, load_arrays(actual))
    return background, expected


def block_metadata(inputs, settings, spin, background_receipt):
    return dict(schema=SCHEMA, spin_one_based=spin, k_spin_local_one_based=1,
                inputs=inputs, settings=settings, background_receipt=background_receipt,
                full_pw_dimension=18613, saved_band_count=346, dyad_bands_one_based=[1, 2],
                physical_grid=[72, 72, 72], field_order='F physical grid',
                coordinate_order='atom-major: original atoms1..37, Cartesian x,y,z',
                positive_gradient_sign='positive dH/du, NOT force', units=UNITS,
                scope=SCOPE, state_scope=STATE_SCOPE, D=None, preparation=None,
                physical_temperature=None, held_sector=None, force_offset=None,
                comparison_order=list(COMPARISONS), rtol=1e-10, atol=1e-12,
                allowance_scope='new transpose/integrated identities only, not material convergence')


def comparison_report(pairs):
    """Report raw complex scalars; neither conjugate bilinear fields nor repair."""
    result = {}
    for name, (left, right) in zip(COMPARISONS, pairs):
        error, allowance = float(abs(left-right)), float(1e-12+1e-10*abs(right))
        result[name] = {'left': [float(left.real), float(left.imag)],
                        'right': [float(right.real), float(right.imag)],
                        'absolute_error': error, 'allowance': allowance,
                        'passed': bool(error <= allowance)}
    return result


def validate_arrays(arrays):
    shapes = {'left_pw': (18613, 2), 'right_pw': (18613, 2), 'coefficient': (2, 2),
              'density_bohr_minus3': (72**3,), 'positive_gradient_ry_per_bohr': (111,),
              'potential_ry': (72**3,), 'potential_action_ry': (18613, 2),
              'direction_dimensionless': (111,), 'ionic_action_ry_per_bohr': (18613, 2),
              'comparison_pairs': (3, 2), 'pw_miller': (3, 18613),
              'positions_bohr': (3, 37), 'cell_bohr': (3, 3),
              'reciprocal_bohr_inv': (3, 3), 'atom_one_based': (37,),
              'atom_symbols': (37,), 'fft_shape': (3,), 'dyad_bands_one_based': (2,)}
    if set(arrays) != set(shapes):
        raise ValueError('saved map array layout differs')
    for name, shape in shapes.items():
        value = arrays[name]
        if value.shape != shape or (value.dtype.kind in 'fc' and not np.isfinite(value).all()):
            raise ValueError('invalid saved map array: '+name)
        if name in tuple(shapes)[:10] and value.dtype != np.dtype('complex128'):
            raise ValueError('raw complex128 maps required: '+name)
    if (not np.array_equal(arrays['atom_one_based'], np.arange(1, 38))
            or not np.array_equal(arrays['atom_symbols'], ['Ni']*31+['H']*6)
            or not np.array_equal(arrays['fft_shape'], [72]*3)
            or not np.array_equal(arrays['dyad_bands_one_based'], [1, 2])):
        raise ValueError('saved map coordinate identities differ')
    reports = comparison_report(arrays['comparison_pairs'])
    return reports


def load_saved_block(target, key, identity, expected_receipt):
    """Caller owns authenticated recipe/receipt; reads arrays, never invokes maps."""
    target = Path(target).resolve()
    if key not in ('spin-1-dyad', 'spin-2-dyad'):
        raise ValueError('unknown completed dyad block')
    expected = {name: expected_receipt[name] for name in RECEIPT_FIELDS}
    actual = recover_block(target, key, identity, expected['metadata'])
    if actual != expected or actual['archive']['path'] != str(target/(key+'.npz')):
        raise ValueError('saved dyad block does not match authenticated receipt')
    if actual['metadata']['schema'] != SCHEMA:
        raise ValueError('unsupported dyad output schema')
    arrays = load_arrays(actual)
    reports = validate_arrays(arrays)
    if reports != actual['metadata']['comparisons']:
        raise ValueError('saved raw comparison scalars differ from receipt')
    for array in arrays.values():
        array.flags.writeable = False
    return MappingProxyType(arrays), reports


def load_saved_maps(terminal_binding):
    """Authenticate complete source/terminal/archive bytes without recomputing.

    Returns {spin: {arrays, comparisons, receipt}}. Manager post-exit success
    remains separate owner authority, as in the existing publication protocol.
    This API does not claim crash recovery or reread/rerun historical producers.
    """
    from material_run_contract import verify_consumers
    terminal = bound_json(terminal_binding)
    recipe, identity = terminal['recipe'], terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or recipe['kind'] != 'material-dyad-maps'
            or hashlib.sha256(canonical(recipe)).hexdigest() != identity):
        raise ValueError('complete authenticated dyad recipe required')
    sources = recipe['implementation_sources']
    if len({item['path'] for item in sources}) != len(sources):
        raise ValueError('duplicate source binding')
    final_sources = terminal['final_implementation_sources']
    verify_consumers(recipe, final_sources)
    all_sources = {item['path']: item for item in sources}
    all_sources.update({item['path']: item for item in final_sources})
    for source in all_sources.values():
        if binding(source['path']) != source:
            raise ValueError('bound implementation source changed: '+source['path'])
    if all_sources.get(str(Path(__file__).resolve())) != binding(__file__):
        raise ValueError('current saved-map consumer is not source-authenticated')
    result = {}
    if len(terminal['outputs']) != 2:
        raise ValueError('both complete spins required')
    for entry in terminal['outputs']:
        spin = entry['metadata']['spin_one_based']
        if spin not in (1, 2) or spin in result or entry['key'] != f'spin-{spin}-dyad':
            raise ValueError('spin identity/coverage differs')
        metadata = entry['metadata']
        base = block_metadata(recipe['inputs'], recipe['settings'], spin,
                              terminal['accepted_background'])
        if any(metadata.get(name) != value for name, value in base.items()):
            raise ValueError('saved block inputs/units/scope differ from authenticated recipe')
        if not metadata['raw_operator_diagnostic']['nonhermitian_witness']:
            raise ValueError('raw operator must remain explicitly non-Hermitian')
        directory = Path(entry['archive']['path']).parent
        arrays, reports = load_saved_block(directory, entry['key'], identity, entry)
        if not all(item['passed'] for item in reports.values()):
            raise ValueError('terminal reports completion with a failed map comparison')
        result[spin] = {'arrays': arrays, 'comparisons': reports, 'receipt': entry}
    return result
