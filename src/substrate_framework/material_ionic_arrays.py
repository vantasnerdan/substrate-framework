"""Typed whole-baseline ionic arrays and rounded-force comparison diagnostics.

No Ewald evaluation occurs here. The comparison is an observation, not an
accuracy gate: Native and the finite scalar may have different cutoffs, and
F14.8 values are rounded. Neither the candidate nor the comparator is centered.
"""
import json
from types import MappingProxyType

import numpy as np
from .qe_ionic_ewald import IonicLinearTerm
from .material_ionic_recipe import COORDINATES, sum_contract

FLOAT_SHAPES = {
    'cell_bohr': (3, 3), 'positions_bohr': (37, 3), 'valence_charges': (37,),
    'energy_ry': (), 'energy_terms_ry': (4,), 'force_ry_bohr': (37, 3),
    'ionic_force_ry_bohr': (111,), 'gradient_ry_bohr': (111,),
    'real_force_ry_bohr': (37, 3), 'reciprocal_force_ry_bohr': (37, 3),
    'native_force_ry_bohr': (37, 3), 'difference_ry_bohr': (37, 3),
}
JSON_FIELDS = ('sum_metadata_json', 'source_authority_json', 'ionic_authority_json')


def baseline_arrays(term, data):
    native = data.native_ledger.blocks['ion_ion']
    arrays = {name: np.asarray(getattr(term, name), dtype=np.float64) for name in (
        'cell_bohr', 'positions_bohr', 'valence_charges', 'energy_ry', 'force_ry_bohr',
        'gradient_ry_bohr', 'real_force_ry_bohr', 'reciprocal_force_ry_bohr')}
    arrays.update(energy_terms_ry=np.asarray([term.energy_terms_ry[k] for k in
                                             ('real', 'reciprocal', 'self', 'background')]),
                  ionic_force_ry_bohr=term.ionic_force,
                  native_force_ry_bohr=native.values_ry_per_bohr,
                  difference_ry_bohr=term.force_ry_bohr-native.values_ry_per_bohr,
                  atom_indices=np.asarray(data.reference.atom_indices, dtype=np.int64),
                  native_row_lines=np.asarray(native.row_lines, dtype=np.int64),
                  native_decimal_values=np.asarray(native.decimal_values),
                  sum_metadata_json=np.asarray(json.dumps(dict(term.sum_metadata), sort_keys=True, allow_nan=False)),
                  source_authority_json=np.asarray(json.dumps(term.source_authority)),
                  ionic_authority_json=np.asarray(term.ionic_authority))
    return arrays


def comparison_report(arrays):
    """Retain full signed residual separately; F14.8 half-step is not a tolerance."""
    difference = arrays['difference_ry_bohr']
    absolute = np.abs(difference)
    half_step = 5e-9
    return dict(max_abs_difference_ry_bohr=float(absolute.max()),
                rms_difference_ry_bohr=float(np.sqrt(np.mean(difference*difference))),
                per_atom_max_abs_difference_ry_bohr=absolute.max(axis=1).tolist(),
                native_decimal_half_step_ry_bohr=half_step,
                coordinates_outside_decimal_half_step=int(np.count_nonzero(absolute > half_step)),
                maximum_distance_beyond_decimal_half_step_ry_bohr=float(np.maximum(absolute-half_step, 0).max()),
                unprocessed_candidate_net_force_ry_bohr=arrays['force_ry_bohr'].sum(axis=0).tolist(),
                rounded_native_net_force_ry_bohr=arrays['native_force_ry_bohr'].sum(axis=0).tolist(),
                claim='observed finite-sum minus rounded original ion-ion force; no accuracy certificate or matching gate')


def validate_arrays(arrays, data):
    """Authenticate typed archived physical data against original decoded inputs."""
    expected = set(FLOAT_SHAPES) | set(JSON_FIELDS) | {'atom_indices', 'native_row_lines', 'native_decimal_values'}
    if set(arrays) != expected:
        raise ValueError('whole ionic archive array membership differs')
    for name, shape in FLOAT_SHAPES.items():
        a = arrays[name]
        if a.shape != shape or a.dtype != np.dtype('float64') or not np.isfinite(a).all():
            raise ValueError('invalid ionic float array: '+name)
    for name in ('atom_indices', 'native_row_lines'):
        if arrays[name].shape != (37,) or arrays[name].dtype != np.dtype('int64'):
            raise ValueError('invalid integer identity array: '+name)
    if arrays['native_decimal_values'].shape != (37, 3) or arrays['native_decimal_values'].dtype.kind != 'U':
        raise ValueError('all original Native decimal strings required')
    for name in JSON_FIELDS:
        if arrays[name].shape != () or arrays[name].dtype.kind != 'U':
            raise ValueError('scalar Unicode provenance field required: '+name)
    native = data.native_ledger.blocks['ion_ion']
    originals = dict(cell_bohr=data.reference.cell_bohr, positions_bohr=data.reference.positions_bohr.T,
                     valence_charges=data.valence_charges, atom_indices=data.reference.atom_indices,
                     native_force_ry_bohr=native.values_ry_per_bohr, native_row_lines=native.row_lines,
                     native_decimal_values=native.decimal_values)
    for name, expected_value in originals.items():
        if not np.array_equal(arrays[name], np.asarray(expected_value)):
            raise ValueError('archive differs from original geometry/UPF/Native source: '+name)
    equalities = (
        (arrays['force_ry_bohr'], arrays['real_force_ry_bohr']+arrays['reciprocal_force_ry_bohr']),
        (arrays['ionic_force_ry_bohr'], arrays['force_ry_bohr'].ravel()),
        (arrays['gradient_ry_bohr'], -arrays['ionic_force_ry_bohr']),
        (arrays['difference_ry_bohr'], arrays['force_ry_bohr']-arrays['native_force_ry_bohr']),
        (arrays['energy_ry'], sum(float(value) for value in arrays['energy_terms_ry'])),
    )
    if any(not np.array_equal(a, b) for a, b in equalities):
        raise ValueError('force sign/component, energy decomposition or comparison arrays disagree')
    return comparison_report(arrays)


def restore_term(arrays, data, settings):
    """Restore the committed IonicLinearTerm; never recompute its finite sum."""
    validate_arrays(arrays, data)
    metadata = json.loads(str(arrays['sum_metadata_json']))
    if (metadata['energy_units'] != 'Ry' or metadata['force_units'] != 'Ry/bohr'
            or metadata['coordinate_order'] != 'atom-major Cartesian xyz'
            or metadata['cell_vectors'] != 'columns' or metadata['e2_ry_bohr'] != 2.
            or metadata['processing'] != 'none: no net-force subtraction, symmetrization or fitted support'
            or not np.isfinite(metadata['real_cutoff_gap_bohr']) or metadata['real_cutoff_gap_bohr'] <= 0):
        raise ValueError('archived finite-sum conventions differ')
    for name in ('real_ordered_pair_terms', 'reciprocal_signed_vectors'):
        if type(metadata[name]) is not int or metadata[name] < 0:
            raise ValueError('invalid finite sum count')
    for value in arrays.values():
        value.flags.writeable = False
    term = IonicLinearTerm(
        float(arrays['energy_ry']), arrays['force_ry_bohr'], arrays['real_force_ry_bohr'],
        arrays['reciprocal_force_ry_bohr'], arrays['cell_bohr'], arrays['positions_bohr'],
        arrays['valence_charges'], COORDINATES, json.loads(str(arrays['source_authority_json'])),
        tuple(MappingProxyType(dict(source)) for source in data.sources), sum_contract(settings),
        MappingProxyType(dict(zip(('real', 'reciprocal', 'self', 'background'),
                                 (float(value) for value in arrays['energy_terms_ry'])))), MappingProxyType(metadata))
    if term.ionic_authority != str(arrays['ionic_authority_json']):
        raise ValueError('restored ionic authority differs from saved geometry/source/sum contract')
    return term
