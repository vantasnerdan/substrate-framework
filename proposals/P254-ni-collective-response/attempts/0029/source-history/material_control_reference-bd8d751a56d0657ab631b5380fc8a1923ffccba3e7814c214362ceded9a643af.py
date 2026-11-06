"""Exact per-case saved-density backgrounds; reuse the admitted A25 codec.

The caller makes A24/A25 serialization helpers importable. No loaded-material
production authority is a generic input gate. A loaded archive may have other
path labels only when every original input buffer has identical bytes/hash.
"""
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from .qe_saved_native import read_bound_bytes
from .qe_material_reference import load_material_inputs, reconstruct_reference


def bound_json(binding):
    raw, _ = read_bound_bytes(binding, max_bytes=32 << 20)
    return json.loads(raw)


def bound_arrays(binding):
    raw, _ = read_bound_bytes(binding, max_bytes=1 << 30)
    with np.load(io.BytesIO(raw), allow_pickle=False) as saved:
        return {name: saved[name] for name in saved.files}


def _content(source):
    return int(source['bytes']), source['sha256']


def restore_control_background(preparation_binding, receipt_binding):
    """Reusable matched-response intake; authenticate originals and raw arrays."""
    from material_reference_output import restore_background
    from material_pilot_io import array_record
    preparation = bound_json(preparation_binding)
    if preparation['status'] != 'completed':
        raise ValueError('completed original control preparation required')
    provenance = preparation['provenance']
    inputs = load_material_inputs(provenance['xml'], provenance['density'], provenance['upfs'])
    receipt = bound_json(receipt_binding)
    expected = [_content(s) for s in inputs.sources]
    if [_content(s) for s in receipt['metadata']['sources']] != expected:
        raise ValueError('background is not this exact XML/density/UPF content')
    arrays = bound_arrays(receipt['archive'])
    if set(arrays) != set(receipt['arrays']):
        raise ValueError('background array membership mismatch')
    for name, array in arrays.items():
        if array_record(array) != receipt['arrays'][name]:
            raise ValueError('background raw array digest mismatch: ' + name)
    background = restore_background(inputs, arrays)
    if [_content(s) for s in background.sources[:len(inputs.sources)]] != expected:
        raise ValueError('restored background sources differ from original inputs')
    return background, preparation, receipt


def retain_control_background(preparation_binding, library_binding, array_root, *,
                              loaded_receipt_binding=None):
    """Reconstruct only missing controls; loaded is always verbatim archive reuse.

    Return a receipt binding usable by restore_control_background. A completed
    nonloaded cache is recovered, not regenerated. Partial archives fail closed.
    """
    from material_pilot_io import binding, canonical, commit_block, recover_block
    from material_reference_output import background_arrays, raw_reference_metadata
    from .qe_internal_pbe import NativePBE
    preparation = bound_json(preparation_binding)
    if preparation['status'] != 'completed':
        raise ValueError('completed original control preparation required')
    p = preparation['provenance']
    if p['case'] == 'loaded-ni31h6':
        if loaded_receipt_binding is None:
            raise ValueError('loaded case requires the existing A25 receipt; reconstruction forbidden')
        # The restore path authenticates original byte equality, not path aliases.
        return dict(loaded_receipt_binding), 'existing-loaded-archive-verbatim'
    if loaded_receipt_binding is not None:
        raise ValueError('a nonloaded case cannot consume the loaded archive')
    inputs = load_material_inputs(p['xml'], p['density'], p['upfs'])
    metadata = dict(schema='material-control-reference-background-v1', case=p['case'],
        preparation=preparation_binding, sources=list(inputs.sources), library_json=library_binding,
        raw_reference_metadata=raw_reference_metadata(inputs.reference),
        physical_grid=list(inputs.reference.fft_grid),
        field_order='F physical grid; explicit up/down except rho_total_magnetization',
        units='density bohr^-3; potential Ry; gradient physical bohr derivatives')
    identity = hashlib.sha256(canonical(metadata)).hexdigest()
    target = Path(array_root) / p['case']
    target.mkdir(parents=True, exist_ok=True)
    key = 'background-coordinates'
    receipt = recover_block(target, key, identity, metadata)
    disposition = 'completed-control-cache-reuse'
    if receipt is None:
        backend = NativePBE(bound_json(library_binding))
        background = reconstruct_reference(inputs, backend)
        commit_block(target, key, background_arrays(background), metadata, identity, 1 << 30)
        disposition = 'new-original-control-background'
    return binding(target / (key + '.json')), disposition


def selected_control_hamiltonian(background, preparation):
    """Construct actual full-PW H0; authenticate the full selected original frame/D.

    No all-column H0 replay, new equilibrium, or state repair is performed.
    """
    from .material_control_inputs import load_control_inputs
    from .qe_saved_native import load_wavefunctions
    from .qe_real_bloch_grid import transport_real_frame, build_real_pw_factors
    from .qe_material_hamiltonian import MaterialHamiltonian
    p, selected = preparation['provenance'], preparation['selected']
    control = load_control_inputs(p['case_manifest'], p['native_receipt'], p['primary_admission'])
    for name in ('xml', 'density', 'upfs'):
        if control.provenance()[name] != p[name]:
            raise ValueError('preparation differs from original admitted ' + name)
    slot = next(s for s in control.grid.slots if s.ispin == selected['ispin']
                and s.mesh_index == tuple(selected['mesh_index']))
    wave = load_wavefunctions(slot.source_binding)
    frame = transport_real_frame(control.grid, slot, wave)
    arrays = bound_arrays(preparation['selected_arrays'])
    if (frame.basis.identity != selected['basis_identity']
            or frame.gauge_identity != selected['frame_identity']
            or frame.coefficients.shape != (selected['rows'], selected['columns'])):
        raise ValueError('selected full-frame identity mismatch')
    for name, expected in (('coefficients', frame.coefficients), ('miller', frame.basis.miller),
            ('reciprocal_bohr_inv', frame.basis.reciprocal_bohr_inv),
            ('k_cart_bohr_inv', frame.basis.k_cart_bohr_inv),
            ('cell_bohr', background.inputs.reference.cell_bohr),
            ('positions_bohr', background.inputs.reference.positions_bohr)):
        if not np.array_equal(arrays[name], expected):
            raise ValueError('original selected frame mismatch: ' + name)
    factors = build_real_pw_factors(background.inputs, frame.basis)
    return MaterialHamiltonian(background, factors, slot.ispin, 1), arrays
