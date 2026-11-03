"""Authenticate original case, complete slot lineage and existing source states."""
import numpy as np
import matched_reference_io as io
from substrate_framework.material_control_reference import bound_json
from substrate_framework.material_matched_reference import slot_identity
from substrate_framework.streamed_material_state import authenticate


def registration():
    record = io.binding(io.HERE/'matched-joint-case-registration.json')
    data = bound_json(record)
    for source in data['sources']:
        authenticate(source)
    paths = [io.HERE/name for name in ('matched_joint_case.py', 'matched_joint_sources.py',
        'matched_joint_state.py', 'run_matched_joint_case.py', 'run_matched_joint_oracle.py')]
    # The response owner may improve execution without changing this constructor.
    # Its exact final bytes are bound by the executed oracle/model, not frozen by
    # the earlier inventory of unchanged numerical primitives.
    paths += [io.ROOT/'src/substrate_framework/material_joint_response.py']
    implementation = [record, *data['sources'], *[io.binding(path) for path in paths]]
    return record, data, implementation


def validate_aggregate(value, context, rows):
    if (value['status'] != 'complete' or value['complete_case'] is not True
            or value['context'] != context or value['slot_count'] != 128
            or len(rows) != 128 or value['positive_time_symmetry'] is not False
            or value['model'] != 'matched-initial-source-orbit-declared-geometric-v1'
            or value['source_coefficients_recomputed'] is not False
            or value['coefficient_defect_upper_bound'] is not None):
        raise ValueError('complete original matched case with unchanged approximation required')
    expected = [slot_identity(row[1]) for row in rows]
    if ([entry['slot'] for entry in value['coefficient_lineage']] != expected
            or [entry['receipt'] for entry in value['coefficient_lineage']] != value['source_receipts']):
        raise ValueError('complete ordered original128-slot coefficient lineage differs')
    if value['original_source_count'] != len({(s.ispin, s.source_ik) for _, s, _ in rows}):
        raise ValueError('original independent source coverage differs')


def validate_atoms(control, background, arrays):
    ref = control.reference
    expected = dict(cell_bohr=ref.cell_bohr, positions_bohr=ref.positions_bohr,
        atom_indices=np.asarray(ref.atom_indices), atom_names=np.asarray(ref.atom_names),
        saved_density_bohr3=background.rho_updown)
    masses = {s.name: s.mass_amu for s in ref.species}
    expected['atomic_masses_amu'] = np.asarray([masses[n] for n in ref.atom_names])
    expected['valence_charges'] = np.asarray([background.inputs.tables[n].upf.z_valence for n in ref.atom_names])
    for key, value in expected.items():
        if not np.array_equal(arrays[key], value):
            raise ValueError('original case material identity differs: '+key)
    mesh = tuple(ref.monkhorst_pack['nk'+str(i)] for i in (1, 2, 3))
    if mesh != (4, 4, 4):
        raise ValueError('original BvK4^3 normalization required')


def source_catalog(aggregate, context, rows, ledger):
    from matched_joint_case import slot_identity as disk_identity
    states, metadata = {}, {}
    for entry, (index, slot, _) in zip(aggregate['coefficient_lineage'], rows):
        receipt = bound_json(entry['receipt'])
        c, m = receipt['context'], receipt['metadata']
        if (receipt['status'] != 'complete' or any(c.get(k) != v for k, v in context.items())
                or m['slot'] != slot_identity(slot)
                or m['basis_identity'] != str(ledger['basis_identity'][index])
                or c['original_source'] != dict(slot.source_binding)
                or c['registration'] != aggregate['original_direct_registration']
                or c['transport_registration'] != aggregate['registration']):
            raise ValueError('slot lineage is not this original case/ledger/source')
        initial = m['initial_state']
        if (initial['gram_conjugated'] != slot.orbital_conjugation
                or initial['density_coefficients_conjugated'] != slot.orbital_conjugation):
            raise ValueError('initial Gram AND D conjugation policy differs')
        key = slot.ispin, slot.source_ik
        record = c['state_receipt']
        state = bound_json(record)
        sc = state['context']
        if (state['status'] != 'complete' or any(sc.get(k) != v for k, v in context.items())
                or sc['original_source'] != dict(slot.source_binding)
                or (sc['spin'], sc['source_ik']) != key
                or sc['registration'] != aggregate['original_direct_registration']):
            raise ValueError('original source state authority differs')
        if key in states and states[key] != record:
            raise ValueError('same original source has conflicting state receipts')
        states[key] = record
        metadata[disk_identity(slot)] = m
    return states, metadata
