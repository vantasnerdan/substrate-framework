#!/usr/bin/env python3
"""Aggregate exact original slots with declared initial-coefficient reuse.

Uses the unchanged matched numerical reducer and frozen ionic assembly. The
original DIRECT runner remains unchanged and intentionally rejects transported
receipts. This bridge never labels geometric predictions direct/exact.
"""
import argparse
import json
from pathlib import Path

from matched_reference_io import (binding, save_json, save_arrays, runtime, registration, load_case)
from matched_transport_io import transport_registration, slot_catalog, authenticated_arrays


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', required=True)
    p.add_argument('--source-dir', action='append', type=Path, required=True)
    p.add_argument('--array-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--accept-declared-geometry', action='store_true')
    p.add_argument('--max-array-bytes', type=int, default=1 << 30)
    args = p.parse_args()
    if not args.accept_declared_geometry:
        p.error('--accept-declared-geometry required; no inherited material error budget')
    runtime(args.array_root)
    import numpy as np
    from substrate_framework.material_matched_reference import (aggregate_complete,
        assemble_frozen, finite_ionic_reference, slot_identity, KEYS)
    from substrate_framework.material_matched_transport import (GEOMETRY_POLICY, identity_representative,
                                                               require_arrays)
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.material_force_lattice import lattice_residuals
    reg_binding, reg = registration()
    transport_binding, _ = transport_registration()
    control, background, ledger, rows, context = load_case(args.case, reg)
    catalog = slot_catalog(args.source_dir, context, reg_binding, transport_registration=transport_binding)
    expected = {(s.ispin, s.mesh_index): (i, s, pop) for i, s, pop in rows}
    if set(catalog) != set(expected):
        raise ValueError('complete nonduplicated original target set required: '+str(dict(
            missing=sorted(set(expected)-set(catalog)), extra=sorted(set(catalog)-set(expected)))))
    receipts, lineage, comparisons = [], [], []
    def records():
        for key in sorted(expected):
            record, value = catalog[key]
            row = expected[key]
            slot, meta, c = row[1], value['metadata'], value['context']
            arrays = authenticated_arrays(catalog[key], row, ledger)
            require_arrays(arrays, control.reference.nat, int(np.prod(control.reference.fft_smooth)))
            transformed = meta['geometry_policy'] == GEOMETRY_POLICY
            if transformed:
                representative = identity_representative(control.grid, slot.ispin, slot.source_ik)
                source = bound_json(c['representative_receipt'])
                if (meta['original_representative'] != slot_identity(representative)
                        or source['metadata']['slot'] != slot_identity(representative)
                        or source['metadata']['geometry_policy'] != 'direct-original-target'
                        or any(source['context'].get(k) != v for k, v in context.items())
                        or source['context']['registration'] != reg_binding
                        or meta['initial_state']['gram_conjugated'] != slot.orbital_conjugation
                        or meta['initial_state']['density_coefficients_conjugated'] != slot.orbital_conjugation):
                    raise ValueError('original direct representative or initial orbital lineage differs')
            lineage.append(dict(slot=slot_identity(slot), receipt=record,
                coefficient_origin='transformed-original-representative' if transformed else 'direct-original-target',
                representative_receipt=c.get('representative_receipt'),
                geometry_policy=meta['geometry_policy']))
            if 'direct_comparison' in meta:
                comparisons.append(meta['direct_comparison'])
            receipts.append(record)
            yield meta['slot'], arrays
    aggregate = aggregate_complete(control.grid, records())
    ref = control.reference
    integrated = abs(np.linalg.det(ref.cell_bohr))/len(aggregate['density_bohr3'])*aggregate['density_bohr3'].sum(axis=0)
    expected_number = np.zeros(2)
    for _, slot, pop in rows:
        expected_number[slot.ispin-1] += slot.weight*pop.sum()
    charge_error = max(float(np.max(abs(integrated-expected_number))),
                       float(np.max(abs(aggregate['weighted_spin_number']-expected_number))))
    if charge_error > 1e-8*(1+expected_number.max()):
        raise ValueError('exact original spin/mesh weighted charge mismatch')
    linear, ionic_q, masses = finite_ionic_reference(background.inputs, reg['ionic_contract'])
    arrays = dict(aggregate)
    arrays.update(assemble_frozen(aggregate, linear.ionic_force, ionic_q))
    arrays.update(ionic_q_grid_ry_bohr2=ionic_q, atomic_masses_amu=masses,
        valence_charges=linear.valence_charges, cell_bohr=ref.cell_bohr,
        positions_bohr=ref.positions_bohr, atom_indices=np.asarray(ref.atom_indices),
        atom_names=np.asarray(ref.atom_names), integrated_spin_number=integrated,
        original_ledger_spin_number=expected_number, saved_density_bohr3=background.rho_updown,
        kernel_potential_ry=background.apply_kernel(aggregate['density_bohr3']-background.rho_updown))
    args.output.mkdir(parents=True, exist_ok=False)
    archive = save_arrays(args.array_root.resolve()/(args.case+'-matched-transport.npz'), arrays, args.max_array_bytes)
    transformed_count = sum(row['coefficient_origin']=='transformed-original-representative' for row in lineage)
    result = dict(status='complete', complete_case=True,
        model='matched-initial-source-orbit-declared-geometric-v1', direct_model=reg['model'],
        registration=transport_binding, original_direct_registration=reg_binding,
        context=context, archive=archive, source_receipts=receipts, coefficient_lineage=lineage,
        slot_count=len(receipts), original_source_count=len({(s.ispin, s.source_ik) for _, s, _ in rows}),
        direct_slot_count=len(receipts)-transformed_count, transformed_slot_count=transformed_count,
        geometry_policy=GEOMETRY_POLICY, direct_target_comparisons=comparisons,
        coefficient_defect_upper_bound=None,
        unbounded_coefficient_defect_scope=[k for k in KEYS if k not in ('density_bohr3', 'number')],
        approximation='original rounded geometry unchanged; no full-case coefficient error budget assigned',
        exact_scope='complete original slots/weights/populations; kinematic density and number transport',
        diagnostics=dict(number_residual_max=charge_error,
            force_imaginary_max=float(np.max(abs(arrays['force_at_origin_ry_bohr'].imag))),
            onsite_hermiticity_max=float(np.max(abs(arrays['electronic_onsite_matrix_ry_bohr2']-
                                                   arrays['electronic_onsite_matrix_ry_bohr2'].conj().T))),
            field_contact_hermiticity_max=float(np.max(abs(arrays['contact']-arrays['contact'].conj().T))),
            force_norm=float(np.linalg.norm(arrays['force_at_origin_ry_bohr'])),
            lattice=lattice_residuals(arrays['common_force_lattice_ry_bohr2'])),
        ionic_reference=dict(authority=linear.ionic_authority,
            quadratic='same original full signed-density-G finite ionic construction as direct matched aggregator',
            finite_sum_difference=reg['ionic_finite_sum_qualification']),
        scientific_limits=reg['scientific_limits'], positive_time_symmetry=False,
        response_comparison=False, source_coefficients_recomputed=False,
        direct_method='unchanged run_matched_reference_source.py and run_matched_reference_aggregate.py')
    print(json.dumps(dict(result=save_json(args.output/'result.json', result)), sort_keys=True))


if __name__ == '__main__':
    main()
