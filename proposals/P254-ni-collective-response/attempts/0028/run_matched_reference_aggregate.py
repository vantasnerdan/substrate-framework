#!/usr/bin/env python3
"""Deterministically join every original target, then construct full ionic q lattice.

Supply each original-source output directory once. Incomplete or duplicate
coverage is rejected BEFORE ionic construction; no missing source gets weight.
"""
import argparse
from pathlib import Path
from matched_reference_io import (binding, save_json, save_arrays, runtime, registration,
                                  load_case, receipt_arrays)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--source-dir', action='append', type=Path, required=True)
    parser.add_argument('--array-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-array-bytes', type=int, default=1 << 30)
    args = parser.parse_args()
    runtime(args.array_root)
    import numpy as np
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.material_matched_reference import (aggregate_complete,
        assemble_frozen, finite_ionic_reference, slot_identity)
    from substrate_framework.material_force_lattice import lattice_residuals
    reg_binding, reg = registration()
    control, background, ledger, rows, context = load_case(args.case, reg)
    catalog = {}
    for directory in args.source_dir:
        for path in sorted(directory.glob('slot-*.json')):
            record = binding(path)
            value = bound_json(record)
            if value['status'] != 'complete':
                raise ValueError('incomplete target receipt')
            c = value['context']
            if any(c.get(k) != v for k, v in context.items()) or c['registration'] != reg_binding:
                raise ValueError('mixed case/model/implementation/preparation receipts')
            identity = value['metadata']['slot']
            key = (identity['spin'], tuple(identity['mesh_index']))
            if key in catalog:
                raise ValueError('duplicated target input')
            catalog[key] = (record, value)
    expected = {(s.ispin, tuple(s.mesh_index)): (i, s, pop) for i, s, pop in rows}
    if set(catalog) != set(expected):
        raise ValueError('complete nonduplicated original target set required: '+
                         str(dict(missing=sorted(set(expected)-set(catalog)),
                                  extra=sorted(set(catalog)-set(expected)))))
    receipts = []
    def records():
        for key in sorted(catalog):
            record, value = catalog[key]
            i, slot, pop = expected[key]
            meta = value['metadata']
            c = value['context']
            if (meta['slot'] != slot_identity(slot) or meta['basis_identity'] != str(ledger['basis_identity'][i])
                    or meta['geometry_policy'] != 'direct-original-target'
                    or meta['normalization'] != 'unweighted-one-explicit-spin-slot'
                    or c['original_source'] != dict(slot.source_binding)
                    or c['spin'] != slot.ispin or c['source_ik'] != slot.source_ik
                    or meta['rows'] != slot.expected_rows or meta['columns'] != len(pop)):
                raise ValueError('target source/basis/normalization mismatch')
            _, arrays = receipt_arrays(record)
            if max(abs(arrays['number']-pop.sum()), np.max(abs(arrays['number_pair']-pop.sum()))) > 1e-8*(1+pop.sum()):
                raise ValueError('full raw-column target charge mismatch')
            receipts.append(record)
            yield meta['slot'], arrays
    aggregate = aggregate_complete(control.grid, records())
    ref = control.reference
    volume = abs(np.linalg.det(ref.cell_bohr))
    integrated = volume/len(aggregate['density_bohr3'])*aggregate['density_bohr3'].sum(axis=0)
    expected_number = np.zeros(2)
    for _, slot, pop in rows:
        expected_number[slot.ispin-1] += slot.weight*pop.sum()
    if max(np.max(abs(integrated-expected_number)),
           np.max(abs(aggregate['weighted_spin_number']-expected_number))) > 1e-8*(1+expected_number.max()):
        raise ValueError('exact original weighted spin charge mismatch')
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
    archive = save_arrays(args.array_root.resolve()/(args.case+'-matched-reference.npz'), arrays, args.max_array_bytes)
    diagnostics = dict(number_residual_max=float(np.max(abs(integrated-expected_number))),
        force_imaginary_max=float(np.max(abs(arrays['force_at_origin_ry_bohr'].imag))),
        onsite_hermiticity_max=float(np.max(abs(arrays['electronic_onsite_matrix_ry_bohr2']-
                                               arrays['electronic_onsite_matrix_ry_bohr2'].conj().T))),
        field_contact_hermiticity_max=float(np.max(abs(arrays['contact']-arrays['contact'].conj().T))),
        force_norm=float(np.linalg.norm(arrays['force_at_origin_ry_bohr'])),
        lattice=lattice_residuals(arrays['common_force_lattice_ry_bohr2']))
    result = dict(status='complete', model=reg['model'], registration=reg_binding, context=context,
        archive=archive, source_receipts=receipts, slot_count=len(receipts),
        original_source_count=len({(s.ispin,s.source_ik) for s in control.grid.slots}),
        ionic_reference=dict(authority=linear.ionic_authority,
            quadratic='qe_ionic_curvature on original full signed density G; same alpha/real cutoff; no imported Native Hessian',
            finite_sum_difference=reg['ionic_finite_sum_qualification']),
        diagnostics=diagnostics, scientific_limits=reg['scientific_limits'],
        loaded_independent_evidence=reg['loaded_independent_evidence'],
        loaded_aggregate_coefficients_reused=False,
        exact_loaded_identity_field_cache_slots=sum(
            value['metadata'].get('field_reuse') is not None for _,value in catalog.values()),
        response_comparison=False,
        scalar='per cell: (g_el-F_II).u + 1/2 u.C(q).u + field_gradient.a + 1/2 a.contact.a + u.mixed.a; fixed prepared electrons',
        geometry_policy='direct-original-target; no moment covariance approximation')
    save_json(args.output/'result.json', result)
    print(__import__('json').dumps(dict(result=binding(args.output/'result.json')), sort_keys=True))


if __name__ == '__main__':
    main()
