#!/usr/bin/env python3
"""Reuse a completed original representative; never recompute source coefficients.

--plan writes exact original identity-representative direct commands. The normal
mode consumes direct slot/source-state receipts, authenticates RealBlochBasis
and physical UPF factors, then produces one source orbit. Supplied nonidentity
direct receipts are comparisons only: they NEVER silently replace predictions.
"""
import argparse
import json
from pathlib import Path
import shlex
import time

from matched_reference_io import (HERE, binding, save_json, save_arrays, runtime,
                                  registration, load_case)
from matched_transport_io import (transport_registration, slot_catalog,
                                  authenticated_arrays, source_state_receipt)


def write_plan(args, control, rows, reg_binding, transport_binding):
    from substrate_framework.material_matched_transport import identity_representative
    commands = []
    prefix = ('PYTHONPATH=/home/administrator/.cache/dft/p254-a28-runtime '
              'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 '
              'taskset -c 6 /usr/bin/python3 ')
    for spin, ik in sorted({(s.ispin, s.source_ik) for _, s, _ in rows}):
        s = identity_representative(control.grid, spin, ik)
        directory = args.output.parent/(args.case+'-representatives')/f'source-{spin}-{ik}'
        command = prefix+shlex.quote(str(HERE/'run_matched_reference_source.py'))
        command += (f' --case {args.case} --spin {spin} --source-ik {ik} --target '
                    + ','.join(map(str, s.mesh_index))+' --array-root '
                    + shlex.quote(str(args.array_root))+' --output '+shlex.quote(str(directory)))
        commands.append(dict(spin=spin, source_ik=ik, mesh_index=list(s.mesh_index),
            output=str(directory), command=command,
            loaded_cache='existing exact optical/mixed/state/spatial cache' if args.case == 'loaded-ni31h6' else None))
    save_json(args.output, dict(status='planned-not-executed', case=args.case,
        direct_registration=reg_binding, transport_registration=transport_binding,
        representatives=commands, instruction='Reuse already completed identity receipts; do not replay them.'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', required=True)
    p.add_argument('--spin', type=int, choices=(1, 2))
    p.add_argument('--source-ik', type=int)
    p.add_argument('--target', action='append')
    p.add_argument('--input-dir', action='append', type=Path, default=[])
    p.add_argument('--array-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--oracle', type=Path)
    p.add_argument('--accept-declared-geometry', action='store_true')
    p.add_argument('--plan', action='store_true')
    p.add_argument('--max-array-bytes', type=int, default=1 << 30)
    args = p.parse_args()
    runtime(args.array_root)
    import numpy as np
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.material_matched_reference import slot_identity
    from substrate_framework.material_matched_transport import (identity_representative,
        transport_slot, coefficient_residuals, GEOMETRY_POLICY, require_arrays)
    from substrate_framework.material_matched_transport_certificate import (
        authenticated_basis, initial_state_diagnostics, physical_diagnostics)
    from substrate_framework.material_physical_operators import build_physical_factors
    from substrate_framework.qe_saved_native import load_wavefunctions
    from substrate_framework.qe_bloch_frame import admit_source
    reg_binding, reg = registration()
    transport_binding, _ = transport_registration()
    control, background, ledger, rows, context = load_case(args.case, reg, background=not args.plan)
    if args.plan:
        write_plan(args, control, rows, reg_binding, transport_binding)
        return
    if not args.accept_declared_geometry or args.oracle is None:
        p.error('--accept-declared-geometry and a completed --oracle are required')
    oracle_binding = binding(args.oracle)
    oracle = bound_json(oracle_binding)
    if oracle['status'] != 'passed' or oracle['registration'] != transport_binding:
        raise ValueError('passed independent source-bound transport oracle required')
    orbit = [row for row in rows if (row[1].ispin, row[1].source_ik) == (args.spin, args.source_ik)]
    if not orbit or any(not np.array_equal(orbit[0][2], row[2]) for row in orbit):
        raise ValueError('one full original source orbit with identical populations required')
    representative = identity_representative(control.grid, args.spin, args.source_ik)
    row_map = {(s.ispin, s.mesh_index): (i, s, pop) for i, s, pop in rows}
    rep_key = (representative.ispin, representative.mesh_index)
    catalog = slot_catalog(args.input_dir, context, reg_binding)
    if rep_key not in catalog:
        raise ValueError('completed DIRECT ORIGINAL representative receipt required; use --plan')
    rep_record, rep_value = catalog[rep_key]
    coefficients = authenticated_arrays(catalog[rep_key], row_map[rep_key], ledger)
    require_arrays(coefficients, control.reference.nat, int(np.prod(control.reference.fft_smooth)))
    selected = orbit
    if args.target and args.target != ['all']:
        targets = [tuple(map(int, item.split(','))) for item in args.target]
        selected = [row for row in orbit if row[1].mesh_index in targets]
        if len(set(targets)) != len(targets) or len(selected) != len(targets):
            raise ValueError('unique original orbit targets required')
    source = load_wavefunctions(representative.source_binding)
    admission = admit_source(control.reference, representative, source)
    state_binding, state = source_state_receipt(args.input_dir, context, reg_binding,
        args.spin, args.source_ik, orbit[0][2], source)
    rep_basis = authenticated_basis(control, representative, source, admission,
                                    str(ledger['basis_identity'][row_map[rep_key][0]]))
    source_factors = build_physical_factors(background.inputs, rep_basis)
    context = dict(context, registration=reg_binding, transport_registration=transport_binding,
        original_source=dict(representative.source_binding), spin=args.spin, source_ik=args.source_ik,
        representative_receipt=rep_record, state_receipt=state_binding, oracle=oracle_binding)
    args.output.mkdir(parents=True, exist_ok=False)
    array_dir = args.array_root.resolve()/args.case/f'source-{args.spin}-{args.source_ik}'
    array_dir.mkdir(parents=True, exist_ok=True)
    receipts, comparisons = [], []
    for i, slot, population in sorted(selected, key=lambda row: row[1].mesh_index):
        began = time.monotonic()
        key = (slot.ispin, slot.mesh_index)
        basis = authenticated_basis(control, slot, source, admission, str(ledger['basis_identity'][i]))
        if key == rep_key:
            values, meta = coefficients, dict(rep_value['metadata'])
            meta.update(coefficient_origin='direct-original-representative', direct_receipt=rep_record)
        else:
            values, meta = transport_slot(control.grid, slot, representative, coefficients,
                                           geometry_policy=GEOMETRY_POLICY)
            meta.update(basis_identity=basis.identity, frame_identity=basis.identity+':raw-saved-columns',
                        rows=slot.expected_rows, columns=len(population))
        meta['initial_state'] = initial_state_diagnostics(state, slot)
        meta['physical_diagnostics'] = physical_diagnostics(control, background.inputs,
            slot, source, basis, source_factors)
        if key in catalog and key != rep_key:
            direct = authenticated_arrays(catalog[key], row_map[key], ledger)
            comparison = dict(slot=slot_identity(slot), direct_receipt=catalog[key][0],
                residuals=coefficient_residuals(values, direct),
                tolerance_license=False, transported_values_replaced_by_direct=False)
            meta['direct_comparison'] = comparison
            comparisons.append(comparison)
        name = 'slot-'+'-'.join(map(str, slot.mesh_index))
        archive = save_arrays(array_dir/(name+'.npz'), values, args.max_array_bytes)
        receipt = save_json(args.output/(name+'.json'), dict(status='complete', context=context,
            metadata=meta, archive=archive, elapsed_seconds=time.monotonic()-began))
        receipts.append(receipt)
    result = dict(status='complete-selected-source-orbit', context=context, slots=receipts,
        original_representative=slot_identity(representative), complete_source_orbit=len(selected)==len(orbit),
        complete_case=False, direct_comparisons=comparisons,
        source_coefficients_recomputed=False, positive_time_symmetry=False,
        geometry_policy=GEOMETRY_POLICY, approximation_error_budget=None)
    print(json.dumps(dict(result=save_json(args.output/'result.json', result)), sort_keys=True))


if __name__ == '__main__':
    main()
