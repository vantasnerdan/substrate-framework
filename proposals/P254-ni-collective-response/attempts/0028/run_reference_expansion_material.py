#!/usr/bin/env python3
"""Reuse actual A28 spin1/sourceik2; compare one transported slot and curvature.

Root allocates execution. All 18569 physical PW rows and 346 raw columns are
retained. No H0/SCF recalculation, no full-material or converged-reference claim.
A fresh --output directory and explicit archive capacity are required.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--output-max-bytes', type=int, required=True)
    parser.add_argument('--block-size', type=int, default=8)
    parser.add_argument('--target-mesh-slot', type=int)
    parser.add_argument('--oracle', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    began = time.monotonic()
    summary = dict(status='failed', full_material_total=False, timings_seconds={})
    try:
        import numpy as np
        from material_run_contract import binding, bound_json, canonical, recover_block, implementation, runtime_context
        from material_pilot_io import commit_block
        from material_reference_output import load_arrays
        from substrate_framework.material_bloch_output import preload_parsers
        from substrate_framework.material_preparation_output import load_saved_material_preparation
        from substrate_framework.material_bloch_admission import restore_accepted, wave_bindings
        from substrate_framework.qe_saved_reference import load_saved_reference
        from substrate_framework.qe_saved_native import load_wavefunctions
        from substrate_framework.qe_bloch_grid import build_bloch_grid
        from substrate_framework.qe_bloch_frame import transport_bloch_frame
        from substrate_framework.qe_bloch_dyads import BlochDyadCoupling
        from substrate_framework.material_reference_state import (
            source_population_indices, restore_source, sector_density_moments)
        from substrate_framework import material_reference_transport as transport
        from substrate_framework import material_frozen_state_curvature as curvature
        if args.output_max_bytes <= 65536:
            raise ValueError('explicit output capacity above receipt reserve required')
        preload_parsers()
        summary['runtime_context'] = runtime_context()
        with (args.output/'runtime-context.json').open('x') as stream:
            json.dump(summary['runtime_context'], stream, indent=2)
        registration_path = HERE/'reference-expansion-registration.json'
        registration_binding = binding(registration_path)
        registration = bound_json(registration_binding)
        records = {item['path']: item for item in registration['read_source_bindings']}
        for record in records.values():
            if binding(record['path']) != record:
                raise ValueError('registered dependency changed: '+record['path'])
        oracle_binding = binding(args.oracle.resolve())
        oracle = bound_json(oracle_binding)
        if oracle['status'] != 'passed' or oracle['material_execution'] is not False:
            raise ValueError('completed independent reference-expansion oracle required')
        for record in oracle['implementation_sources']:
            if binding(record['path']) != record:
                raise ValueError('passed oracle belongs to changed implementation: '+record['path'])
        code_paths = (Path(__file__), Path(transport.__file__), Path(curvature.__file__), registration_path)
        sources = implementation(code_paths)
        state_binding = records[str(HERE/'material-reference-source-001/source-state.json')]
        expected = bound_json(state_binding)
        metadata = expected['metadata']
        if (metadata['spin'], metadata['source_ik'], metadata['all_raw_columns']) != (1, 2, True):
            raise ValueError('exact completed full-column spin1/sourceik2 state required')
        archive = Path(expected['archive']['path'])
        if recover_block(archive.parent, archive.stem, expected['recipe_sha256'], metadata) != expected:
            raise ValueError('source-state archive differs from bound completed checkpoint')
        source_arrays = load_arrays(expected)
        saved = load_saved_material_preparation(metadata['preparation_terminal'])
        bloch = bound_json(saved.recipe['inputs']['bloch_terminal'])['recipe']
        original = bloch['inputs']
        reference = load_saved_reference(original['xml'])
        grid = build_bloch_grid(reference, wave_bindings(original),
                                symmetry_policy=bloch['settings']['symmetry_policy'])
        populations, pairs = source_population_indices(grid, saved.population_ledger, 1, 2)
        representatives = [(i, s) for i, s in pairs if np.array_equal(s.operation.rotation, np.eye(3, dtype=int))
                           and not np.any(s.operation.fractional_translation) and not np.any(s.reciprocal_shift)]
        if len(representatives) != 1:
            raise ValueError('unique exact identity representative required')
        population_index, representative = representatives[0]
        source = load_wavefunctions(representative.source_binding)
        frame = transport_bloch_frame(reference, representative, source)
        authority = ('A27 neutral source recipe '+saved.recipe['settings']['preparation_recipe']
                     +'; '+saved.coverage_identity)
        if (dict(source.source) != metadata['source']
                or frame.gauge_identity != metadata['original_frame_identity']
                or frame.basis.identity != metadata['original_basis_identity']
                or saved.coverage_identity != metadata['coverage_identity']
                or authority != metadata['state_authority']
                or not np.array_equal(populations, source_arrays['populations'])):
            raise ValueError('source, state, original basis, populations or authority differs')
        state = restore_source(source_arrays, source_identity=source.source['sha256'])
        gradient_binding = records[str(HERE/'material-reference-source-002/slot-1-density-gradient.json')]
        checkpoint = bound_json(gradient_binding)
        gm = checkpoint['metadata']
        if (any(gm.get(key) != value for key, value in metadata.items())
                or gm['basis_identity'] != frame.basis.identity or gm['gauge_identity'] != frame.gauge_identity
                or gm['mesh_index'] != list(representative.mesh_index)
                or gm['measure'] != representative.weight or gm['population_index'] != population_index):
            raise ValueError('density/gradient checkpoint differs from exact source and state')
        ga = Path(checkpoint['archive']['path'])
        if recover_block(ga.parent, ga.stem, checkpoint['recipe_sha256'], gm) != checkpoint:
            raise ValueError('density/gradient archive differs from completed checkpoint')
        moments = load_arrays(checkpoint)
        candidates = [(i, s) for i, s in pairs if s is not representative
                      and (args.target_mesh_slot is None or s.mesh_slot == args.target_mesh_slot)]
        if not candidates:
            raise ValueError('requested distinct same-source slot is unavailable')
        target_index, slot = candidates[0]
        transform = transport.slot_transport(grid, slot)
        target_frame, target_state = transport.transport_source(
            grid, slot, source, state, source_admission=frame.basis.source_admission)
        background, background_receipt = restore_accepted(original)
        summary['timings_seconds']['authenticate_restore'] = time.monotonic()-began
        started = time.monotonic()
        coupling = BlochDyadCoupling(background.inputs, frame.basis, block_size=args.block_size)
        correction = curvature.conditional_minus_mv(frame.coefficients, state, coupling,
            frame.raw_occupations, conditional_density=moments['density_bohr3'])
        # Keep actual local reference separate. A selected source cannot replace
        # full conditional density and must never subtract the full saved local
        # reference once per source.
        saved_local = curvature.local_density_curvature(background.rho_updown.sum(axis=1), coupling)
        source_local = curvature.local_density_curvature(moments['density_bohr3'], coupling)
        summary['timings_seconds']['source_full_column_curvature'] = time.monotonic()-started
        # The curvature is small; release representative projectors before the
        # full target density/gradient contraction allocates target projectors.
        del coupling
        started = time.monotonic()
        target_coupling = BlochDyadCoupling(background.inputs, target_frame.basis,
                                           block_size=args.block_size)
        direct = sector_density_moments(target_frame, target_state, target_coupling)
        predicted_density = transform.density(moments['density_bohr3'])
        predicted_gradient = transform.gradient(moments['positive_gradient_ry_bohr'])
        density_residual = direct['density_bohr3']-predicted_density
        gradient_residual = direct['positive_gradient_ry_bohr']-predicted_gradient
        # Measure density-G support closure independently from geometry. The
        # inherited PW target support already follows checked frame transport.
        gmiller = background.inputs.density.miller
        support = {tuple(row) for row in gmiller.T}
        rotated_support = {tuple(row) for row in (slot.operation.rotation @ gmiller).T}
        diagnostic = transform.diagnostics(background)
        diagnostic['density_G_support_symmetric_difference_count'] = len(support ^ rotated_support)
        def residual_record(value, reference_value):
            return dict(max_abs=float(np.max(np.abs(value))), frobenius=float(np.linalg.norm(value)),
                        reference_frobenius=float(np.linalg.norm(reference_value)))
        diagnostic['density_direct_minus_transport'] = residual_record(density_residual, direct['density_bohr3'])
        diagnostic['gradient_direct_minus_transport'] = residual_record(gradient_residual, direct['positive_gradient_ry_bohr'])
        summary['timings_seconds']['target_direct_density_gradient_and_discrepancies'] = time.monotonic()-started
        # Density is exactly the declared finite-frame pullback; force covariance
        # residual is retained, not clipped, fitted or declared equal by a tolerance.
        if np.max(np.abs(density_residual)) > 1e-8*(1+np.max(np.abs(direct['density_bohr3']))):
            raise ValueError('checked finite-frame density pullback failed')
        arrays = dict(correction,
            source_conditional_local_ry_bohr2=source_local,
            saved_full_density_local_ry_bohr2=saved_local,
            source_positive_gradient_ry_bohr=moments['positive_gradient_ry_bohr'],
            target_direct_density_bohr3=direct['density_bohr3'],
            target_transport_density_bohr3=predicted_density,
            target_density_residual_bohr3=density_residual,
            target_direct_positive_gradient_ry_bohr=direct['positive_gradient_ry_bohr'],
            target_transport_positive_gradient_ry_bohr=predicted_gradient,
            target_gradient_residual_ry_bohr=gradient_residual,
            target_covariance_conditional_nonlocal_ry_bohr2=transform.onsite(correction['conditional_nonlocal_ry_bohr2']),
            original_source_xml_weight=np.asarray(reference.select(representative.source_ik, representative.ispin)[0].weight),
            source_density_coefficients=state.density_coefficients, source_gram=state.gram,
            raw_mv_numerical_reference_coefficients=frame.raw_occupations,
            geometry_residual_bohr=transform.geometry_residual_bohr,
            source_mesh_index=np.asarray(representative.mesh_index), target_mesh_index=np.asarray(slot.mesh_index),
            source_slot_measure=np.asarray(representative.weight), target_slot_measure=np.asarray(slot.weight))
        recipe = dict(registration=registration_binding, oracle=oracle_binding,
                      source_checkpoint=state_binding, moments_checkpoint=gradient_binding,
                      implementation_sources=sources, block_size=args.block_size,
                      source_basis_identity=frame.basis.identity, target_basis_identity=target_frame.basis.identity)
        result_metadata = dict(schema='material-reference-expansion-selected-source-v1', recipe=recipe,
            source_metadata=metadata, background_receipt=background_receipt,
            source_population_index=population_index, target_population_index=target_index,
            all_raw_columns=True, full_physical_pw_rows=True,
            raw_column_count=source.nbnd, physical_pw_rows=source.igwx,
            source_slot_measure_applied=False, target_slot_measure_applied=False,
            diagnostics=diagnostic, full_material_total=False,
            curvature_scope='source unweighted conditional and original numerical-MV nonlocal blocks; source delta and orbital local difference diagnostic only',
            common_reference_rule='local[sum_complete_weighted_conditional_density - saved_full_QE_density] + sum_full_slot_weighted_conditional_nonlocal - sum_original_source_XML_weighted_MV_nonlocal',
            ionic_scope='A22 ionic-q sum/alpha/cutoffs/phases unchanged; no ionic addition in this output',
            h0_reused=False, units=dict(density='bohr^-3', positive_gradient='Ry/bohr', curvature='Ry/bohr^2'))
        for record in sources:
            if binding(record['path']) != record:
                raise ValueError('implementation changed during selected-source execution')
        identity = hashlib.sha256(canonical(recipe)).hexdigest()
        receipt = commit_block(args.output, 'reference-expansion', arrays, result_metadata,
                               identity, args.output_max_bytes)
        summary.update(status='complete', output=receipt, diagnostics=diagnostic,
            source_checkpoint=state_binding, moments_checkpoint=gradient_binding,
            source_slot_measure=representative.weight, target_slot_measure=slot.weight,
            scope=result_metadata['curvature_scope'], common_reference_rule=result_metadata['common_reference_rule'])
    except BaseException as error:
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    summary['timings_seconds']['total'] = time.monotonic()-began
    with (args.output/'summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=summary['status'], summary=str(args.output/'summary.json'))), flush=True)
    return 0 if summary['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
