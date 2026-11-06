#!/usr/bin/env python3
"""Actual A28 spin1/sourceik2 initial electronic force noise; Root allocates runtime.

Restores accepted A25 background and the exact checkpoint's WFC/frame. No SCF,
new population recipe, H0 stationarity assumption or whole-material weighting.
Use a fresh --output directory. Actions use temporary disk-backed full frames.
"""
import time
BOOT = time.monotonic()
import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scratch-directory', type=Path)
    parser.add_argument('--block-size', type=int, default=8)
    parser.add_argument('--row-block-size', type=int, default=2048)
    parser.add_argument('--output-max-bytes', type=int, required=True)
    args = parser.parse_args()
    target = args.output.resolve()
    target.mkdir(parents=True, exist_ok=False)
    summary = dict(status='failed', scope='one unweighted electronic source sector, equal time only',
                   timings_seconds={}, outputs=[])
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
        from substrate_framework.material_reference_state import source_population_indices, restore_source
        from substrate_framework import material_force_fluctuations as noise
        if args.output_max_bytes <= 65536:
            raise ValueError('explicit output capacity greater than receipt reserve required')
        preload_parsers()
        summary['runtime_context'] = runtime_context()
        with (target/'runtime-context.json').open('x') as stream:
            json.dump(summary['runtime_context'], stream, indent=2)
        began = time.monotonic()
        registration_path = HERE/'force-noise-registration.json'
        registration_binding = binding(registration_path)
        registration = bound_json(registration_binding)
        records = {item['path']: item for item in registration['read_source_bindings']}
        # Bind original sources before the first material operation. These are
        # read-only scientific dependencies, not a replacement runtime grant.
        for record in records.values():
            if binding(record['path']) != record:
                raise ValueError('registered original source changed: '+record['path'])
        receipt_binding = records[str(HERE/'material-reference-source-001/source-state.json')]
        expected = bound_json(receipt_binding)
        metadata = expected['metadata']
        if (metadata['spin'], metadata['source_ik'], metadata['all_raw_columns']) != (1, 2, True):
            raise ValueError('exact full-column spin1/sourceik2 checkpoint required')
        archive = Path(expected['archive']['path'])
        actual = recover_block(archive.parent, archive.stem, expected['recipe_sha256'], metadata)
        if actual != expected:
            raise ValueError('source-state archive differs from original-byte receipt')
        source_arrays = load_arrays(actual)
        saved = load_saved_material_preparation(metadata['preparation_terminal'])
        bloch = bound_json(saved.recipe['inputs']['bloch_terminal'])['recipe']
        original = bloch['inputs']
        reference = load_saved_reference(original['xml'])
        grid = build_bloch_grid(reference, wave_bindings(original),
                                symmetry_policy=bloch['settings']['symmetry_policy'])
        populations, pairs = source_population_indices(grid, saved.population_ledger, 1, 2)
        matches = [(i, s) for i, s in pairs if s.operation.index == 1 and not np.any(s.reciprocal_shift)]
        if len(matches) != 1:
            raise ValueError('unique original source representative required')
        population_index, slot = matches[0]
        source = load_wavefunctions(slot.source_binding)
        frame = transport_bloch_frame(reference, slot, source)
        authority = ('A27 neutral source recipe '+saved.recipe['settings']['preparation_recipe']
                     +'; '+saved.coverage_identity)
        if (dict(source.source) != metadata['source']
                or frame.gauge_identity != metadata['original_frame_identity']
                or frame.basis.identity != metadata['original_basis_identity']
                or saved.coverage_identity != metadata['coverage_identity']
                or authority != metadata['state_authority']
                or not np.array_equal(populations, source_arrays['populations'])
                or frame.coefficients.shape != (source.igwx, len(populations))):
            raise ValueError('checkpoint source, raw frame, basis, population or state authority differs')
        state = restore_source(source_arrays, source_identity=source.source['sha256'])
        gradient_binding = records[str(HERE/'material-reference-source-002/slot-1-density-gradient.json')]
        gradient_receipt = bound_json(gradient_binding)
        gm = gradient_receipt['metadata']
        if (any(gm.get(key) != value for key, value in metadata.items())
                or gm['basis_identity'] != frame.basis.identity
                or gm['gauge_identity'] != frame.gauge_identity
                or gm['mesh_index'] != list(slot.mesh_index)
                or gm['measure'] != slot.weight or gm['population_index'] != population_index):
            raise ValueError('completed gradient source/state/slot differs')
        ga = Path(gradient_receipt['archive']['path'])
        if recover_block(ga.parent, ga.stem, gradient_receipt['recipe_sha256'], gm) != gradient_receipt:
            raise ValueError('completed gradient archive differs')
        with np.load(ga, allow_pickle=False) as stored:
            gradient = stored['positive_gradient_ry_bohr']
            number_pair = stored['number_pair']
        summary['timings_seconds']['authenticate_restore_state_frame'] = time.monotonic()-began
        began = time.monotonic()
        background, background_receipt = restore_accepted(original)
        summary['timings_seconds']['restore_accepted_A25_background'] = time.monotonic()-began
        began = time.monotonic()
        coupling = BlochDyadCoupling(background.inputs, frame.basis, block_size=args.block_size)
        if reference.atom_names[31] != 'H' or reference.atom_names[0] != 'Ni':
            raise ValueError('requested H32-x/Ni1-y physical atom identities differ')
        directions = {}
        for name, index in (('H32-x', 3*31), ('Ni1-y', 1)):
            direction = np.zeros(coupling.ndof)
            direction[index] = 1.
            directions[name] = direction
        coordinate_identity = 'atom-major-Cartesian-bohr:'+original['xml']['sha256']
        original_sources = implementation((Path(__file__), Path(noise.__file__), registration_path))
        summary['timings_seconds']['construct_checked_ionic_coupling'] = time.monotonic()-began
        began = time.monotonic()
        result = noise.initial_force_fluctuations(frame.coefficients, state.density_coefficients,
            coupling, directions, factor_identity=frame.gauge_identity,
            pw_basis_identity=frame.basis.identity, state_authority=authority,
            operator_authority='accepted-A25:'+background_receipt['archive']['sha256']
                +'; XML:'+original['xml']['sha256']+'; BlochDyadCoupling positive ionic gradient',
            coordinate_identity=coordinate_identity, gram=state.gram,
            gram_authority=receipt_binding['sha256']+':gram:'+expected['arrays']['gram']['sha256_C_bytes'],
            scratch_directory=args.scratch_directory, block_size=args.block_size,
            row_block_size=args.row_block_size)
        arrays, result_metadata = result.archive_payload()
        arrays.update(populations=state.populations,
                      checkpoint_density_eigenvalues=source_arrays['density_eigenvalues'],
                      mesh_index=np.asarray(slot.mesh_index, dtype=np.int64),
                      source_slot_measure=np.asarray(slot.weight))
        expected_means = np.asarray([-directions[name] @ gradient for name in directions]
                                   + [number_pair[1]])
        arrays['independent_saved_gradient_number_means'] = expected_means
        arrays['saved_moment_pairing_residual'] = arrays['means']-expected_means
        # Same trace-pair allowance as the admitted A28 material producer,
        # retained as a diagnostic; no real-part extraction or state repair.
        if np.any(np.abs(arrays['saved_moment_pairing_residual']) > 1e-8*(1+np.abs(expected_means))):
            raise ValueError('force sign/number means disagree with completed independent gradient moments')
        if not np.array_equal(arrays['density_coefficients'], source_arrays['density_coefficients']):
            raise ValueError('consumer changed exact checkpoint density coefficients')
        summary['timings_seconds']['full_column_actions_and_noise'] = time.monotonic()-began
        final_sources = implementation((Path(__file__), Path(noise.__file__), registration_path))
        initial_map = {entry['path']: entry for entry in original_sources}
        if any(binding(path) != record for path, record in initial_map.items()):
            raise ValueError('implementation source changed during computation')
        recipe = dict(registration=registration_binding, source_state=receipt_binding,
            gradient_receipt=gradient_binding, source=metadata['source'],
            preparation_terminal=metadata['preparation_terminal'],
            block_size=args.block_size, row_block_size=args.row_block_size,
            directions={key: value.tolist() for key, value in directions.items()},
            implementation_sources=final_sources)
        identity = hashlib.sha256(canonical(recipe)).hexdigest()
        result_metadata.update(schema='material-equal-time-force-noise-v1', recipe=recipe,
            source_checkpoint_metadata=metadata, background_receipt=background_receipt,
            mesh_index=list(slot.mesh_index), population_index=population_index,
            source_slot_measure=slot.weight, source_slot_measure_applied=False,
            all_raw_columns=True, full_physical_pw_rows=True,
            raw_column_count=int(frame.coefficients.shape[1]),
            full_material_total=False, timings_seconds=dict(summary['timings_seconds']))
        began = time.monotonic()
        receipt = commit_block(target, 'force-noise', arrays, result_metadata, identity,
                               args.output_max_bytes)
        summary['timings_seconds']['archive_exact_readback_and_commit'] = time.monotonic()-began
        def packed(value):
            value = np.asarray(value)
            return dict(real=value.real.tolist(), imag=value.imag.tolist())
        summary.update(status='complete', outputs=[receipt], operator_order=result_metadata['operator_order'],
            means=packed(arrays['means']), ordered_AB=packed(arrays['ordered_connected']),
            ordered_BA=packed(arrays['reverse_ordered_connected']),
            symmetrized=packed(arrays['symmetrized_connected']), commutator=packed(arrays['commutator']),
            empty_complement_contribution=packed(arrays['empty_complement_contribution']),
            source_state=receipt_binding, gradient_receipt=gradient_binding,
            implementation_sources=final_sources,
            source_slot_measure=slot.weight, source_slot_measure_applied=False,
            scope=result_metadata['scope'])
    except BaseException as error:
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    summary['timings_seconds']['total_process'] = time.monotonic()-BOOT
    with (target/'summary.json').open('x') as stream:
        json.dump(summary, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=summary['status'], summary=str(target/'summary.json'),
                         elapsed_seconds=summary['timings_seconds']['total_process'])), flush=True)
    return 0 if summary['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
