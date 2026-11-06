#!/usr/bin/env python3
"""Combine completed original-source checkpoints into one complete conditional reference.

No WFC/Gram/density/H0 recalculation. One original frame supplies the already
checked local-curvature adapter. All128 slots contribute; H0 transport is not
inferred. The geometric force/onsite approximation and ionic-sum mismatch remain
explicit. This is not equilibrium, a time trajectory or a prepared ionic state.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--output-max-bytes', type=int, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    summary = dict(status='failed', full_material_total=False)
    def interrupted(signum, frame):
        raise InterruptedError('signal '+str(signum))
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, interrupted)
    try:
        import numpy as np
        from material_run_contract import binding, bound_json, canonical, implementation, runtime_context
        from material_pilot_io import commit_block
        from substrate_framework.material_bloch_output import preload_parsers
        from substrate_framework.material_preparation_output import load_saved_material_preparation
        from substrate_framework.material_ionic_output import load_saved_ionic_linear_term
        from substrate_framework.material_bloch_admission import restore_accepted, wave_bindings
        from substrate_framework.qe_saved_reference import load_saved_reference
        from substrate_framework.qe_saved_native import read_bound_bytes, load_wavefunctions
        from substrate_framework.qe_bloch_grid import build_bloch_grid
        from substrate_framework.qe_bloch_frame import transport_bloch_frame
        from substrate_framework.qe_bloch_dyads import BlochDyadCoupling
        from substrate_framework.material_reference_transport import aggregate_complete_slots
        from substrate_framework.material_frozen_state_curvature import common_reference_correction
        from substrate_framework.material_force_lattice import q_grid_to_force_lattice, force_lattice_to_q_grid
        import full_reference_inputs as source_reader
        preload_parsers()
        manifest_binding = binding(args.manifest)
        if manifest_binding['sha256'] != args.manifest_sha256:
            raise ValueError('explicit complete-manifest digest differs')
        manifest = bound_json(manifest_binding)
        if manifest['gradient_policy'] != 'declared-geometric':
            raise ValueError('this aggregate requires the explicitly measured geometric approximation')
        oracle = bound_json(manifest['oracle'])
        if oracle['status'] != 'passed':
            raise ValueError('completed independent source/transport/curvature oracle required')
        for record in oracle['implementation_sources']:
            if binding(record['path']) != record:
                raise ValueError('reference-expansion oracle implementation changed')
        sources = implementation((Path(__file__), Path(source_reader.__file__)))
        context = runtime_context()
        with (args.output/'runtime-context.json').open('x') as stream:
            json.dump(context, stream, indent=2)
        saved = load_saved_material_preparation(manifest['preparation_terminal'])
        bloch = bound_json(saved.recipe['inputs']['bloch_terminal'])['recipe']
        original = bloch['inputs']
        reference = load_saved_reference(original['xml'])
        grid = build_bloch_grid(reference, wave_bindings(original),
                                symmetry_policy=bloch['settings']['symmetry_policy'])
        diagnostics = []
        aggregate = aggregate_complete_slots(grid, saved.population_ledger,
            source_reader.source_observables(manifest, saved, grid, diagnostics),
            gradient_policy=manifest['gradient_policy'])
        summary['source_reduction_seconds'] = time.monotonic()-started
        background, background_receipt = restore_accepted(original)
        ionic = load_saved_ionic_linear_term(manifest['ionic_terminal'])
        if (not np.array_equal(ionic.arrays['cell_bohr'], reference.cell_bohr)
                or not np.array_equal(ionic.arrays['positions_bohr'], reference.positions_bohr.T)):
            raise ValueError('ionic and electronic expansion geometries differ')
        frozen = bound_json(manifest['frozen_q_result'])
        frozen_manifest = bound_json(frozen['source_manifest'])
        if frozen_manifest != frozen['sources']:
            raise ValueError('frozen-grid manifest differs from its bound original')
        for key in ('xml', 'density', 'upfs'):
            if frozen_manifest['inputs'][key] != original[key]:
                raise ValueError('frozen and conditional references use different '+key)
        if frozen['atom_names'] != list(reference.atom_names):
            raise ValueError('frozen grid atom identities differ')
        raw, _ = read_bound_bytes(frozen['archive'], max_bytes=1 << 30)
        with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
            frozen_grid = archive['frozen_q_grid']
            masses = archive['atom_masses_amu']
        shape = tuple(reference.monkhorst_pack['nk'+str(i)] for i in (1, 2, 3))
        if frozen_grid.shape != shape+(3*reference.nat, 3*reference.nat):
            raise ValueError('complete frozen q grid has a different shape')
        representative = next(slot for slot in grid.slots if slot.ispin == 1
                              and slot.source_ik == 1 and slot.operation.index == 1
                              and not np.any(slot.reciprocal_shift))
        wave = load_wavefunctions(representative.source_binding)
        frame = transport_bloch_frame(reference, representative, wave)
        coupling = BlochDyadCoupling(background.inputs, frame.basis, block_size=8)
        correction = common_reference_correction(aggregate['density_bohr3'], background.rho_updown,
                                                 aggregate['nonlocal_correction_ry_bohr2'], coupling)
        common_grid = frozen_grid+correction.matrix()
        lattice = q_grid_to_force_lattice(common_grid)
        number = coupling.grid_measure_bohr3*np.sum(aggregate['density_bohr3'], axis=0)
        expected_number = np.zeros(2)
        for item in diagnostics:
            expected_number[item['spin']-1] += item['total_slot_measure']*item['population_number']
        if np.max(np.abs(number-expected_number)) > 1e-8*(1+np.max(np.abs(expected_number))):
            raise ValueError('full weighted density does not recover the declared prepared charge')
        kernel = background.apply_kernel(aggregate['density_bohr3']-background.rho_updown)
        arrays = dict(conditional_density_bohr3=aggregate['density_bohr3'],
            saved_density_bohr3=background.rho_updown, mv_density_bohr3=aggregate['mv_density_bohr3'],
            positive_gradient_ry_bohr=aggregate['positive_gradient_ry_bohr'],
            ionic_force_ry_bohr=ionic.arrays['ionic_force_ry_bohr'],
            force_at_origin_ry_bohr=ionic.arrays['ionic_force_ry_bohr']-aggregate['positive_gradient_ry_bohr'],
            kernel_potential_ry=kernel,
            nonlocal_correction_ry_bohr2=aggregate['nonlocal_correction_ry_bohr2'],
            local_correction_ry_bohr2=correction.local, onsite_correction_ry_bohr2=correction.total,
            conditional_nonlocal_ry_bohr2=aggregate['conditional_nonlocal_ry_bohr2'],
            original_mv_nonlocal_ry_bohr2=aggregate['original_mv_nonlocal_ry_bohr2'],
            common_q_grid_ry_bohr2=common_grid, common_force_lattice_ry_bohr2=lattice,
            atomic_masses_amu=masses, cell_bohr=reference.cell_bohr, positions_bohr=reference.positions_bohr,
            integrated_spin_number=number, weighted_population_spin_number=expected_number,
            source_spin=np.asarray([d['spin'] for d in diagnostics], dtype=np.int64),
            source_ik=np.asarray([d['source_ik'] for d in diagnostics], dtype=np.int64),
            representative_energy_ry=np.asarray([d.pop('energy') for d in diagnostics]),
            representative_h0_commutator_norm2_ry2=np.asarray([d.pop('h0_commutator_norm_squared') for d in diagnostics]))
        report = dict(number_residual_max=float(np.max(np.abs(number-expected_number))),
            force_at_origin_norm=float(np.linalg.norm(arrays['force_at_origin_ry_bohr'])),
            force_at_origin_imaginary_max=float(np.max(np.abs(arrays['force_at_origin_ry_bohr'].imag))),
            kernel_potential_norm=float(np.linalg.norm(kernel)),
            density_minus_saved_norm=float(np.linalg.norm(arrays['conditional_density_bohr3']-background.rho_updown)),
            mv_density_minus_saved_norm=float(np.linalg.norm(arrays['mv_density_bohr3']-background.rho_updown)),
            onsite_correction_norm=float(np.linalg.norm(correction.total)),
            full_force_lattice_roundtrip_max=float(np.max(np.abs(force_lattice_to_q_grid(lattice)-common_grid))))
        metadata = dict(schema='material-full-conditional-reference-v1', recipe=manifest_binding,
            source_records=manifest['sources'], coverage_identity=saved.coverage_identity,
            slots=aggregate['slots'], source_diagnostics=diagnostics, runtime_context=context,
            background_receipt=background_receipt, ionic_recipe=dict(ionic.recipe),
            ionic_sum_metadata=dict(ionic.term.sum_metadata),
            frozen_ionic_alpha=frozen['alpha_binary_reconstructed'],
            gradient_policy=aggregate['gradient_policy'], complete_electronic_coverage=True,
            held_density='zero because every admitted spin/mesh slot is active',
            held_gradient='zero because every admitted spin/mesh slot is active',
            support_force='none; no fitted force cancellation',
            h0_transport_supported=False, report=report, implementation_sources=sources,
            limitations=['Conditional initial electrons, not equilibrium or experimental temperature.',
                'Force and onsite blocks use declared geometric covariance with retained saved-geometry residuals.',
                'H0 diagnostics are original-source only; initial common-Hamiltonian drift includes K(delta n) and is not set to zero.',
                'Linear ionic coefficient is the accepted A27 finite scalar; frozen ionic curvature is the A22 finite sum. Their different alpha/cutoffs are retained, not claimed identical finite-Ewald derivatives.',
                'Frozen q curvature is not electronically screened phonon response; no ASR/Hermitian repair.',
                'Ionic and field quantum preparation, coupled trajectory and coherence claim remain separate.'])
        for record in sources:
            if binding(record['path']) != record:
                raise ValueError('aggregate consumer changed during execution')
        recipe_identity = hashlib.sha256(canonical(dict(manifest=manifest_binding, implementation=sources))).hexdigest()
        receipt = commit_block(args.output, 'full-common-reference', arrays, metadata,
                               recipe_identity, args.output_max_bytes)
        summary.update(status='complete', full_material_total=True, output=receipt,
                       source_count=len(diagnostics), slot_count=len(aggregate['slots']), report=report)
    except BaseException as error:
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    summary['elapsed_seconds'] = time.monotonic()-started
    with (args.output/'summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=summary['status'], summary=str(args.output/'summary.json'))), flush=True)
    print('[dft-done] A28 full reference' if summary['status'] == 'complete' else '[dft-fail] A28 full reference', flush=True)
    return 0 if summary['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
