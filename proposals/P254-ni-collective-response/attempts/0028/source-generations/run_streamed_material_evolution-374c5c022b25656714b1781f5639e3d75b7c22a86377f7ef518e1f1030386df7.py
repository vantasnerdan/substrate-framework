#!/usr/bin/env python3
"""Restore the actual complete material and execute streamed mean evolution.

Root supplies the authenticated aggregate receipt and manifest, the initial
ionic state, and explicit numerical/runtime allocation. Inventory mode does
not reconstruct wavefunctions or execute scientific operators. Derivative mode
is an initial RHS only. Mean ions use the common Gamma coefficient; the full
4^3 lattice is preserved in the aggregate, not relabeled as executed dynamics.
"""
import argparse
from contextlib import contextmanager
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0025'), str(HERE.parent/'0024')]

import numpy as np
from material_run_contract import bound_json, recover_block, implementation, runtime_context
from material_reference_output import load_arrays
from substrate_framework.material_preparation_output import load_saved_material_preparation
from substrate_framework.material_bloch_admission import restore_accepted, wave_bindings
from substrate_framework.material_bloch_output import preload_parsers
from substrate_framework.material_reference_state import restore_source, source_population_indices
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.qe_bloch_grid import build_bloch_grid
from substrate_framework.qe_bloch_frame import transport_bloch_frame
from substrate_framework.qe_bloch_dyads import BlochDyadCoupling
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework.prepared_reference import PreparedReference
from substrate_framework.material_modes import AMU_RY, AU_SEC
from substrate_framework.inertial_fermion_dynamics import ionic_masses_ry
from substrate_framework.streamed_material_state import DiskMean, file_binding, authenticate
from substrate_framework.streamed_material_evolution import StreamedMeanEvolution, real_physical


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def explicit_binding(path, expected):
    binding = file_binding(path)
    if binding['sha256'] != expected:
        raise ValueError('explicit input digest differs: '+str(path))
    return binding


def receipt_arrays(binding):
    expected = bound_json(binding)
    archive = Path(expected['archive']['path'])
    actual = recover_block(archive.parent, archive.stem, expected['recipe_sha256'], expected['metadata'])
    if actual != expected:
        raise ValueError('completed receipt/archive differs from authority')
    return load_arrays(actual), actual


def slot_identity(slot):
    return str(slot.ispin)+':'+','.join(str(i) for i in slot.mesh_index)


def restore_inputs(aggregate_binding, manifest_binding):
    arrays, receipt = receipt_arrays(aggregate_binding)
    if receipt['metadata']['recipe'] != manifest_binding:
        raise ValueError('aggregate was not produced from supplied exact manifest')
    manifest = bound_json(manifest_binding)
    saved = load_saved_material_preparation(manifest['preparation_terminal'])
    bloch = bound_json(saved.recipe['inputs']['bloch_terminal'])['recipe']
    original = bloch['inputs']
    reference = load_saved_reference(original['xml'])
    grid = build_bloch_grid(reference, wave_bindings(original), symmetry_policy=bloch['settings']['symmetry_policy'])
    records = manifest['sources']
    expected = set(wave_bindings(original))
    if len(records) != 20 or {(r['spin'], r['source_ik']) for r in records} != expected:
        raise ValueError('all twenty independent source states required')
    ordered, pairs_by_source = [], {}
    for record in records:
        key = (record['spin'], record['source_ik'])
        _, pairs = source_population_indices(grid, saved.population_ledger, *key)
        pairs_by_source[key] = pairs
        ordered.extend(slot for _, slot in pairs)
    if len(ordered) != 128 or len({slot_identity(s) for s in ordered}) != 128:
        raise ValueError('exact complete128-slot coverage required')
    return arrays, receipt, manifest, saved, original, grid, ordered, pairs_by_source


def initial_sectors(manifest, saved, grid, pairs_by_source):
    for record in manifest['sources']:
        key = (record['spin'], record['source_ik'])
        population, _ = source_population_indices(grid, saved.population_ledger, *key)
        pairs = pairs_by_source[key]
        rep = next(s for _, s in pairs if s.operation.index == 1 and not np.any(s.reciprocal_shift))
        source = load_wavefunctions(rep.source_binding)
        original = transport_bloch_frame(grid.reference, rep, source)
        if record['reuse_selected_gamma']:
            if record['state'] is not None or key != (1, 1) or rep.mesh_index != (0, 0, 0):
                raise ValueError('A27 selected Gamma state authority differs')
            state_arrays = {name: saved.selected_gamma[name] for name in
                            ('gram', 'populations', 'density_coefficients', 'density_eigenvalues')}
        else:
            state_arrays, state_receipt = receipt_arrays(record['state'])
            metadata = state_receipt['metadata']
            if (metadata['source']['sha256'] != source.source['sha256']
                    or metadata['original_basis_identity'] != original.basis.identity
                    or metadata['original_frame_identity'] != original.gauge_identity):
                raise ValueError('original source state basis/raw-column gauge differs')
        if not np.array_equal(state_arrays['populations'], population):
            raise ValueError('source populations differ from all explicit ledger slots')
        state = restore_source(state_arrays, source_identity=source.source['sha256'])
        del original
        for _, slot in pairs:
            frame = transport_bloch_frame(grid.reference, slot, source)
            yield (dict(identity=slot_identity(slot), component=slot.ispin-1, weight=slot.weight,
                        basis_identity=frame.basis.identity, gauge_identity=frame.gauge_identity,
                        source_identity=state.source_identity), frame.coefficients, state.density_coefficients)
            del frame
        del source, state, state_arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('aggregate-receipt', 'manifest'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('inventory', 'initialize', 'derivative', 'evolve'), required=True)
    parser.add_argument('--initial-ions', type=Path)
    parser.add_argument('--initial-ions-sha256')
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--checkpoint-sha256')
    parser.add_argument('--oracle', type=Path)
    parser.add_argument('--oracle-sha256')
    parser.add_argument('--physical-interval-ry-inverse', type=float)
    parser.add_argument('--steps', type=int, default=1)
    parser.add_argument('--max-steps', type=int, required=True)
    parser.add_argument('--max-midpoint-iterations', type=int, required=True)
    parser.add_argument('--midpoint-tolerance', type=float, default=1e-11)
    parser.add_argument('--reality-tolerance', type=float, default=1e-10)
    parser.add_argument('--refinement-label', required=True)
    parser.add_argument('--max-h0-actions', type=int, required=True)
    parser.add_argument('--max-output-bytes', type=int, required=True)
    parser.add_argument('--block-size', type=int, default=8)
    args = parser.parse_args()
    if (not 1 <= args.steps <= args.max_steps or args.max_midpoint_iterations < 1
            or args.max_h0_actions < 1 or args.max_output_bytes < 1):
        parser.error('positive bounded numerical/work/storage allocation required')
    if args.mode == 'evolve' and (args.physical_interval_ry_inverse is None
            or not np.isfinite(args.physical_interval_ry_inverse) or args.physical_interval_ry_inverse <= 0):
        parser.error('evolution needs an independently supplied positive physical interval')
    aggregate_binding = explicit_binding(args.aggregate_receipt, args.aggregate_receipt_sha256)
    manifest_binding = explicit_binding(args.manifest, args.manifest_sha256)
    args.output.mkdir(parents=True, exist_ok=False)
    preload_parsers()
    sources_before = implementation((Path(__file__), HERE/'streamed-evolution-registration.json'))
    arrays, receipt, manifest, saved, original, grid, ordered, pairs = restore_inputs(aggregate_binding, manifest_binding)
    frame_bytes = sum(16*s.expected_rows*len(s.eigenvalues_ry) for s in ordered)
    largest = max(16*s.expected_rows*len(s.eigenvalues_ry) for s in ordered)
    d_bytes = sum(16*len(next(s for _, s in pairs[(r['spin'], r['source_ik'])]).eigenvalues_ry)**2
                  for r in manifest['sources'])
    generations = (0 if args.checkpoint else 1) + (args.steps if args.mode == 'evolve' else int(args.mode == 'derivative'))
    # Account for all retained raw frames, unique fixed D, two full density grids
    # in diagnostics, and generous NPY/JSON headers. Not an RSS estimate.
    storage_bound = generations*frame_bytes+d_bytes+4*arrays['saved_density_bohr3'].nbytes+(16 << 20)
    inventory = dict(schema='P254-streamed-evolution-inventory-v1', sectors=128, sources=20,
        frames_payload_bytes=frame_bytes, largest_frame_bytes=largest, fixed_D_payload_bytes=d_bytes,
        planned_output_bound_bytes=storage_bound, output_cap_bytes=args.max_output_bytes,
        max_h0_actions=args.max_h0_actions,
        max_exponential_sector_calls=128*args.steps*(args.max_midpoint_iterations+1) if args.mode == 'evolve' else 0,
        maximum_moment_sector_reductions=128*(1+args.steps*(args.max_midpoint_iterations+2)) if args.mode == 'evolve' else 128,
        maximum_operator_constructions=(128*(1+args.steps*(args.max_midpoint_iterations+3))
                                        if args.mode == 'evolve' else 256 if args.mode == 'derivative' else 0),
        derivative_h0_calls=256 if args.mode == 'derivative' else None,
        memory='one slot full raw frame plus exponential work arrays, D, FFT/projector workspace and global spin grids; never128 resident frames',
        time_units='inverse Ry; physical seconds = inverse_Ry * 2*AU_SEC from accepted material_modes constants',
        physical_interval_ry_inverse=args.physical_interval_ry_inverse,
        physical_interval_fs=None if args.physical_interval_ry_inverse is None else args.physical_interval_ry_inverse*2*AU_SEC*1e15,
        steps=args.steps, max_steps=args.max_steps, refinement_label=args.refinement_label,
        aggregate=aggregate_binding, manifest=manifest_binding,
        limitation='Taylor exponential work depends on actual spectral scale; action cap enforced at runtime; no runtime duration predicted')
    publish(args.output/'inventory.json', inventory)
    if args.mode == 'inventory':
        return
    if storage_bound > args.max_output_bytes:
        raise ValueError('full-coverage payload exceeds explicit output allocation; inventory retained')
    if args.oracle is None or args.oracle_sha256 is None:
        parser.error('an authenticated executed independent oracle is required')
    oracle_binding = explicit_binding(args.oracle, args.oracle_sha256)
    oracle = bound_json(oracle_binding)
    if oracle['status'] != 'passed' or oracle['schema'] != 'P254-streamed-evolution-oracle-v1':
        raise ValueError('streamed evolution independent oracle has not passed')
    for entry in oracle['source_bindings']:
        authenticate(entry)
    required_oracle_paths = {str((ROOT/'src/substrate_framework'/name).resolve()) for name in
                             ('streamed_material_state.py', 'streamed_material_evolution.py')}
    if not required_oracle_paths <= {entry['path'] for entry in oracle['source_bindings']}:
        raise ValueError('oracle does not bind the actual scientific consumers')
    if args.initial_ions is None or args.initial_ions_sha256 is None:
        parser.error('explicit initial ionic coordinates/momenta required even when zero')
    ions_binding = explicit_binding(args.initial_ions, args.initial_ions_sha256)
    ions = bound_json(ions_binding)
    if set(ions)-{'u_bohr', 'p_canonical', 'coordinate_identity', 'support_force_ry_bohr', 'support_authority'}:
        raise ValueError('unknown ionic/external fields; no unimplemented drive accepted')
    ndof = 3*grid.reference.nat
    u, p = np.asarray(ions['u_bohr'], float), np.asarray(ions['p_canonical'], float)
    if u.shape != (ndof,) or p.shape != (ndof,) or not np.isfinite(u).all() or not np.isfinite(p).all():
        raise ValueError('explicit complete atom-major real ionic phase-space required')
    if ions['coordinate_identity'] != 'QE-original-atom-major-xyz-bohr-Ry':
        raise ValueError('wrong original coordinate convention')
    support = np.asarray(ions.get('support_force_ry_bohr', np.zeros(ndof)), float)
    if 'support_force_ry_bohr' in ions and not ions.get('support_authority'):
        raise ValueError('explicit conservative support authority required')
    model_spec = dict(aggregate=aggregate_binding, manifest=manifest_binding,
        preparation=manifest['preparation_terminal'], initial_ions=ions_binding,
        oracle=oracle_binding, implementation=sources_before,
        scalar='PreparedReference joint-second-order; common Gamma coefficient only; actual per-slot H0')
    model_identity = sha256(json.dumps(model_spec, sort_keys=True).encode()).hexdigest()
    if args.checkpoint:
        if not args.checkpoint_sha256:
            parser.error('restart boundary requires its explicit digest')
        checkpoint_binding = explicit_binding(args.checkpoint, args.checkpoint_sha256)
        state = DiskMean(checkpoint_binding['path'], expected_identity=model_identity)
    else:
        state = DiskMean.create(args.output/'initial', initial_sectors(manifest, saved, grid, pairs),
            identity=model_identity, time=0., u=u, p=p,
            evidence=dict(model=model_spec, original_raw_frames_unchanged=True,
                          state='fixed original complex D; all raw columns and rows'))
    publish(args.output/'model.json', dict(identity=model_identity, specification=model_spec,
                                         initial_boundary=file_binding(state.path)))
    if args.mode == 'initialize':
        return
    background, background_receipt = restore_accepted(original)
    if not np.array_equal(arrays['saved_density_bohr3'], background.rho_updown):
        raise ValueError('aggregate expansion density differs from actual saved H0/K background')
    real = lambda value, name: real_physical(value, name, args.reality_tolerance)
    curvature = real(arrays['common_q_grid_ry_bohr2'][0, 0, 0], 'common Gamma curvature').reshape(ndof, ndof)
    atomic_mass = np.asarray(arrays['atomic_masses_amu'])
    masses = ionic_masses_ry(atomic_mass, np.arange(1, len(atomic_mass)+1), AMU_RY)
    reference = PreparedReference(real(arrays['saved_density_bohr3'], 'saved density'),
        real(arrays['ionic_force_ry_bohr'], 'ionic baseline').reshape(ndof),
        np.zeros_like(arrays['saved_density_bohr3'], dtype=float), np.zeros(ndof), curvature, masses,
        support, background.apply_kernel, 64, aggregate_binding['sha256'],
        manifest['ionic_terminal']['sha256'], saved.coverage_identity,
        ions.get('support_authority', 'no supplied support: explicit zero'),
        ions['coordinate_identity'], manifest['preparation_terminal']['sha256'])
    by_identity = {slot_identity(s): s for s in ordered}

    @contextmanager
    def factory(record):
        slot = by_identity[record.identity]
        # This reconstructs the actual slot basis from original bytes. No
        # transported H0 assumption, diagonal saved labels or cached band H.
        frame = grid.load_frame(slot)
        if (frame.basis.identity != record.basis_identity or frame.gauge_identity != record.gauge_identity
                or slot.source_binding['sha256'] != record.source_identity
                or frame.coefficients.shape != record.shape):
            raise ValueError('checkpoint differs from original full raw-column basis/gauge')
        basis = frame.basis
        del frame
        coupling = BlochDyadCoupling(background.inputs, basis, block_size=args.block_size)
        h0 = MaterialHamiltonian(background, coupling.factors, slot.ispin, block_size=args.block_size)
        yield coupling, h0.apply

    measure = abs(float(np.linalg.det(grid.reference.cell_bohr)))/int(np.prod(grid.reference.fft_smooth))
    model = StreamedMeanEvolution(reference, factory, model_identity=model_identity,
        expected_slots=tuple(slot_identity(s) for s in ordered), grid_measure=measure,
        reality_tolerance=args.reality_tolerance, max_h0_actions=args.max_h0_actions)
    model.validate(state)
    if args.mode == 'derivative':
        initial, result = model.derivative(state, args.output/'initial-derivative')
    else:
        initial = model.reduce(state, energy=True)
        state, history = model.evolve(state, args.output/'evolution',
            physical_interval=args.physical_interval_ry_inverse, steps=args.steps,
            max_steps=args.max_steps, max_iterations=args.max_midpoint_iterations,
            tolerance=args.midpoint_tolerance, refinement_label=args.refinement_label)
        result = dict(scope='bounded full-coverage mean evolution; not coherence/equilibrium evidence',
                      final_boundary=file_binding(state.path), history=history)
    diagnostics = dict(initial_density_bohr3=initial.density,
        initial_positive_gradient_ry_bohr=initial.gradient, initial_force_ry_bohr=initial.force,
        initial_kernel_ry=initial.potential,
        density_difference_from_aggregate=initial.density-arrays['conditional_density_bohr3'],
        gradient_difference_from_geometric_aggregate=initial.gradient-arrays['positive_gradient_ry_bohr'])
    with (args.output/'initial-diagnostics.npz').open('xb') as stream:
        np.savez(stream, **diagnostics)
    for entry in sources_before:
        authenticate(entry)
    result.update(schema='P254-streamed-material-evolution-v1', model_identity=model_identity,
        mode=args.mode, initial_energy_ry=initial.energy,
        initial_number=[initial.number.real, initial.number.imag], work_counts=model.counts,
        background_receipt=background_receipt, inventory=inventory,
        diagnostics=file_binding(args.output/'initial-diagnostics.npz'),
        implementation_sources=sources_before, runtime_context=runtime_context(),
        limitations=['All initial raw columns/rows and original D retained; mean need not commute with H0',
            'Curvature retains declared geometric-transport approximation; direct current mean forces recomputed at every stage',
            'No force, Hermitian or ASR repair; non-roundoff imaginary physical fields reject',
            'Finite-time accuracy requires physical-interval time refinement and nonlinear/exponential accuracy studies',
            'Gamma mean ionic motion only; full finite-q fluctuation/coherence evolution not executed here'])
    publish(args.output/'result.json', result)


if __name__ == '__main__':
    main()
