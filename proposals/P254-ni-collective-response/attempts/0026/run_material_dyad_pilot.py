#!/usr/bin/env python3
"""Actual Ni31H6 full-PW raw operator maps, not a material preparation.

CLI: --admission JSON --output DIR [--context-only]. Reuses0025 admission and
publication; material-dyad-maps has a fixed90s total cap. Grant registration
and source_access are the unchanged0025 control authority; settings ALSO binds
exact0026 continuation-registration.json. inputs: xml,density,upfs{Ni,H},
wavefunctions{1,2},reference_terminal (all exact path/bytes/sha256 bindings).
settings: block_size1..16, probe_recipe='saved-bands1-2-complex-all-coordinate-v1',
continuation_registration binding. library_binding must be absent/null.
resume_receipts={} fresh; otherwise spin-1-dyad/spin-2-dyad maps to immutable
receipt bindings. Explicit resume_identity and reduced cumulative time required.
No state D, physical temperature, force offset or held-sector authority inferred.
"""
import time
BOOT = time.monotonic()
import argparse
import json
from pathlib import Path
import signal
import sys
import traceback
import uuid

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))


def preload():
    """Import actual cold consumers BEFORE source inventory, no material reads."""
    from substrate_framework import material_dyad_output
    from substrate_framework.qe_material_dyads import GammaDyadCoupling
    from substrate_framework.qe_material_hamiltonian import build_pw_factors
    from substrate_framework.qe_material_ionic_action import GammaIonicAction
    material_dyad_output.preload_parsers()
    return (Path(__file__), Path(material_dyad_output.__file__))


def load_spin(inputs, background, spin):
    import numpy as np
    from substrate_framework.qe_saved_native import load_wavefunctions
    wave = load_wavefunctions(inputs['wavefunctions'][str(spin)])
    ref = background.inputs.reference
    point, _, _ = ref.select(1, spin)
    if (wave.ik != 1 or wave.ispin != spin or wave.igwx != 18613 or wave.nbnd != 346
            or wave.coefficients.shape != (18613, 346) or wave.gamma_only or wave.npol != 1
            or wave.scalef != 1 or point.npw != 18613
            or np.any(wave.k_cart_bohr_inv != 0) or np.any(point.k_cart_bohr_inv != 0)
            or not np.allclose(wave.reciprocal_bohr_inv, ref.reciprocal_bohr_inv,
                               rtol=0, atol=2e-13)):
        raise ValueError('requires exact original full-complex Gamma spin/346bands/18613PW')
    return wave


def channel_for(background, wave, block_size, cache):
    from substrate_framework.qe_material_dyads import GammaDyadCoupling
    from substrate_framework.qe_material_hamiltonian import build_pw_factors
    from substrate_framework.qe_material_ionic_action import GammaIonicAction
    if cache is not None and cache.factors.matches(wave):
        return cache, True
    inputs, ref = background.inputs, background.inputs.reference
    factors = build_pw_factors(inputs, wave)
    ionic = GammaIonicAction(ref.fft_grid, wave.miller, inputs.density.miller,
                            ref.reciprocal_bohr_inv, ref.positions_bohr,
                            ref.atom_names, inputs.tables)
    return GammaDyadCoupling(factors, ionic, block_size), False


def raw_probes(background, wave, spin):
    """Fixed mathematical probes: no populations, physical displacement or state."""
    import numpy as np
    # Copy only the two named full-PW columns, not all346 or any PW-square matrix.
    left = wave.coefficients[:, :2].copy()
    right = left.copy()
    coefficient = np.array([[.7+.2j, -.3+.4j], [.1-.6j, -.2+.5j]], dtype=np.complex128)
    shape = background.fft.physical
    field = background.vrs[:, spin-1].astype(np.complex128)
    grid = field.reshape(shape, order='F')
    grid += .25+.17j
    for axis, amplitude in enumerate((.31-.23j, .11+.19j, .07-.13j)):
        layout = [1, 1, 1]
        layout[axis] = shape[axis]
        phase = np.exp(2j*np.pi*np.arange(shape[axis])/shape[axis]).reshape(layout)
        grid += amplitude*phase
    index = np.arange(111, dtype=np.float64)
    direction = (1+index/111)+1j*(.25+index/222)
    return left, right, coefficient, field, direction


def compute_maps(background, wave, spin, channel, check):
    import numpy as np
    left, right, coefficient, field, direction = raw_probes(background, wave, spin)
    timings = {}
    start = time.monotonic()
    density, gradient = channel.moments(left, right, coefficient)
    timings['moments_full_grid_all111'] = time.monotonic()-start
    check()
    start = time.monotonic()
    potential_action = channel.potential_action(field, left)
    timings['complex_potential_action'] = time.monotonic()-start
    check()
    start = time.monotonic()
    ionic_action = channel.ionic_action(direction, left)
    timings['one_mixed_all111_ionic_action'] = time.monotonic()-start
    check()
    # Trace(C Rdagger A L); no conjugation of the physical field/direction.
    adjoint = right.conj().T
    pairs = np.array([
        [channel.grid_measure_bohr3*np.sum(density), np.trace(coefficient @ (adjoint @ left))],
        [channel.grid_measure_bohr3*np.dot(field, density), np.trace(coefficient @ (adjoint @ potential_action))],
        [np.dot(direction, gradient), np.trace(coefficient @ (adjoint @ ionic_action))]], dtype=np.complex128)
    ref = background.inputs.reference
    arrays = dict(left_pw=left, right_pw=right, coefficient=coefficient,
                  density_bohr_minus3=density, positive_gradient_ry_per_bohr=gradient,
                  potential_ry=field, potential_action_ry=potential_action,
                  direction_dimensionless=direction, ionic_action_ry_per_bohr=ionic_action,
                  comparison_pairs=pairs, pw_miller=wave.miller, positions_bohr=ref.positions_bohr,
                  cell_bohr=ref.cell_bohr, reciprocal_bohr_inv=ref.reciprocal_bohr_inv,
                  atom_one_based=np.asarray(ref.atom_indices), atom_symbols=np.asarray(ref.atom_names),
                  fft_shape=np.asarray(ref.fft_grid), dyad_bands_one_based=np.array([1, 2]))
    witness = {'imaginary_full_pw_trace': float(pairs[0, 1].imag),
               'nonhermitian_witness': bool(abs(pairs[0, 1].imag) > 1e-12),
               'reason': 'A Hermitian full-space operator has a real trace; no PW-square construction.'}
    return arrays, timings, witness


def reuse_block(run, key, metadata, receipt_binding):
    from material_run_contract import bound_json
    from substrate_framework.material_dyad_output import load_saved_block
    if Path(receipt_binding['path']).resolve() != run.target/(key+'.json'):
        raise ValueError('completed-block grant must name the actual immutable receipt')
    receipt = bound_json(receipt_binding)
    recorded = dict(receipt['metadata'])
    recorded.pop('comparisons')
    recorded.pop('raw_operator_diagnostic')
    if recorded != metadata:
        raise ValueError('completed-block provenance differs from current recipe')
    arrays, reports = load_saved_block(run.target, key, run.identity, receipt)
    witness = receipt['metadata']['raw_operator_diagnostic']
    del arrays
    run.record(key, receipt, 'completed-cache-reuse')
    return reports, witness


def produce_spin(run, grant, background, background_receipt, spin, channel, paths):
    from material_run_contract import implementation, verify_consumers
    from substrate_framework.material_dyad_output import block_metadata, validate_arrays
    key = f'spin-{spin}-dyad'
    metadata = block_metadata(grant['inputs'], grant['settings'], spin, background_receipt)
    started = time.monotonic()
    exists = any((run.target/(key+suffix)).exists() for suffix in ('.json', '.npz'))
    admitted = grant['resume_receipts'].get(key)
    if exists:
        if not grant['resume'] or admitted is None:
            raise ValueError('completed maps need explicit resume and receipt admission')
        reports, witness = reuse_block(run, key, metadata, admitted)
        timing = {'completed_map_readback': time.monotonic()-started, 'map_calls': 0}
    else:
        if admitted is not None:
            raise ValueError('admitted completed maps are missing')
        wave = load_spin(grant['inputs'], background, spin)
        run.check()
        setup_start = time.monotonic()
        channel, shared = channel_for(background, wave, grant['settings']['block_size'], channel)
        setup_seconds = time.monotonic()-setup_start
        run.check(16 << 20)
        arrays, timing, witness = compute_maps(background, wave, spin, channel, run.check)
        timing.update(input_and_setup_seconds=setup_start-started,
                      projector_channel_setup_seconds=setup_seconds,
                      exact_PW_factor_reuse=shared, map_calls=3)
        reports = validate_arrays(arrays)
        metadata.update(comparisons=reports, raw_operator_diagnostic=witness)
        verify_consumers(grant, implementation(paths))
        commit_start = time.monotonic()
        run.commit(key, arrays, metadata)
        timing['archive_commit_and_exact_readback'] = time.monotonic()-commit_start
        del arrays, wave
    timing['complete_spin_elapsed'] = time.monotonic()-started
    run.summary['timing_seconds'][key] = timing
    run.summary['comparisons'][key] = reports
    if not witness['nonhermitian_witness'] or not all(x['passed'] for x in reports.values()):
        raise ValueError('raw non-Hermitian map identity failed; arrays retained, no data repair')
    run.check()
    return channel


def produce(run, grant, paths):
    from material_run_contract import implementation, verify_consumers
    from substrate_framework.material_dyad_output import validate_grant, load_accepted_background
    validate_grant(grant)
    if not grant['resume'] and any(run.target.glob('*.npz')):
        raise ValueError('fresh grant cannot consume or overwrite completed maps')
    run.check()
    started = time.monotonic()
    background, receipt = load_accepted_background(grant['inputs'])
    run.summary['timing_seconds']['accepted_input_radial_setup_background_restore'] = time.monotonic()-started
    run.summary['accepted_background'] = receipt
    run.summary['background_consumption'] = 'Exact saved fields restored; vrs feeds complex probe, no PBE/SCF/PH or H0 producer'
    run.check()
    channel = None
    for spin in (1, 2):
        channel = produce_spin(run, grant, background, receipt, spin, channel, paths)
    sources = implementation(paths)
    verify_consumers(grant, sources)
    run.summary['final_implementation_sources'] = sources
    run.check()


def main():
    run, code = None, 1
    summary = dict(status='failed', generation='0026', run_id=uuid.uuid4().hex,
                   outputs=[], timing_seconds={}, comparisons={}, first_useful_block_seconds=None,
                   scope='Full37atom/111coordinate/18613PW/72cubed selected Gamma; both spins, not complete BvK',
                   state=None, D=None, physical_preparation=None, temperature=None,
                   recovery_scope='Completed whole-spin maps only; unfinished work lost; setup charged; no automatic retry')
    try:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--admission', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument('--context-only', action='store_true')
        args = parser.parse_args()
        def interrupted(signum, frame):
            raise InterruptedError(f'signal {signum}; only completed maps reusable')
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
            signal.signal(signum, interrupted)
        paths = preload()
        from material_run_contract import admit, implementation, runtime_context, binding
        from material_reference_output import OutputRun
        from material_publication import Publication
        target = args.output.resolve()
        target.mkdir(parents=True, exist_ok=True)
        run = OutputRun(target, summary, BOOT)
        context, sources = runtime_context(), implementation(paths)
        summary.update(context=context, implementation_sources=sources)
        grant, recipe, identity = admit(args.admission, 'material-dyad-maps', sources,
                                        context, args.context_only)
        run.publisher = Publication(target, summary, BOOT, grant, context_only=args.context_only)
        run.publisher.arm()
        if not args.context_only:
            summary.update(recipe=recipe, recipe_sha256=identity, admission=binding(args.admission),
                           remaining_seconds_granted=grant['remaining_seconds'], resume=grant['resume'])
            run.configure(grant, identity)
        context_path = target/('context-'+summary['run_id']+'.json')
        run.publisher.write(context_path, dict(context=context, implementation_sources=sources,
                                               generation='0026', run_id=summary['run_id']))
        run.event('material_dyad_context', context=binding(context_path))
        if args.context_only:
            summary['status'] = 'context-only'
        else:
            produce(run, grant, paths)
            summary['status'] = 'complete'
        code = 0
    except BaseException as error:
        summary['status'] = 'failed'
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
        if isinstance(error, SystemExit) and error.code == 0:
            code = 0
    finally:
        if run is not None:
            code = run.finish(code)
        elif code:
            print(json.dumps(dict(event='pre_output_failure', run_id=summary['run_id'],
                                  exit_code=code, error=str(summary.get('failure'))[:240])), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
