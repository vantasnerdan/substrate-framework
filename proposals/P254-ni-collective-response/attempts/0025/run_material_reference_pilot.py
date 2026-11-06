#!/usr/bin/env python3
"""Actual saved-density Ni31H6 Gamma H0/K; no physical preparation or eigengate.

CLI: --admission JSON --output DIR [--context-only]. Common grant fields are
owned by material_run_contract.admit. Additional explicit schema:
inputs={xml:binding,density:binding,upfs:{Ni:binding,H:binding},
        wavefunctions:{'1':binding,'2':binding},probe_modes:binding,
        probe_coefficients:binding}; each binding is {path,bytes,sha256}.
settings={block_size:integer1..16,probe_amplitude_bohr_minus3:0.0001,
          probe_q_bohr_inv:[0,0,0]}. There is no implicit block-size default;
Root can explicitly grant 8. Probe binary layouts are little-endian int32
(3,3), complex128 (2,3), both F-order, original Native bound arrays.
library_binding binds the JSON consumed by NativePBE, not the ELF directly.

Six immutable NPZ/JSON pairs: background-coordinates, spin-1-diagnostics,
spin-2-diagnostics, kernel-probe-1/2/3. Fresh grants refuse existing pairs.
Resume needs identical immutable recipe and explicit remaining time, restores
completed arrays, and rebuilds charged radial/projector setup only as needed.
An interrupted whole-spin action loses its unfinished band blocks and derived
matrices; interrupted background/probe loses that whole uncommitted unit.
No automatic retries, interrupted-operation continuation or zero-work claim.
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


def validate_settings(grant):
    settings, inputs = grant['settings'], grant['inputs']
    if set(settings) != {'block_size', 'probe_amplitude_bohr_minus3', 'probe_q_bohr_inv'}:
        raise ValueError('settings require exactly block_size, probe_amplitude_bohr_minus3, probe_q_bohr_inv')
    if type(settings['block_size']) is not int or not 1 <= settings['block_size'] <= 16:
        raise ValueError('explicit block_size integer1..16 required; grant8 for default-size work')
    if type(settings['probe_amplitude_bohr_minus3']) not in (int, float) or settings['probe_amplitude_bohr_minus3'] != 1e-4:
        raise ValueError('original Native probe amplitude must be explicitly 1e-4 bohr^-3')
    q = settings['probe_q_bohr_inv']
    if not isinstance(q, list) or len(q) != 3 or any(type(x) not in (int, float) or x != 0 for x in q):
        raise ValueError('current material grant is explicit Gamma q=[0,0,0] only')
    if set(inputs) != {'xml', 'density', 'upfs', 'wavefunctions', 'probe_modes', 'probe_coefficients'}:
        raise ValueError('missing or extra material-reference input authority')
    if set(inputs['wavefunctions']) != {'1', '2'} or set(inputs['upfs']) != {'Ni', 'H'}:
        raise ValueError('both spin1/2 WFC and Ni/H UPF bindings required')
    if not isinstance(grant.get('library_binding'), dict):
        raise ValueError('bound backend library JSON required')


def probe_parameters(inputs):
    import numpy as np
    from substrate_framework.qe_saved_native import read_bound_bytes
    result = []
    for name, dtype, shape, count in (('probe_modes', '<i4', (3, 3), 36),
                                      ('probe_coefficients', '<c16', (2, 3), 96)):
        raw, source = read_bound_bytes(inputs[name], max_bytes=count)
        if len(raw) != count:
            raise ValueError('wrong original Native array byte count: '+name)
        array = np.frombuffer(raw, dtype=dtype).reshape(shape, order='F')
        if not np.isfinite(array).all():
            raise ValueError('nonfinite probe parameters')
        result.append((array, dict(source)))
    return result


def spin_diagnostics(background, factors, wave, spin, block_size):
    import numpy as np
    from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
    diagnostic = MaterialHamiltonian(background, factors, spin, block_size).saved_eigenpair_diagnostics(wave)
    point, _, occupations = background.inputs.reference.select(wave.ik, spin)
    arrays = {name: diagnostic[name] for name in ('action_ry', 'residual_ry', 'gram',
              'projected_h0_ry', 'residual_norm_ry', 'eigenvalues_ry')}
    arrays.update(raw_mv_occupations=occupations, eigenvalues_ha=point.eigenvalues_ha[spin-1],
                  pw_miller=wave.miller, reciprocal_bohr_inv=wave.reciprocal_bohr_inv,
                  k_bohr_inv=wave.k_cart_bohr_inv, kinetic_ry=factors.kinetic_ry,
                  source_bands_one_based=np.arange(1, wave.nbnd+1, dtype=np.int64),
                  raw_native_k_weight=np.asarray(point.weight))
    return arrays


def produce(run, grant, summary):
    import numpy as np
    from material_run_contract import bound_json, implementation, verify_consumers
    from material_reference_output import background_arrays, restore_background, raw_reference_metadata
    from substrate_framework.qe_internal_pbe import NativePBE
    from substrate_framework.qe_material_reference import load_material_inputs, reconstruct_reference
    from substrate_framework.qe_material_hamiltonian import build_pw_factors
    from substrate_framework.qe_saved_native import load_wavefunctions
    validate_settings(grant)
    if not grant['resume'] and any(run.target.glob('*.npz')):
        raise ValueError('fresh execution cannot use existing immutable blocks')
    run.check()
    start = time.monotonic()
    backend_binding = bound_json(grant['library_binding'])
    backend = NativePBE(backend_binding)
    summary['backend'] = dict(library_json=grant['library_binding'], bound_library=backend_binding,
                              loaded_source=dict(backend.source),
                              loaded_dependencies=[dict(item) for item in backend.dependencies])
    inputs = load_material_inputs(grant['inputs']['xml'], grant['inputs']['density'], grant['inputs']['upfs'])
    ref = inputs.reference
    if (ref.nat != 37 or ref.atom_names.count('Ni') != 31 or ref.atom_names.count('H') != 6
            or ref.nbnd_by_spin != (346, 346)):
        raise ValueError('requires actual Ni31H6 and both346saved bands')
    summary['raw_reference_metadata'] = raw_reference_metadata(ref)
    summary['input_sources'] = list(inputs.sources)
    summary['timing_seconds']['backend_input_radial_setup'] = time.monotonic()-start
    run.check()
    meta = dict(schema='material-reference-background-v1', physical_grid=list(ref.fft_grid),
                sources=list(inputs.sources), backend=summary['backend'],
                raw_reference_metadata=summary['raw_reference_metadata'],
                field_order='F physical grid; explicit up/down except rho_total_magnetization',
                units='density bohr^-3; potential Ry; gradient physical bohr derivatives')
    start = time.monotonic()
    restored = run.recover('background-coordinates', meta)
    if restored is None:
        background = reconstruct_reference(inputs, backend)
        run.commit('background-coordinates', background_arrays(background), meta)
    else:
        background = restore_background(inputs, restored)
        del restored
    summary['timing_seconds']['background_compute_or_restore_and_readback'] = time.monotonic()-start
    summary['background_sources'] = list(background.sources)
    factor_cache = None
    for spin in (1, 2):
        run.check()
        start = time.monotonic()
        key = f'spin-{spin}-diagnostics'
        meta = dict(schema='material-reference-spin-v1', spin_one_based=spin,
                    bands_one_based=[1, 346], source=grant['inputs']['wavefunctions'][str(spin)],
                    scope='Gamma full original PW rows; all346saved columns',
                    raw_mv='Raw saved XML occupations, width and k weight; never probabilities or temperature.',
                    qualification='Output-density H0 saved-eigenpair residuals, no convergence gate, shift or vector repair.',
                    raw_reference_metadata=summary['raw_reference_metadata'])
        restored = run.recover(key, meta)
        if restored is not None:
            # Arrays really are restored/read back; completed diagnostics never recompute.
            del restored
            summary['timing_seconds'][key] = time.monotonic()-start
            continue
        wave = load_wavefunctions(meta['source'])
        point, energies, _ = ref.select(wave.ik, spin)
        if (wave.ispin != spin or wave.ik != 1 or wave.nbnd != 346
                or wave.gamma_only or wave.npol != 1 or wave.scalef != 1
                or point.npw != wave.igwx or energies.size != 346
                or not np.array_equal(wave.k_cart_bohr_inv, np.zeros(3))
                or not np.allclose(point.k_cart_bohr_inv, wave.k_cart_bohr_inv, atol=2e-13, rtol=0)
                or not np.allclose(ref.reciprocal_bohr_inv, wave.reciprocal_bohr_inv, atol=2e-13, rtol=0)):
            raise ValueError('bound WFC must be matching full-complex Gamma spin/all346 bands')
        summary['input_sources'].append(dict(wave.source))
        setup = time.monotonic()
        shared = factor_cache is not None and factor_cache.matches(wave)
        if not shared:
            factor_cache = build_pw_factors(inputs, wave)
        summary.setdefault('projector_setup', []).append(dict(spin=spin, shared_exact_match=shared,
                    elapsed_seconds=time.monotonic()-setup, sources=list(factor_cache.sources)))
        # Raw actions and residuals dominate archive size; check before calculation.
        run.check(2*wave.igwx*346*16 + 2*346*346*16 + wave.miller.nbytes + 131072)
        arrays = spin_diagnostics(background, factor_cache, wave, spin, grant['settings']['block_size'])
        run.check()
        run.commit(key, arrays, meta)
        del arrays, wave
        summary['timing_seconds'][key] = time.monotonic()-start
    del factor_cache
    (modes, mode_source), (coefficients, coefficient_source) = probe_parameters(grant['inputs'])
    summary['input_sources'].extend((mode_source, coefficient_source))
    for probe in range(3):
        run.check()
        start = time.monotonic()
        key = f'kernel-probe-{probe+1}'
        meta = dict(schema='material-reference-kprobe-v1', probe_one_based=probe+1,
                    mode_source=mode_source, coefficient_source=coefficient_source,
                    amplitude_bohr_minus3=grant['settings']['probe_amplitude_bohr_minus3'],
                    q_bohr_inv=grant['settings']['probe_q_bohr_inv'],
                    formula='delta_n[r,s]=1e-4*coefficients[s,p]*exp(2*pi*i*mode[:,p].fractional_r)',
                    input_units='bohr^-3', output_units='Ry', spin_order=['up', 'down'],
                    qualification='Explicit mathematical K probe, not a physical preparation.')
        restored = run.recover(key, meta)
        if restored is not None:
            del restored
        else:
            found = np.flatnonzero(np.all(background.fft.miller == modes[:, probe, None], axis=0))
            if found.size != 1:
                raise ValueError('original probe G must occur exactly once on actual density grid')
            run.check(background.fft.nnr*2*16*2 + 65536)
            fourier = np.zeros((background.fft.miller.shape[1], 2), dtype=np.complex128)
            fourier[found[0]] = grant['settings']['probe_amplitude_bohr_minus3']*coefficients[:, probe]
            field = background.fft.from_g(fourier)
            action = background.apply_kernel(field, q_bohr_inv=grant['settings']['probe_q_bohr_inv'])
            run.commit(key, dict(delta_n_updown_bohr_minus3=field, kernel_action_ry=action,
                                 mode=modes[:, probe], coefficients=coefficients[:, probe]), meta)
            del fourier, field, action
        summary['timing_seconds'][key] = time.monotonic()-start
    run.check()
    sources = implementation((Path(__file__), Path(__file__).with_name('material_reference_output.py')))
    verify_consumers(grant, sources)
    summary['implementation_sources_after_lazy_consumers'] = sources
    summary['source_coverage_after_lazy_consumers'] = 'verified against admitted dependency universe and explicit bindings'
    summary['completion_mode'] = ('completed-cache-reuse-only' if all(x['disposition'] == 'completed-cache-reuse'
                                  for x in summary['outputs']) else 'new-completed-work-with-optional-cache-reuse')
    run.check()




def main():
    run = None
    target = None
    code = 1
    summary = dict(status='failed', generation='0025', run_id=uuid.uuid4().hex,
                   outputs=[], first_useful_block_seconds=None, timing_seconds={},
                   scope='Ni31H6 Gamma saved-density H0/K only; no full BvK, state or Native complement claim',
                   recovery_scope='Completed background/spin/probe arrays only. In-flight unit lost; setup charged again, no automatic retry.')
    try:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--admission', required=True, type=Path)
        parser.add_argument('--output', required=True, type=Path)
        parser.add_argument('--context-only', action='store_true')
        args = parser.parse_args()
        target = args.output.resolve()
        target.mkdir(parents=True, exist_ok=True)
        from material_run_contract import admit, implementation, runtime_context, binding
        from material_reference_output import OutputRun
        from material_publication import Publication
        run = OutputRun(target, summary, BOOT)
        def interrupted(signum, frame):
            raise InterruptedError(f'signal {signum}; only immutable completed checkpoints recover')
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
            signal.signal(signum, interrupted)
        summary['context'] = runtime_context()
        summary['implementation_sources'] = implementation((Path(__file__), Path(__file__).with_name('material_reference_output.py')))
        context_path = target/('context-'+summary['run_id']+'.json')
        context = dict(generation=summary['generation'], run_id=summary['run_id'],
                       context=summary['context'], implementation_sources=summary['implementation_sources'])
        grant, recipe, identity = admit(args.admission, 'material-reference',
                                      summary['implementation_sources'], summary['context'], args.context_only)
        run.publisher = Publication(target, summary, BOOT, grant, context_only=args.context_only)
        run.publisher.arm()
        if not args.context_only:
            summary.update(recipe=recipe, recipe_sha256=identity, admission=binding(args.admission),
                           remaining_seconds_granted=grant['remaining_seconds'], resume=grant['resume'])
            run.configure(grant, identity)
        run.publisher.write(context_path, context)
        if args.context_only:
            print(json.dumps(dict(event='material_reference_context', **context), allow_nan=False), flush=True)
            summary['status'] = 'context-only'
        else:
            run.event('material_reference_context', path=str(context_path), context=summary['context'])
            produce(run, grant, summary)
            summary['status'] = 'complete'
        code = 0
    except BaseException as error:
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
        if isinstance(error, SystemExit) and error.code == 0:
            code = 0
    finally:
        if run is not None:
            code = run.finish(code)
        elif code:
            # No admitted authority: journal only, never an unbounded file fallback.
            print(json.dumps({'event': 'pre_output_failure', 'run_id': summary['run_id'],
                              'exit_code': code, 'error': str(summary.get('failure'))[:240]}), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
