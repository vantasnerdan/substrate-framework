#!/usr/bin/env python3
"""Compare bounded full-PW probes; time only the NEW full346-column action.

Original old full-frame output was not retained and is deliberately not rerun.
No trajectory, new physical displacement, or whole-material timing is inferred.
"""
import argparse
import json
import os
from pathlib import Path
import resource
from time import perf_counter
from types import SimpleNamespace as NS

from run_fused_action_oracle import (
    ROOT, HERE, np, binding, source_bindings, publish, relative,
    FrozenMaterialAction, MaterialHamiltonian)
from run_streamed_material_evolution import (
    explicit_binding, bound_json, preload_parsers, restore_inputs,
    build_evolution_model, authenticate, DiskMean, physical_identity, implementation)
from substrate_framework.streamed_sector_executor import (
    physical_specification, KNOWN_T0_SHA256, KNOWN_T0_IDENTITY)
from substrate_framework.streamed_material_state import save_array
from scipy.linalg.blas import zgemm
from threadpoolctl import threadpool_limits, threadpool_info


def memory():
    fields = {}
    for path in ('/proc/self/status', '/proc/self/smaps_rollup'):
        for line in Path(path).read_text().splitlines():
            name, _, value = line.partition(':')
            if name in ('VmRSS', 'VmHWM', 'VmSize', 'Rss', 'Pss', 'Private_Dirty', 'Shared_Clean'):
                fields[name+'_kib'] = int(value.split()[0])
    fields['ru_maxrss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return fields


def probes(psi):
    """Two deterministic complement columns plus first/last original columns."""
    rng = np.random.default_rng(22528062)
    candidate = rng.normal(size=(psi.shape[0], 2))+1j*rng.normal(size=(psi.shape[0], 2))
    # Original raw columns are not assumed orthonormal. No QR/SVD truncation,
    # saved-label occupation selection or modified state enters this challenge.
    gram = zgemm(1., psi, psi, trans_a=2)
    overlaps = zgemm(1., psi, candidate, trans_a=2)
    complement = candidate-zgemm(1., psi, np.linalg.solve(gram, overlaps))
    norms = np.linalg.norm(complement, axis=0)
    if not np.isfinite(norms).all() or np.min(norms) <= 0:
        raise ValueError('deterministic full-PW complement does not exist')
    complement /= norms
    residual = float(np.max(np.abs(zgemm(1., psi, complement, trans_a=2))))
    if residual > 2e-11:
        raise ValueError('all-saved-column complement residual exceeds tolerance')
    return np.column_stack((psi[:, 0], psi[:, -1], complement)), dict(
        seed=22528062, saved_column_indices=[0, psi.shape[1]-1],
        complement_columns=[2, 3], complement_against_all_saved_columns=psi.shape[1],
        maximum_saved_overlap=residual,
        interpretation='full-PW complement of the entire saved span, not a fitted occupation sector')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('legacy-result', 'oracle'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fft-block-size', type=int, default=8)
    parser.add_argument('--projector-block-size', type=int, default=64)
    args = parser.parse_args()
    if not 1 <= args.fft_block_size <= 16 or args.projector_block_size < 1:
        parser.error('FFT block must be1..16 and projector block positive')
    os.sched_setaffinity(0, {args.cpu})
    sources = source_bindings()
    legacy_binding = explicit_binding(args.legacy_result, args.legacy_result_sha256)
    registration = json.loads((HERE/'fused-action-registration.json').read_text())
    if legacy_binding['sha256'] != registration['legacy_timing']['sha256']:
        raise ValueError('only the preregistered original full-frame action is accepted')
    old = bound_json(legacy_binding)
    if old['status'] != 'completed' or old['phase'] != 'action':
        raise ValueError('completed original action pilot required')
    oracle_binding = explicit_binding(args.oracle, args.oracle_sha256)
    oracle = bound_json(oracle_binding)
    if oracle.get('schema') != 'P254-fused-action-oracle-v1' or oracle.get('status') != 'passed':
        raise ValueError('passing sensitive fused dense oracle required')
    if oracle['source_bindings'] != sources:
        raise ValueError('oracle does not bind these exact old/new APIs and consumers')
    authenticate(oracle['arrays'])
    # Authenticate old API versions from the historical record, not just current
    # files. Other historical runtime modules are not silently relabeled current.
    historical_sources = {entry['path']: entry for entry in old['implementation_sources']}
    for path in registration['source_authority']:
        entry = historical_sources[str((ROOT/path).resolve())]
        authenticate(entry)
        if entry['sha256'] != registration['source_authority'][path]:
            raise ValueError('historical action API differs from preregistration')
    inputs = old['inputs']
    for entry in inputs.values():
        authenticate(entry)
    preload_parsers()
    arrays, receipt, manifest, saved, original, grid, ordered, pairs = restore_inputs(
        inputs['aggregate_receipt'], inputs['manifest'])
    ions = bound_json(inputs['initial_ions'])
    spec = dict(aggregate=inputs['aggregate_receipt'], manifest=inputs['manifest'],
        preparation=manifest['preparation_terminal'], initial_ions=inputs['initial_ions'],
        scalar='PreparedReference joint-second-order; common Gamma coefficient only; actual per-slot H0')
    state = DiskMean(inputs['checkpoint']['path'], verify=False)
    known = inputs['checkpoint']['sha256'] == KNOWN_T0_SHA256 and state.identity == KNOWN_T0_IDENTITY
    if not known and state.identity != physical_identity(spec):
        raise ValueError('boundary is not the original physical model')
    if (state.time != 0. or not state.evidence.get('original_raw_frames_unchanged')
            or physical_specification(state.evidence['model']) != spec
            or not np.array_equal(state.u, ions['u_bohr'])
            or not np.array_equal(state.p, ions['p_canonical'])):
        raise ValueError('unchanged original physical initial boundary required')
    selected = [slot for slot in state.slots if slot.identity == old['sector']['identity']]
    if len(selected) != 1:
        raise ValueError('historical sector absent from exact boundary')
    slot = selected[0]
    if (slot.frame != old['frame'] or slot.coefficient != old['original_D']
            or slot.shape[1] != 346 or slot.basis_identity != old['sector']['basis_identity']
            or slot.gauge_identity != old['sector']['gauge_identity']):
        raise ValueError('full346 raw frame/D or exact basis/gauge differs')
    authenticate(slot.frame)
    authenticate(slot.coefficient)
    settings = NS(block_size=args.fft_block_size, reality_tolerance=1e-10, max_h0_actions=None)
    model, background_receipt = build_evolution_model(settings, arrays, manifest, saved,
        original, grid, ordered, inputs['aggregate_receipt'], ions, state.identity)
    model.validate(state)
    stage = model.stage(arrays['conditional_density_bohr3'],
                        arrays['positive_gradient_ry_bohr'].reshape(-1), state.u)
    args.output.mkdir(parents=True, exist_ok=False)
    psi, raw_d = state.arrays(slot)
    probe, probe_receipt = probes(psi)
    probe_binding = save_array(args.output/'probe-columns.npy', probe)
    calls = []
    def measured(label, function, vectors):
        before, start = memory(), perf_counter()
        result = function(vectors)
        calls.append(dict(label=label, seconds=perf_counter()-start,
                          columns=vectors.shape[1], memory_before=before, memory_after=memory()))
        return result
    errors, outputs, directions = {}, {}, {}
    start = perf_counter()
    with model.factory(slot) as (coupling, legacy_h0):
        construction = dict(seconds=perf_counter()-start, memory=memory())
        h0 = getattr(legacy_h0, '__self__', None)
        if not isinstance(h0, MaterialHamiltonian) or h0.factors is not coupling.factors:
            raise ValueError('unchanged factory must expose its actual bound MaterialHamiltonian.apply')
        def frozen(field=None, u=None):
            return FrozenMaterialAction(h0, coupling, potential=field, direction=u,
                fft_block_size=args.fft_block_size, projector_block_size=args.projector_block_size)
        real_u = np.zeros(coupling.ndof)
        real_u[0], real_u[-1], real_u[3*(grid.reference.nat//2)+1] = .017, -.011, .007
        complex_u = real_u.astype(complex)
        complex_u[1], complex_u[-3] = .013j, -.009j
        x = np.arange(coupling.nr, dtype=float)
        challenge_field = stage.potential[:, slot.component]+.019*np.cos(.013*x)+.011j*np.sin(.017*x)
        cases = (
            ('h0', None, None),
            ('actual_initial_stage', stage.potential[:, slot.component], stage.u),
            ('nonzero_real_original_atoms', stage.potential[:, slot.component], real_u),
            ('complex_original_atoms_and_field', challenge_field, complex_u))
        for name, field, u in cases:
            current = frozen(field, u)
            def legacy(v):
                result = legacy_h0(v)
                if field is not None:
                    result += coupling.potential_action(field, v)
                if u is not None:
                    result += coupling.ionic_action(u, v)
                return result
            reference = measured(name+'_legacy_probe', legacy, probe)
            actual = measured(name+'_fused_probe', current.apply, probe)
            adjoint = measured(name+'_fused_adjoint_probe', current.adjoint_apply, probe)
            error = relative(actual, reference)
            bilinear_left = np.vdot(probe, reference[:, ::-1])
            bilinear_right = np.vdot(adjoint, probe[:, ::-1])
            errors[name] = dict(relative=error,
                maximum_absolute=float(np.max(np.abs(actual-reference))),
                adjoint_bilinear_relative=relative(bilinear_right, bilinear_left),
                legacy_norm=float(np.linalg.norm(reference)), fused_norm=float(np.linalg.norm(actual)))
            path = args.output/(name+'-probe-arrays.npz')
            np.savez(path, legacy=reference, fused=actual, difference=actual-reference, adjoint=adjoint,
                field=np.zeros(coupling.nr) if field is None else field,
                direction=np.zeros(coupling.ndof) if u is None else u,
                bilinear_left=bilinear_left, bilinear_right=bilinear_right)
            outputs[name] = binding(path)
            directions[name] = dict(array_file=outputs[name], field_key='field', direction_key='direction',
                semantics='actual original stage' if name == 'actual_initial_stage' else
                          'algebraic equivalence challenge only; not a physical initial condition',
                atom_order=list(grid.reference.atom_names), coordinate_order='original atom-major xyz')
        passed = max(max(v['relative'], v['adjoint_bilinear_relative']) for v in errors.values()) <= 2e-11
        publish(args.output/'probe-result.json', dict(status='passed' if passed else 'failed',
            errors=errors, arrays=outputs, probes=probe_binding, probe_receipt=probe_receipt,
            directions=directions, source_bindings=sources, inputs=inputs, original_D=slot.coefficient))
        if not passed:
            raise AssertionError('material probes failed; no full action launched; arrays retained')
        full = frozen(stage.potential[:, slot.component], stage.u)
        publish(args.output/'started-new-full-action.json', dict(sector=slot.metadata(), memory=memory(),
            cpu_affinity=sorted(os.sched_getaffinity(0)), source_bindings=sources))
        applied = measured('new_full346_initial_stage', full.apply, psi)
        new_norm = float(np.linalg.norm(applied))
        full_binding = save_array(args.output/'new-full346-action.npy', applied)
        # Preserve the original raw D ordering in the optional scalar receipt.
        projected = zgemm(1., psi, applied, trans_a=2)
        raw_trace = np.einsum('ij,ji->', raw_d, projected)
        del applied
    current_implementation = implementation((Path(__file__), HERE/'run_fused_action_oracle.py',
                                            HERE/'fused-action-registration.json'))
    if source_bindings() != sources:
        raise ValueError('source changed during material consumer')
    for entry in (slot.frame, slot.coefficient, *inputs.values()):
        authenticate(entry)
    new_seconds = calls[-1]['seconds']
    old_seconds = old['measurements']['full_frame_action']['seconds']
    publish(args.output/'result.json', dict(schema='P254-fused-action-material-v1', status='passed',
        scope='one exact full-PW sector action equivalence and new full346 timing; no evolution',
        source_bindings=sources, implementation_sources=current_implementation, inputs=inputs,
        oracle=oracle_binding, historical_result=legacy_binding, sector=slot.metadata(),
        frame=slot.frame, original_D=slot.coefficient, background_receipt=background_receipt,
        probe_receipt=probe_receipt, probe_columns=probe_binding, probe_errors=errors,
        probe_outputs=outputs, named_directions=directions, full_action=full_binding,
        construction=construction, measurements=calls, new_full_action_norm=new_norm,
        raw_D_trace=dict(real=float(raw_trace.real), imag=float(raw_trace.imag)),
        historical_full_action_norm=old['full_frame_action_norm'],
        historical_vs_new_time_ratio=old_seconds/new_seconds,
        timing_comparison_scope='historical single sample, not paired benchmark or equivalence proof',
        fft_block_size=args.fft_block_size, projector_block_size=args.projector_block_size,
        cpu_affinity=sorted(os.sched_getaffinity(0)), native_threadpools=threadpool_info(),
        h0_equivalent_work=dict(calls=len(calls), vector_columns=sum(c['columns'] for c in calls),
                               scope='explicit standalone consumer ledger; no production counter bypass integration'),
        limitations=['Original legacy full-frame output not retained and not rerun; elementwise comparison is on retained full-PW probes',
            'Full new output retained, including all346 original raw columns',
            'Probe-only synthetic nonzero directions/fields are not physical nuclear amplitudes',
            'No global128-sector reduction, nonlinear midpoint, trajectory or full rollout speedup',
            'Global H0 reservation cutover is a separate Root-owned integration',
            'Only selected frame/D authenticated; complete boundary and source inputs authenticated'] ))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
