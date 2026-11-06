#!/usr/bin/env python3
"""Compare execution FFTs against exact DFT and actual full-PW material blocks.

This changes no live producer, physical support, source array or model. Reports
measured changed-path costs, not a forecast for the full propagated response.
"""
import argparse
from dataclasses import fields, replace
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import sympy as sp
from matched_joint_case import load_case
from substrate_framework.material_joint_response import RealBlochEdge
from substrate_framework.material_physical_operators import physical_hamiltonian, physical_finiteq
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.q_response_fft import ResponseRhoFFT
from substrate_framework.streamed_material_state import file_binding


def digest_array(value):
    return hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()


def exact_dft_check(check):
    miller = np.array(list(np.ndindex(2, 2, 2)), dtype=int).T
    fft = ResponseRhoFFT.create((2, 2, 2), miller, np.eye(3), 7.)
    coefficients = sp.Matrix([sp.Integer(i-3)+sp.I*(i % 3-1) for i in range(8)])
    points = [np.unravel_index(i, (2, 2, 2), order='F') for i in range(8)]
    transform = sp.Matrix([[(-1)**sum(int(miller[axis, j])*point[axis] for axis in range(3))
                            for j in range(8)] for point in points])
    check('exact_DFT_inverse_forward_and_pairing', transform.T*transform == 8*sp.eye(8))
    values = np.asarray(coefficients, dtype=complex).reshape(-1)
    reference = np.asarray(transform*coefficients, dtype=complex).reshape(-1)
    actual = fft.from_g(values)
    check('executed_F_order_inverse_matches_exact_DFT', np.array_equal(actual, reference))
    check('executed_forward_normalization_matches_exact_DFT', np.array_equal(fft.to_g(reference), values))
    check('missing_forward_norm_and_inverse_norm_mutations_rejected',
          not np.array_equal(fft.to_g(reference)*8, values)
          and not np.array_equal(actual/8, reference))
    expected_pairing = complex((coefficients.conjugate().T*coefficients)[0])
    check('independent_exact_Parseval_pairing',
          np.vdot(actual, actual)/8 == expected_pairing)


def execute(args, budget):
    checks, timings, errors = [], [], []
    def check(name, passed, **details):
        checks.append(dict(name=name, passed=bool(passed), **details))
        if not passed:
            raise AssertionError(name)
    exact_dft_check(check)
    context = load_case(args.case)
    source_slot = context.ordered[0]
    target_index = tuple((np.asarray(source_slot.mesh_index)+np.array([1, 0, 0])) % 4)
    target_slot = next(slot for slot in context.ordered
        if slot.ispin == source_slot.ispin and tuple(slot.mesh_index) == target_index)
    source, target = context.grid.load_frame(source_slot), context.grid.load_frame(target_slot)
    ref = context.control.reference
    volume = float(abs(np.linalg.det(ref.cell_bohr)))
    q = np.array([.25, 0., 0.])
    edge = RealBlochEdge.create(ref, source.basis, target.basis, q)
    source_fft = PhysicalRhoFFT.create(ref.fft_smooth, source.basis.miller, ref.reciprocal_bohr_inv, volume)
    target_fft = PhysicalRhoFFT.create(ref.fft_smooth, target.basis.miller, ref.reciprocal_bohr_inv, volume)
    density_fft = context.background.fft.at_q(ref.reciprocal_bohr_inv @ q)
    source_columns, target_columns = source.coefficients[:, :8], target.coefficients[:, :8]
    width = min(source_columns.shape[1], target_columns.shape[1])
    source_columns, target_columns = source_columns[:, :width], target_columns[:, :width]
    rho = context.background.rho_updown.astype(complex)
    rho += 1j*np.roll(rho.real, 1, axis=0)
    before = [digest_array(value) for value in (source_columns, target_columns, rho)]
    def compare(name, old, new, *inputs):
        budget.check()
        started = perf_counter(); expected = old(*inputs); old_seconds = perf_counter()-started
        started = perf_counter(); actual = new(*inputs); new_seconds = perf_counter()-started
        left, right = (expected, actual) if isinstance(expected, tuple) else ((expected,), (actual,))
        differences = [float(np.max(np.abs(a-b), initial=0.)) for a, b in zip(left, right)]
        check(name, len(left) == len(right) and all(a.shape == b.shape and a.dtype == b.dtype for a, b in zip(left, right))
              and all(value <= args.absolute_target for value in differences),
              max_abs_differences=differences, absolute_target=args.absolute_target)
        timings.append(dict(operation=name, original_seconds=old_seconds, execution_seconds=new_seconds))
        errors.append(dict(operation=name, max_abs_differences=differences))
        del expected, actual, left, right
        budget.check_processes(os.getpid(), [])
    for name, fft, columns in (('source', source_fft, source_columns), ('target', target_fft, target_columns)):
        fast = ResponseRhoFFT.from_physical(fft)
        check(name+'_all_physical_metadata_shared_exactly', all(getattr(fast, field.name) is getattr(fft, field.name) for field in fields(fft)))
        compare(name+'_full_PW_rows_inverse_columns', fft.from_g, fast.from_g, columns)
        wave = fft.from_g(columns)
        compare(name+'_forward_C_order', fft.to_g, fast.to_g, np.array(wave, order='C'))
        compare(name+'_forward_F_order', fft.to_g, fast.to_g, np.array(wave, order='F'))
        del wave
    fast_density = ResponseRhoFFT.from_physical(density_fft)
    compare('actual_density_nonzero_q_forward', density_fft.to_g, fast_density.to_g, rho)
    compare('actual_density_nonzero_q_gradient', density_fft.gradient, fast_density.gradient, rho)
    h0 = physical_hamiltonian(context.background, source.basis, source_slot.ispin, block_size=8)
    fast_h0 = replace(h0, factors=replace(h0.factors, fft=ResponseRhoFFT.from_physical(h0.factors.fft)))
    def original_h0(value):
        budget.action(value)
        return h0.apply(value)
    def execution_h0(value):
        budget.action(value)
        return fast_h0.apply(value)
    compare('actual_full_PW_H0_action', original_h0, execution_h0, source_columns)
    htarget = physical_hamiltonian(context.background, target.basis, target_slot.ispin, block_size=8)
    coupling = physical_finiteq(context.background.inputs, edge, h0.factors, htarget.factors, block_size=8)
    fast_coupling = physical_finiteq(context.background.inputs, edge, h0.factors, htarget.factors, block_size=8)
    bridge = fast_coupling.bridge
    for name in ('source_fft', 'target_fft', 'density_fft'):
        setattr(bridge, name, ResponseRhoFFT.from_physical(getattr(bridge, name)))
    compare('actual_rectangular_density_and_all_ionic_gradients',
        lambda a, b: coupling.moments(a, b, weight=source_slot.weight),
        lambda a, b: fast_coupling.moments(a, b, weight=source_slot.weight), target_columns, source_columns)
    direction = np.zeros((ref.nat, 3), complex); direction[0, 0] = 1.+.25j
    old_action, new_action = coupling.prepare(direction), fast_coupling.prepare(direction)
    compare('actual_complex_ionic_direction_forward', old_action.apply, new_action.apply, source_columns)
    compare('actual_independent_ionic_adjoint', old_action.adjoint_apply, new_action.adjoint_apply, target_columns)
    check('all_original_input_bytes_preserved', before == [digest_array(value) for value in (source_columns, target_columns, rho)])
    return dict(schema='P254-physical-response-FFT-check-v1', status='passed', checks=checks,
        case=args.case, physical_grid=list(ref.fft_smooth), q=q.tolist(), column_block=width,
        source_basis=source.basis.identity, target_basis=target.basis.identity,
        source_native=dict(source_slot.source_binding), target_native=dict(target_slot.source_binding),
        timings=timings, errors=errors, mathematical_model_unchanged=True,
        live_producer_changed=False, material_response=False, numerical_accuracy_certified=False,
        source_bindings=[file_binding(Path(__file__).resolve()),
            *[file_binding(ROOT/'src/substrate_framework'/name) for name in (
                'q_response_fft.py', 'qe_material_grid.py', 'qe_material_hamiltonian.py',
                'material_physical_operators.py', 'qe_bloch_finiteq.py', 'qe_bloch_finiteq_ionic.py',
                'material_joint_response.py', 'qe_saved_native.py')]],
        context=context.specification,
        scope='New FFT primitive/execution comparisons, not a replay of settled q/mean oracles. Exact small DFT, real source/target full-PW rows and actual density/H0/nonzero-q force/directional/adjoint blocks. Timings are these changed-path calls only, not a full-column or whole-q speedup claim. A genuine full-edge cost is still required before promoting the backend.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--absolute-target', type=float, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=16, max_seconds=args.max_seconds, max_bytes=1048576,
        max_memory_bytes=args.max_memory_bytes, worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(524288)
    started = perf_counter()
    result = dict(schema='P254-physical-response-FFT-check-v1', status='failed')
    try:
        result = execute(args, budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(seconds=perf_counter()-started, budget=budget.snapshot(),
        process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result['budget']['time_cap'] = 'Cooperative checks and external systemd hard wall cap; standalone comparison has no worker watchdog.'
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
