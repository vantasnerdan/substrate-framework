#!/usr/bin/env python3
"""Focused paired-field algebra, mutation and original Native material check.

Run only by the integration owner, with an external hard wall/memory allocation.
Timings cover the old Cartesian preparation/actions versus the new paired call;
they are not whole-edge, response-stage or propagation speedup measurements.
"""
import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import sympy as sp
from matched_joint_case import load_case
from substrate_framework.joint_bloch_response import CartesianFieldActions
from substrate_framework.material_physical_operators import PhysicalFieldVertices, physical_hamiltonian
from substrate_framework.paired_field_moments import paired_field_moments_and_drive
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.streamed_material_state import file_binding


def digest_array(value):
    return hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()


def exact_kinetic_check(check):
    p = sp.Symbol('p', real=True)
    l, r = sp.symbols('l r', complex=True)
    check('SymPy_kinetic_bilinear_orientation_and_factor_two',
          sp.expand(sp.conjugate(2*p*r)*l-2*p*sp.conjugate(r)*l) == 0)


def exact_projector_check(check):
    b = sp.Matrix([[1+sp.I, 2-sp.I], [-sp.I, 3+2*sp.I]])
    t = sp.Matrix([[2-sp.I, -1+sp.I], [1+2*sp.I, sp.I]])
    d0, d1, u, v = sp.symbols('d0 d1 u v', real=True)
    d = sp.Matrix([[d0, u+sp.I*v], [u-sp.I*v, d1]])
    l = sp.Matrix(sp.symbols('l0:2', complex=True))
    r = sp.Matrix(sp.symbols('r0:2', complex=True))
    action = t*d*b.H+b*d*t.H
    direct = ((action*r).H*l)[0]
    overlaps = ((t.H*r).H*(d*(b.H*l))+(d*(b.H*r)).H*(t.H*l))[0]
    check('SymPy_Hermitian_offdiagonal_D_direct_projector_bilinear',
          sp.expand(direct-overlaps) == 0)


def exact_complex_drive_check(check):
    b = sp.Matrix([[1+sp.I, 2], [2-sp.I, -sp.I]])
    jets = (sp.Matrix([[1, sp.I], [2-sp.I, 3]]),
            sp.Matrix([[sp.I, 2], [-1, 2+sp.I]]),
            sp.Matrix([[2, -sp.I], [sp.I, 1-sp.I]]))
    d = sp.Matrix([[2, 1+sp.I], [1-sp.I, 3]])
    a = sp.symbols('a0:3', complex=True)
    zero = sp.zeros(2)
    axes = [t*d*b.H+b*d*t.H for t in jets]
    directional = sum((a[j]*jets[j] for j in range(3)), zero)
    first_overlap = sum((a[j]*jets[j].H for j in range(3)), zero)
    fused = directional*d*b.H+b*d*first_overlap
    reference = sum((a[j]*axes[j] for j in range(3)), zero)
    check('SymPy_complex_forward_uses_unconjugated_second_overlap_coefficients',
          all(sp.expand(x) == 0 for x in fused-reference))
    adjoint = sum((sp.conjugate(a[j])*axes[j] for j in range(3)), zero)
    check('SymPy_independent_adjoint_conjugates_direction_only',
          all(sp.expand(x) == 0 for x in fused.H-adjoint))


def exact_pair_and_placement_check(check):
    p, tau = sp.symbols('p tau', real=True)
    phase = sp.exp(-sp.I*p*tau)
    check('SymPy_original_atom_phase_cancels_only_in_same_row_scalar_product',
          sp.simplify(sp.conjugate(phase)*phase-1) == 0)
    axis = sp.Matrix([[2, 1+sp.I], [1-sp.I, 3]])
    left = sp.Matrix([[1+sp.I], [2-sp.I]])
    right = sp.Matrix([[2-sp.I], [-1+2*sp.I]])
    left2 = sp.Matrix([[sp.I, 2-sp.I], [1+sp.I, -2]])
    right2 = sp.Matrix([[2+sp.I, -1], [sp.I, 1-2*sp.I]])
    expected = sp.trace((axis*right).H*left)+sp.trace((axis*right2).H*left2)
    compact = sp.trace(right.H*axis*left)+sp.trace(right2.H*axis*left2)
    reversed_second = sp.trace((axis*right).H*left)+sp.trace((axis*left2).H*right2)
    check('SymPy_two_unequal_width_pairs_exact_order', sp.expand(expected-compact) == 0)
    check('SymPy_reversing_second_pair_mutation_rejected', sp.expand(expected-reversed_second) != 0)


def original_response(vertex, left, right, left2, right2, direction, weight):
    axes = CartesianFieldActions.from_vertex(vertex)
    first = np.asarray([np.vdot(axis.apply(right), left) for axis in axes.axes])
    second = np.asarray([np.vdot(axis.apply(right2), left2) for axis in axes.axes])
    return ((weight*(first+second), axes.apply(direction, right),
             axes.apply(direction, left2, adjoint=True)), (weight*first, weight*second))


def row_fixture(h0, *, complex_d):
    """Small algebra fixture using real original radial tables and atom phases.

    Only this separate fixture optionally replaces D by a complex Hermitian
    offdiagonal probe. Actual material comparisons below never replace D.
    """
    f = h0.factors
    rows = np.unique(np.linspace(0, f.fft.miller.shape[1]-1, 11, dtype=int))
    fft = SimpleNamespace(miller=f.fft.miller[:, rows], reciprocal_bohr_inv=f.fft.reciprocal_bohr_inv)
    species = []
    changed = False
    for item in f.species:
        d = np.array(item.active_d, dtype=complex, copy=True)
        if complex_d and d.shape[0] >= 2:
            d[0, 1] += .7j
            d[1, 0] -= .7j
            changed = True
        species.append(replace(item, beta=np.asfortranarray(item.beta[rows]), active_d=d))
    if complex_d and not changed:
        raise ValueError('actual multichannel UPF required for complex Hermitian D fixture')
    factors = SimpleNamespace(fft=fft, k_bohr_inv=f.k_bohr_inv, species=tuple(species), sources=f.sources)
    return PhysicalFieldVertices(SimpleNamespace(factors=factors, background=h0.background),
                                 projector_block_size=2)


def numerical_fixture_checks(h0, check, compare, reject):
    vertex = row_fixture(h0, complex_d=True)
    rng = np.random.default_rng(281903)
    def random(width):
        shape = (len(vertex.momenta), width)
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    left, right, left2, right2 = random(2), random(2), random(3), random(3)
    direction = np.array([.31+.27j, -.18+.41j, .23-.35j])
    weight = .37
    inputs = (left, right, left2, right2)
    before = [digest_array(value) for value in inputs]
    expected, parts = original_response(vertex, *inputs, direction, weight)
    actual = paired_field_moments_and_drive(vertex, *inputs, direction, weight=weight)
    for name, a, b in zip(('ell', 'drive', 'adjoint_drive'), actual, expected):
        compare('small_complex_unequal_width_'+name, a, b)
    check('fixture_all_outputs_nonzero_complex', all(np.linalg.norm(x.real) > 0
          and np.linalg.norm(x.imag) > 0 for x in expected))
    reject('ell_conjugation', actual[0].conj(), expected[0])
    reject('ell_missing_weight', actual[0]/weight, expected[0])
    reject('ell_omitted_second_dyad', parts[0], expected[0])
    reject('ell_omitted_first_dyad', parts[1], expected[0])
    reject('ell_spurious_sqrtN', np.sqrt(64)*actual[0], expected[0])
    reject('ell_wrong_dual_sign', -actual[0], expected[0])
    reject('drive_spurious_weight', weight*actual[1], expected[1])
    axes = CartesianFieldActions.from_vertex(vertex)
    reversed_second = weight*np.asarray([np.vdot(axis.apply(right), left)
        + np.vdot(axis.apply(left2), right2) for axis in axes.axes])
    reject('second_pair_reversed', reversed_second, expected[0])
    reject('complex_forward_conjugated_direction', axes.apply(direction.conj(), right), expected[1])
    reject('complex_adjoint_treated_as_forward', axes.apply(direction, left2), expected[2])
    omitted = direction.copy(); omitted[2] = 0
    reject('complex_drive_omitted_z_axis', axes.apply(omitted, right), expected[1])
    matrices = tuple(axis.apply(np.eye(len(vertex.momenta), dtype=complex)) for axis in axes.axes)
    dense = sum(a*m for a, m in zip(direction, matrices))
    compare('independent_dense_forward', dense @ right, actual[1])
    compare('independent_dense_adjoint', dense.conj().T @ left2, actual[2])
    compare('independent_adjoint_pairing', np.vdot(left2[:, :2], actual[1]),
            np.vdot(actual[2][:, :2], right))
    kinetic = 2*(vertex.momenta @ direction)[:, None]*right
    reject('drive_omitted_nonlocal_projectors', kinetic, expected[1])
    reject('drive_omitted_kinetic', actual[1]-kinetic, expected[1])
    wrong_drive = kinetic.copy()
    unplaced_drive = kinetic.copy()
    transposed_d_drive = kinetic.copy()
    positions = vertex.hamiltonian.background.inputs.reference.positions_bohr
    for index, species in enumerate(vertex.hamiltonian.factors.species):
        nh = species.active_d.shape[0]
        db = sum(direction[j]*axes.axes[j].terms[index][1] for j in range(3))
        unplaced = np.empty_like(db)
        for slot, atom in enumerate(species.atom_indices):
            sl = slice(slot*nh, (slot+1)*nh)
            unplaced[:, sl] = db[:, sl]*np.exp(1j*(vertex.momenta @ positions[:, atom-1]))[:, None]
        for slot, _ in enumerate(species.atom_indices):
            sl = slice(slot*nh, (slot+1)*nh)
            b, first = species.beta[:, sl], db[:, sl]
            d = species.active_d
            weighted = d @ (b.conj().T @ right)
            wrong_drive += first @ weighted+b @ (d @ (first.conj().T @ right))
            # Both terms get the correct independent complex coefficients;
            # only placement is removed, isolating the original atom phase.
            first_unplaced = unplaced[:, sl]
            unplaced_second = sum(direction[j]*(axes.axes[j].terms[index][1][:, sl]
                *np.exp(1j*(vertex.momenta @ positions[:, species.atom_indices[slot]-1]))[:, None]).conj().T
                @ right for j in range(3))
            unplaced_drive += first_unplaced @ weighted+b @ (d @ unplaced_second)
            second = sum(direction[j]*axes.axes[j].terms[index][1][:, sl].conj().T @ right
                         for j in range(3))
            transposed_d_drive += first @ (d.T @ (b.conj().T @ right))+b @ (d.T @ second)
    reject('complex_directional_bank_wrong_second_conjugation', wrong_drive, expected[1])
    reject('original_atom_placement_omitted', unplaced_drive, expected[1])
    reject('complex_channel_D_transposed', transposed_d_drive, expected[1])
    diagonal_species = tuple(replace(s, active_d=np.diag(np.diag(s.active_d)))
                             for s in vertex.hamiltonian.factors.species)
    factors = SimpleNamespace(**vars(vertex.hamiltonian.factors)); factors.species = diagonal_species
    diagonal_vertex = PhysicalFieldVertices(SimpleNamespace(factors=factors, background=h0.background))
    diagonal, _ = original_response(diagonal_vertex, *inputs, direction, weight)
    reject('offdiagonal_channel_D_omitted', diagonal[1], expected[1])
    vector_inputs = tuple(value[:, 0] for value in inputs)
    vector_expected, _ = original_response(vertex, *vector_inputs, direction, weight)
    vector_actual = paired_field_moments_and_drive(vertex, *vector_inputs, direction, weight=weight)
    for name, a, b in zip(('ell', 'drive', 'adjoint_drive'), vector_actual, vector_expected):
        compare('single_vector_'+name, a, b)
    check('all_small_fixture_factor_bytes_preserved', before == [digest_array(value) for value in inputs])


def material_probes(columns):
    """Keep source columns exact; construct two independent nonzero complex legs."""
    row = np.arange(len(columns), dtype=float)[:, None]
    left = ((.43+.29j)+(.17-.11j)*np.cos(.013*row))*columns
    right2 = ((-.21+.37j)+(.09+.23j)*np.sin(.017*row))*columns
    # Distinct column permutations make these bilinears genuinely off-diagonal.
    left += (.19-.13j)*np.roll(columns, 1, axis=1)
    right2 += (.14+.08j)*np.roll(columns, -1, axis=1)
    return left, columns, columns, right2


def execute(args, budget):
    checks, timings = [], []
    def check(name, passed, **details):
        checks.append(dict(name=name, passed=bool(passed), **details))
        if not passed:
            raise AssertionError(name)
    def compare(name, actual, expected):
        a, b = np.asarray(actual), np.asarray(expected)
        error = float(np.max(np.abs(a-b), initial=0.))
        check(name, a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
              and error <= args.absolute_target, max_abs_difference=error,
              absolute_target=args.absolute_target)
    def reject(name, mutated, expected):
        error = float(np.max(np.abs(np.asarray(mutated)-np.asarray(expected)), initial=0.))
        check('mutation_rejected_'+name, np.isfinite(error) and error > args.absolute_target,
              max_abs_difference=error, absolute_target=args.absolute_target)
    exact_kinetic_check(check)
    exact_projector_check(check)
    exact_complex_drive_check(check)
    exact_pair_and_placement_check(check)
    budget.check_processes(os.getpid(), [])
    context = load_case(args.case)
    slot = context.ordered[args.slot_index]
    frame = context.grid.load_frame(slot)
    h0 = physical_hamiltonian(context.background, frame.basis, slot.ispin, block_size=8)
    check('material_keeps_original_PhysicalRhoFFT', type(h0.factors.fft) is PhysicalRhoFFT)
    check('material_retains_all_Native_PW_rows', len(frame.coefficients) == frame.basis.miller.shape[1])
    numerical_fixture_checks(h0, check, compare, reject)
    width = frame.coefficients.shape[1] if args.columns == 'full' else args.block_columns
    if width > frame.coefficients.shape[1]:
        raise ValueError('requested block exceeds original Native column count')
    columns = frame.coefficients[:, :width]
    check('full_mode_selects_every_original_column' if args.columns == 'full' else 'block_selects_explicit_original_columns',
          columns.shape == (frame.coefficients.shape[0], width)
          and (args.columns != 'full' or width == frame.coefficients.shape[1]))
    inputs = material_probes(columns)
    direction = np.array([.31+.27j, -.18+.41j, .23-.35j])
    vertex = PhysicalFieldVertices(h0)
    preserved = [*inputs, frame.coefficients, frame.basis.miller,
                 *[s.active_d for s in h0.factors.species], *[s.beta for s in h0.factors.species]]
    before = [digest_array(value) for value in preserved]
    for repeat in range(args.timing_repeats):
        budget.check_processes(os.getpid(), [])
        started = perf_counter()
        expected, parts = original_response(vertex, *inputs, direction, slot.weight)
        old_seconds = perf_counter()-started
        budget.check_processes(os.getpid(), [])
        started = perf_counter()
        actual = paired_field_moments_and_drive(vertex, *inputs, direction, weight=slot.weight)
        new_seconds = perf_counter()-started
        budget.check_processes(os.getpid(), [])
        for name, a, b in zip(('both_dyad_ell', 'complex_forward_drive', 'independent_complex_adjoint_drive'), actual, expected):
            compare('actual_material_'+name+'_repeat_'+str(repeat), a, b)
        reject('actual_material_omitted_first_pair_'+str(repeat), parts[1], expected[0])
        reject('actual_material_omitted_second_pair_'+str(repeat), parts[0], expected[0])
        reject('actual_material_second_pair_conjugated_'+str(repeat), parts[0]+parts[1].conj(), expected[0])
        check('actual_material_all_three_outputs_nonzero_'+str(repeat),
              all(float(np.linalg.norm(x)) > args.absolute_target for x in expected))
        compare('actual_material_independent_adjoint_pairing_'+str(repeat),
                np.vdot(inputs[2], actual[1]), np.vdot(actual[2], inputs[1]))
        timings.append(dict(repeat=repeat, original_cartesian_prepare_and_calls_seconds=old_seconds,
            paired_shared_jets_call_seconds=new_seconds, columns=width,
            PW_rows=len(columns), preparation_included=True,
            original_path='CartesianFieldActions.from_vertex; two ordered dyad reductions; independent forward and adjoint calls',
            changed_path='paired_field_moments_and_drive; one shared species Cartesian jet bank'))
        del expected, actual, parts
    check('all_original_material_factor_and_Native_input_bytes_preserved',
          before == [digest_array(value) for value in preserved])
    return dict(schema='P254-paired-field-moments-check-v1', status='passed', checks=checks,
        case=args.case, slot_index=args.slot_index, slot=dict(ispin=slot.ispin,
            mesh_index=list(slot.mesh_index), source_ik=slot.source_ik, weight=float(slot.weight)),
        original_native=dict(slot.source_binding), original_frame_gauge=frame.gauge_identity,
        basis_identity=frame.basis.identity, physical_grid=list(context.control.reference.fft_smooth),
        PW_rows=len(columns), selected_columns=width, original_columns=frame.coefficients.shape[1],
        column_selection=args.columns, full_original_columns=args.columns == 'full',
        complex_direction_real=direction.real.tolist(), complex_direction_imag=direction.imag.tolist(),
        absolute_target=args.absolute_target, timings=timings, context=context.specification,
        source_bindings=[file_binding(Path(__file__).resolve()), file_binding(HERE/'matched_joint_case.py'),
            *[file_binding(ROOT/'src/substrate_framework'/name) for name in (
                'paired_field_moments.py', 'material_physical_operators.py', 'material_mixed_field_ionic.py',
                'qe_material_field_derivatives.py', 'qe_nc_projectors.py', 'qe_material_grid.py',
                'qe_saved_native.py', 'qe_real_bloch_grid.py', 'joint_bloch_response.py')]],
        original_model_unchanged=True, original_FFT=True, production_integration_changed=False,
        material_response=False, numerical_accuracy_certified=False,
        scope='Small exact SymPy identities, nonzero complex unequal-width algebra fixture using original radial tables/positions with explicit independent complex-Hermitian-D probes, then genuine original-Native full-PW material columns with untouched physical channel D. Timings cover only the compared paired field path, not a whole edge, stage or propagated response.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--columns', choices=('block', 'full'), required=True)
    parser.add_argument('--block-columns', type=int, default=8)
    parser.add_argument('--slot-index', type=int, default=0)
    parser.add_argument('--timing-repeats', type=int, default=1)
    parser.add_argument('--absolute-target', type=float, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (args.block_columns < 1 or args.slot_index < 0 or args.timing_repeats < 1
            or not np.isfinite(args.absolute_target) or args.absolute_target <= 0):
        raise ValueError('positive block/repeat/absolute-target and nonnegative slot required')
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1, max_seconds=args.max_seconds, max_bytes=1048576,
        max_memory_bytes=args.max_memory_bytes, worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(524288)
    started = perf_counter()
    result = dict(schema='P254-paired-field-moments-check-v1', status='failed')
    try:
        result = execute(args, budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(seconds=perf_counter()-started, budget=budget.snapshot(),
        cpu=args.cpu, process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result['budget']['time_cap'] = 'Cooperative boundary checks; invoke under an external hard wall cap.'
    result['budget']['memory_accounting'] = 'Cooperative boundary RSS and process peak RSS; external memory allocation required for native-call hard enforcement.'
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
