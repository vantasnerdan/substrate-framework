#!/usr/bin/env python3
"""Focused paired-ionic oracle: exact identities, finite probes, real material.

--columns full drives EVERY original saved source and target column, including
unequal original band counts. It does not mean a full-PW-row saved-band subset.
Only original PhysicalRhoFFT is used. Run under an external hard wall/RSS cap;
QBudget checks are cooperative between changed-path calls, not a watchdog.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0027'),
                str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import sympy as sp
from bloch_oracle_fixture import fixture
from run_finiteq_oracle import TwoProjectorNC, literal_beta, literal_local
from matched_joint_case import load_case
from substrate_framework.material_joint_response import RealBlochEdge
from substrate_framework.material_physical_operators import physical_hamiltonian, physical_finiteq
from substrate_framework.paired_finiteq_moments import paired_moments_and_drive
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.qe_bloch_edge import BlochEdge
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic
from substrate_framework.qe_bloch_grid import build_bloch_grid, UNITARY_POLICY
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.streamed_material_state import file_binding


class Checks:
    def __init__(self, target, budget):
        self.target, self.budget, self.records = target, budget, []

    def require(self, name, passed, **details):
        self.budget.check_processes(os.getpid(), [])
        self.records.append(dict(name=name, passed=bool(passed), **details))
        if not passed:
            raise AssertionError(name)

    def close(self, name, actual, expected):
        a, b = np.asarray(actual), np.asarray(expected)
        same_shape = a.shape == b.shape
        difference = float(np.max(np.abs(a-b), initial=0.)) if same_shape else None
        self.require(name, same_shape and np.isfinite(a).all() and np.isfinite(b).all()
                     and difference <= self.target,
                     max_abs_difference=difference, absolute_target=self.target)

    def mutation(self, name, mutant, expected):
        difference = float(np.max(np.abs(np.asarray(mutant)-expected), initial=0.))
        threshold = max(10*self.target, 1e-9)
        self.require(name, np.isfinite(difference) and difference > threshold,
                     mutant_max_abs_difference=difference, sensitivity_threshold=threshold)


def symbolic(check):
    # Folding/volume/weight belong to the dyad, not the local action.
    l, r, v = sp.symbols('l r v', complex=True)
    theta = sp.symbols('theta', real=True)
    omega, w = sp.symbols('omega w', positive=True)
    phase = sp.exp(-sp.I*theta)
    n = w*sp.conjugate(phase)*l*sp.conjugate(r)/omega
    check.require('sympy_fold_volume_weight_pairing',
                  sp.simplify(omega*sp.conjugate(n)*v
                              - w*sp.conjugate(l)*phase*v*r) == 0)

    # The independent adjoint conjugates the entire local multiplier.
    a, b = sp.symbols('a b', complex=True)
    check.require('sympy_local_independent_adjoint',
                  sp.simplify(sp.conjugate(a)*phase*v*b
                              - sp.conjugate(sp.conjugate(phase*v)*a)*b) == 0)

    # Local finite-q support differentiates G+q, never orbital k+G.
    reciprocal, q, position, displacement = sp.symbols('G q tau displacement', real=True)
    amplitude, local_direction = sp.symbols('amplitude local_direction', complex=True)
    local = amplitude*sp.exp(-sp.I*(reciprocal+q)*(position+displacement*local_direction))
    check.require('sympy_local_qplusG_derivative',
                  sp.simplify(sp.diff(local, displacement)
                              + sp.I*(reciprocal+q)*local_direction*local) == 0)

    # The reverse transpose has a positive i(G+q) and conjugated potential.
    coefficient = sp.symbols('density_coefficient', complex=True)
    coordinate_vertex = -sp.I*(reciprocal+q)*amplitude*sp.exp(-sp.I*(reciprocal+q)*position)
    positive_local = sp.I*(reciprocal+q)*coefficient*sp.conjugate(amplitude)*sp.exp(sp.I*(reciprocal+q)*position)
    check.require('sympy_local_positive_gradient_transpose',
                  sp.simplify(sp.conjugate(coordinate_vertex)*coefficient-positive_local) == 0)

    # Both projector placements differentiate, including a complex direction.
    eps, ps, pt, tau = sp.symbols('eps ps pt tau', real=True)
    u, beta = sp.symbols('u beta', complex=True)
    placed = beta*sp.exp(-sp.I*(pt-ps)*(tau+eps*u))
    check.require('sympy_both_projector_derivatives',
                  sp.simplify(sp.diff(placed, eps)
                              - sp.I*(ps-pt)*u*placed) == 0)

    # Positive gradient pairing is conjugate-linear in direction, not a force.
    u1, u2, g1, g2, x = sp.symbols('u1 u2 g1 g2 x', complex=True)
    gradient_pair = sp.conjugate(u1)*sp.conjugate(g1)*x + sp.conjugate(u2)*sp.conjugate(g2)*x
    check.require('sympy_positive_gradient_direction_pairing',
                  sp.simplify(gradient_pair-sp.conjugate(u1*g1+u2*g2)*x) == 0)

    # Independent ranks: a rank-one and a rank-two dyad sum before reduction.
    ell = sp.Matrix([[1+sp.I], [2-sp.I]])
    rr = sp.Matrix([[3-sp.I], [1+2*sp.I], [-2]])
    ell2 = sp.Matrix([[1, sp.I], [2-sp.I, -1]])
    rr2 = sp.Matrix([[2, -sp.I], [1+sp.I, 3], [-1, 2]])
    g = sp.Matrix([[1, sp.I, 2], [2-sp.I, -1, 3*sp.I]])
    dyad = ell*rr.conjugate().T + ell2*rr2.conjugate().T
    check.require('sympy_rectangular_unequal_rank_trace_additivity',
                  sp.simplify(sp.trace(g.conjugate().T*dyad)
                    - sp.trace(g.conjugate().T*ell*rr.conjugate().T)
                    - sp.trace(g.conjugate().T*ell2*rr2.conjugate().T)) == 0)


def digest_array(value):
    """Bounded scratch hashing, including all columns of mapped native frames."""
    array = np.asarray(value)
    digest = hashlib.sha256()
    digest.update(str((array.shape, array.dtype.str)).encode())
    if array.ndim == 2:
        for start in range(0, array.shape[1], 8):
            digest.update(np.ascontiguousarray(array[:, start:start+8]).view(np.uint8))
    else:
        digest.update(np.ascontiguousarray(array).view(np.uint8))
    return digest.hexdigest()


def original(coupling, left, right, left2, right2, direction, weight):
    n, g = coupling.moments(left, right, weight=weight)
    n2, g2 = coupling.moments(left2, right2, weight=weight)
    n += n2
    g += g2
    action = coupling.prepare(direction)
    return n, g, action.apply(right), action.adjoint_apply(left2)


def literal_vertex(inputs, edge, direction, *, lost=None, diagonal_d=False,
                   wrong_target_k=False):
    """Tiny dense analytic Gaussian model; independent of action/FFT code."""
    ref = inputs.reference
    ps = (ref.reciprocal_bohr_inv @ edge.source_basis.miller
          + edge.source_basis.k_cart_bohr_inv[:, None]).T
    kt = edge.source_basis.k_cart_bohr_inv if wrong_target_k else edge.target_basis.k_cart_bohr_inv
    pt = (ref.reciprocal_bohr_inv @ edge.target_basis.miller+kt[:, None]).T
    gq = (ref.reciprocal_bohr_inv @ inputs.density.miller+edge.q_cart_bohr_inv[:, None]).T
    local = np.zeros(len(gq), complex)
    nonlocal_part = np.zeros((len(pt), len(ps)), complex)
    for atom, name in enumerate(ref.atom_names):
        table, tau, u = inputs.tables[name], ref.positions_bohr[:, atom], direction[atom]
        local += (-1j*(gq @ u)*table.vloc(np.linalg.norm(gq, axis=1))
                  * np.exp(-1j*(gq @ tau)))
        bs, bt = literal_beta(ps, tau), literal_beta(pt, tau)
        d = table.upf.expanded_d
        if diagonal_d:
            d = np.diag(np.diag(d))
        base = bt @ d @ bs.conj().T
        if lost != 'source':
            nonlocal_part += 1j*base*(ps @ u)[None, :]
        if lost != 'target':
            nonlocal_part -= 1j*(pt @ u)[:, None]*base
    return literal_local(edge.source_basis.miller, edge.target_basis.miller,
                         edge.reciprocal_fold, inputs.density.miller, local)+nonlocal_part


def literal_moments(inputs, edge, left, right, left2, right2, weight):
    ref = inputs.reference
    x = np.indices(ref.fft_grid).reshape(3, -1, order='F').T/np.array(ref.fft_grid)
    us = np.exp(2j*np.pi*(x @ edge.source_basis.miller))
    ut = np.exp(2j*np.pi*(x @ edge.target_basis.miller))
    channel = (np.sum((ut @ left)*(us @ right).conj(), axis=1)
               + np.sum((ut @ left2)*(us @ right2).conj(), axis=1))
    channel *= weight*np.exp(2j*np.pi*(x @ edge.reciprocal_fold))/abs(np.linalg.det(ref.cell_bohr))
    density = np.zeros((len(x), 2), complex)
    density[:, edge.ispin-1] = channel
    dyad = left @ right.conj().T + left2 @ right2.conj().T
    gradient = []
    for atom in range(ref.nat):
        for axis in range(3):
            unit = np.zeros((ref.nat, 3)); unit[atom, axis] = 1.
            gradient.append(weight*np.vdot(literal_vertex(inputs, edge, unit), dyad))
    return density, np.asarray(gradient)


def finite_probes(check):
    rng = np.random.default_rng(2542829)
    def random(shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    with TemporaryDirectory(prefix='p254-paired-ionic-') as directory:
        inputs, bindings, _, _ = fixture(directory)
        omega = abs(np.linalg.det(inputs.reference.cell_bohr))
        inputs.tables['X'] = TwoProjectorNC(omega)
        # Exercise full complex Hermitian off-diagonal D, not only real D.
        inputs.tables['X'].upf.expanded_d = np.array([[1.3, .41+.27j], [.41-.27j, -.7]])
        grid = build_bloch_grid(inputs.reference, bindings, symmetry_policy=UNITARY_POLICY)
        cases = ((1, (0, 0, 1), (0, 0, 1)),
                 (2, (1, 0, 0), (1, 1, 0)), (1, (0, 0, 1), (0, 0, 0)))
        for spin, index, qindex in cases:
            target_index = tuple((np.array(index)+qindex) % 4)
            source = grid.load_frame(next(s for s in grid.slots if s.ispin == spin and s.mesh_index == index))
            target = grid.load_frame(next(s for s in grid.slots if s.ispin == spin and s.mesh_index == target_index))
            edge = BlochEdge.create(inputs.reference, source.basis, target.basis, np.array(qindex)/4)
            c = BlochFiniteQIonic(inputs, edge, block_size=2)
            ns, nt = edge.source_dimension, edge.target_dimension
            left, right = random((nt, 3)), random((ns, 3))
            left2, right2 = random((nt, 5)), random((ns, 5))
            direction, weight = random((c.nat, 3)), .137
            tag = f'finite_spin{spin}_q{qindex}'
            arrays = (left, right, left2, right2, direction)
            before = [digest_array(a) for a in arrays]
            actual = paired_moments_and_drive(c, *arrays, weight=weight)
            expected = original(c, *arrays, weight)
            vertex = literal_vertex(inputs, edge, direction)
            literal = (*literal_moments(inputs, edge, left, right, left2, right2, weight),
                       vertex @ right, vertex.conj().T @ left2)
            for name, a, b, d in zip(('density', 'gradient', 'forward', 'adjoint'), actual, expected, literal):
                check.close(tag+'_original_'+name, a, b)
                check.close(tag+'_literal_'+name, a, d)
            check.require(tag+'_input_bytes_preserved', before == [digest_array(a) for a in arrays])
            if spin == 1 and any(qindex):
                check.require(tag+'_rectangular_PW_rows_and_nonzero_fold', ns != nt and np.any(edge.reciprocal_fold))
            check.close(tag+'_complex_positive_gradient_pairing',
                        np.vdot(direction.ravel(), actual[1]),
                        weight*(np.vdot(vertex @ right, left)+np.vdot(vertex @ right2, left2)))
            check.mutation(tag+'_wrong_gradient_sign', -actual[1], literal[1])
            check.mutation(tag+'_missing_weight', actual[0]/weight, literal[0])
            check.mutation(tag+'_missing_volume', actual[0]*omega, literal[0])
            single = c.moments(left, right, weight=weight)
            check.mutation(tag+'_omitted_second_density', single[0], literal[0])
            check.mutation(tag+'_omitted_second_gradient', single[1], literal[1])
            wrong_adjoint = literal_vertex(inputs, edge, direction.conj()).conj().T @ left2
            check.mutation(tag+'_unconjugated_adjoint_direction', wrong_adjoint, literal[3])
            check.mutation(tag+'_conjugated_forward_direction',
                           literal_vertex(inputs, edge, direction.conj()) @ right, literal[2])
            for side in ('source', 'target'):
                check.mutation(tag+'_omitted_'+side+'_nonlocal_derivative',
                               literal_vertex(inputs, edge, direction, lost=side) @ right, literal[2])
            check.mutation(tag+'_diagonalized_species_D',
                           literal_vertex(inputs, edge, direction, diagonal_d=True) @ right, literal[2])
            missing_atom = direction.copy(); missing_atom[-1] = 0
            check.mutation(tag+'_omitted_last_atom', literal_vertex(inputs, edge, missing_atom) @ right, literal[2])
            x = np.indices(inputs.reference.fft_grid).reshape(3, -1, order='F').T/np.array(inputs.reference.fft_grid)
            us = np.exp(2j*np.pi*(x @ edge.source_basis.miller))
            ut = np.exp(2j*np.pi*(x @ edge.target_basis.miller))
            unconjugated = (np.sum((ut @ left)*(us @ right), axis=1)
                           + np.sum((ut @ left2)*(us @ right2), axis=1))
            unconjugated *= weight*np.exp(2j*np.pi*(x @ edge.reciprocal_fold))/omega
            check.mutation(tag+'_unconjugated_density_right', unconjugated, literal[0][:, spin-1])
            if np.any(edge.reciprocal_fold):
                check.mutation(tag+'_wrong_fold_sign', actual[0][:, spin-1]
                               * np.exp(-4j*np.pi*(x @ edge.reciprocal_fold)), literal[0][:, spin-1])
            if any(qindex):
                check.mutation(tag+'_wrong_target_k',
                               literal_vertex(inputs, edge, direction, wrong_target_k=True) @ right, literal[2])


def response_columns(frame, width):
    """Independent complex opposite-endpoint probe; never truncate driven legs."""
    rows = frame.shape[0]
    result = np.empty((rows, width), dtype=np.complex128)
    for start in range(0, width, 8):
        stop = min(start+8, width)
        indices = np.arange(start, stop)
        result[:, start:stop] = frame[:, indices % frame.shape[1]]
        result[:, start:stop] *= (0.63+0.29j)*np.exp(.071j*indices)[None, :]
    return result


def material(args, check, budget):
    context = load_case(args.case)
    source_slot = context.ordered[0]
    target_index = tuple((np.asarray(source_slot.mesh_index)+np.array([1, 0, 0])) % 4)
    target_slot = next(slot for slot in context.ordered
        if slot.ispin == source_slot.ispin and tuple(slot.mesh_index) == target_index)
    source, target = context.grid.load_frame(source_slot), context.grid.load_frame(target_slot)
    ref, q = context.control.reference, np.array([.25, 0., 0.])
    edge = RealBlochEdge.create(ref, source.basis, target.basis, q)
    source_count, target_count = source.coefficients.shape[1], target.coefficients.shape[1]
    if args.columns == 'full':
        right, left2 = source.coefficients, target.coefficients
    else:
        right = source.coefficients[:, :args.block_columns]
        left2 = target.coefficients[:, :args.block_columns]
    check.require('material_all_PW_rows', right.shape[0] == edge.source_dimension
                  and left2.shape[0] == edge.target_dimension)
    if args.columns == 'full':
        check.require('material_every_original_saved_source_and_target_column',
                      right.shape[1] == source_count and left2.shape[1] == target_count)
    left = response_columns(target.coefficients, right.shape[1])
    right2 = response_columns(source.coefficients, left2.shape[1])
    rng = np.random.default_rng(254280128)
    direction = .17*(rng.normal(size=(ref.nat, 3))+1j*rng.normal(size=(ref.nat, 3)))
    h0 = physical_hamiltonian(context.background, source.basis, source_slot.ispin, block_size=args.block_size)
    ht = physical_hamiltonian(context.background, target.basis, target_slot.ispin, block_size=args.block_size)
    c = physical_finiteq(context.background.inputs, edge, h0.factors, ht.factors, block_size=args.block_size)
    check.require('material_original_PhysicalRhoFFT_only', all(type(fft) is PhysicalRhoFFT
                  for fft in (c.bridge.source_fft, c.bridge.target_fft, c.bridge.density_fft)))
    arrays = (left, right, left2, right2, direction)
    frozen = (*arrays, source.coefficients, target.coefficients,
              *[s.beta for pair in c.species for s in pair],
              *[s.active_d for pair in c.species for s in pair])
    before = [digest_array(a) for a in frozen]
    budget.check_processes(os.getpid(), [])
    started = perf_counter()
    expected = original(c, *arrays, source_slot.weight)
    original_seconds = perf_counter()-started
    budget.check_processes(os.getpid(), [])
    started = perf_counter()
    actual = paired_moments_and_drive(c, *arrays, weight=source_slot.weight)
    paired_seconds = perf_counter()-started
    budget.check_processes(os.getpid(), [])
    for name, a, b in zip(('density', 'positive_gradient', 'forward_drive', 'independent_adjoint_drive'), actual, expected):
        check.close('material_'+name, a, b)
    check.require('material_original_frames_factors_and_probes_preserved',
                  before == [digest_array(a) for a in frozen])
    return dict(case=args.case, context=context.specification, q=q.tolist(),
        source_basis=source.basis.identity, target_basis=target.basis.identity,
        source_native=dict(source_slot.source_binding), target_native=dict(target_slot.source_binding),
        physical_grid=list(ref.fft_smooth), column_selection=args.columns,
        original_source_columns=source_count, original_target_columns=target_count,
        consumed_source_columns=right.shape[1], consumed_target_columns=left2.shape[1],
        first_pair_rank=right.shape[1], second_pair_rank=left2.shape[1],
        source_PW_rows=right.shape[0], target_PW_rows=left2.shape[0],
        atom_count=ref.nat, spin=edge.ispin, sector_weight=source_slot.weight,
        block_size=args.block_size, response_probe='complex modulations of cyclic opposite-endpoint saved columns; original driven frames unmodified',
        timing=dict(original_seconds=original_seconds, paired_seconds=paired_seconds,
            scope='Two original moments calls plus prepare/forward/adjoint versus paired call, including preparation; not a whole-edge or whole-response speedup.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--absolute-target', type=float, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--columns', choices=('block', 'full'), required=True)
    parser.add_argument('--block-columns', type=int, default=8)
    parser.add_argument('--block-size', type=int, choices=range(1, 17), default=8)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (not np.isfinite(args.absolute_target) or args.absolute_target <= 0
            or not np.isfinite(args.max_seconds) or args.max_seconds <= 0
            or args.max_memory_bytes <= 0 or args.block_columns <= 0):
        parser.error('positive finite targets/budgets and positive column counts required')
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1, max_seconds=args.max_seconds, max_bytes=1048576,
        max_memory_bytes=args.max_memory_bytes, worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(524288)
    check = Checks(args.absolute_target, budget)
    started = perf_counter()
    result = dict(schema='P254-paired-finiteq-moments-check-v1', status='failed')
    try:
        symbolic(check)
        finite_probes(check)
        result['material'] = material(args, check, budget)
        result['source_bindings'] = [file_binding(Path(__file__).resolve()),
            *[file_binding(ROOT/'src/substrate_framework'/name) for name in (
                'paired_finiteq_moments.py', 'qe_bloch_finiteq_ionic.py', 'qe_bloch_finiteq.py',
                'qe_material_grid.py', 'material_physical_operators.py', 'qe_material_hamiltonian.py')],
            file_binding(HERE/'run_finiteq_oracle.py'), file_binding(HERE.parent/'0027/bloch_oracle_fixture.py')]
        result['status'] = 'passed'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(checks=check.records, absolute_target=args.absolute_target, cpu=args.cpu,
        seconds=perf_counter()-started, budget=budget.snapshot(),
        process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Focused changed-path checks; original FFT/model/raw data and producer unchanged. No material-response completion or full-edge performance claim.')
    result['budget']['time_cap'] = 'Cooperative between calls; external hard wall/RSS cap required for native calls.'
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
