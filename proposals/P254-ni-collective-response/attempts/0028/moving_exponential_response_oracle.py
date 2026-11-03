#!/usr/bin/env python3
"""Executable finite-system oracle for moving_exponential_response.py.

Run from the repository root (no Native resources):
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python \
    proposals/P254-ni-collective-response/attempts/0028/moving_exponential_response_oracle.py \
    --output /tmp/p254-moving-exponential-response.json

Requires the existing NumPy/SciPy/SymPy environment. Output is exclusively
created; omit --output for stdout only. No production artifacts are read or
rewritten. The original oracle supplies its unequal-rank, complex raw-D finite
fixture and packing convention, NOT its saved verdict or fundamental matrix.
Here the reference independently evolves ORIGINAL FULL DENSITIES, both mean
and coupled response, by DOP853; the method evolves independent rectangular
factor legs about a separately integrated moving frame mean. Dense matrices
exist only in this finite oracle. This tests a finite numerical method, not
Ni/H material convergence, Native feasibility, useful long-time dynamics,
covariance pullbacks, damping, or speedup.
"""
import argparse
import copy
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import sympy as sp

from substrate_framework.joint_bloch_response import JointStage, JointTangent
from moving_exponential_response import (
    NAMES, JointH0Flow, JointResidual, SourceLocalMidpoint,
    SplitJointStage, refinement_difference)
from run_joint_response_oracle import FiniteSystem, DenseTangent


CHECKS = []


def require(name, predicate, **evidence):
    if not bool(predicate):
        raise AssertionError(name+': '+repr(evidence))
    CHECKS.append(dict(name=name, **evidence))


def max_error(a, b):
    return float(np.max(np.abs(np.asarray(a)-np.asarray(b)), initial=0))


def symbolic_blocks():
    # Original-density identity: no Hermiticity of the complex drive V is used.
    ns, nt, rs, rt = 3, 2, 2, 1
    hs, ht = sp.MatrixSymbol('Hs', ns, ns), sp.MatrixSymbol('Ht', nt, nt)
    ps, pt = sp.MatrixSymbol('Ps', ns, rs), sp.MatrixSymbol('Pt', nt, rt)
    ds, dt = sp.MatrixSymbol('Ds', rs, rs), sp.MatrixSymbol('Dt', rt, rt)
    left, right = sp.MatrixSymbol('L', nt, rs), sp.MatrixSymbol('R', rt, ns)
    drive = sp.MatrixSymbol('V', nt, ns)
    x = left*sp.Adjoint(ps)+pt*right
    ld = -sp.I*(ht*left+drive*ps*ds)
    rd = sp.I*(right*hs+dt*sp.Adjoint(pt)*drive)
    actual = ld*sp.Adjoint(ps)+sp.I*left*sp.Adjoint(ps)*hs-sp.I*ht*pt*right+pt*rd
    expected = -sp.I*(ht*x-x*hs+drive*ps*ds*sp.Adjoint(ps)-pt*dt*sp.Adjoint(pt)*drive)
    residual = sp.expand(actual-expected)
    require('exact_rectangular_moving_connection_identity', residual == sp.ZeroMatrix(nt, ns),
            source_pw=ns, target_pw=nt, source_rank=rs, target_rank=rt,
            drive_hermiticity_assumed=False, residual=str(residual))

    # Noncommuting series: A0, B(t), and dB/dt need not commute.
    a, b, c = sp.symbols('A0 B Bdot', commutative=False)

    def multiply(p, q):
        return [sp.expand(sum(p[j]*q[k-j] for j in range(k+1))) for k in range(3)]

    e, bm, kick = [1, a/2, a*a/8], [b, c/2, 0], [1, b/2, 0]
    base = multiply(e, e)
    feedback = multiply(multiply(multiply(e, bm), e), kick)
    coefficients = [base[0], base[1]+feedback[0], base[2]+feedback[1]]
    target = [1, a+b, ((a+b)*(a+b)+c)/2]
    defects = [sp.expand(x-y) for x, y in zip(coefficients, target)]
    require('exact_interaction_picture_midpoint_order_two', all(d == 0 for d in defects),
            coefficients=[str(sp.expand(x)) for x in coefficients],
            defects=[str(x) for x in defects])


class OriginalDensityMean:
    """Independent nonlinear mean: rho_dot=-i[H(rho,u,a),rho], no frames."""
    def __init__(self, system, end):
        self.s = system
        self.sizes = [n*n for n in system.sizes]
        rho0 = [p @ d @ p.conj().T for p, d in zip(system.frames, system.ds)]
        _, u, p, a, pi = system.unpack_mean(system.mean0)
        initial = self.pack(rho0, u, p, a, pi)
        self.solution = solve_ivp(self.rhs, (0, end), initial, method='DOP853',
                                  rtol=3e-13, atol=3e-14, dense_output=True)
        require('original_density_mean_completed', self.solution.success,
                message=self.solution.message, rhs_evaluations=self.solution.nfev)

    def pack(self, rho, u, p, a, pi):
        return np.concatenate([*(x.ravel() for x in rho), u, p, a, pi])

    def unpack(self, vector):
        rho, start = [], 0
        for n, size in zip(self.s.sizes, self.sizes):
            rho.append(vector[start:start+size].reshape(n, n))
            start += size
        return rho, vector[start:start+3], vector[start+3:start+6], vector[start+6:start+7], vector[start+7:start+8]

    def values(self, vector):
        s = self.s
        rho, u, p, a, pi = self.unpack(vector)
        # The mean is physical/real. No real projection is used for tangent X.
        density = s.w*sum(np.trace(v @ x) for v, x in zip(s.v, rho)).real
        hs = [h+s.k*(density-.04)*v+np.einsum('a,aij->ij', u.real, g)
              +np.einsum('a,aij->ij', (s.e.T @ a).real, ell)
              for h, v, g, ell in zip(s.h0, s.v, s.g, s.ell)]
        return rho, u, p, a, pi, hs

    def rhs(self, time, vector):
        s = self.s
        rho, u, p, a, pi, hs = self.values(vector)
        drho = [-1j*(h @ x-x @ h) for h, x in zip(hs, rho)]
        gradient = s.w*sum(np.einsum('aij,ji->a', g, x) for g, x in zip(s.g, rho))
        ell = s.w*sum(np.einsum('aij,ji->a', l, x) for l, x in zip(s.ell, rho))
        velocity = (p-s.charge @ a)/s.mass
        force = -(s.c+s.dc) @ u-gradient-s.mixed @ s.e.T @ a
        field_force = -s.n*(s.e @ ell+s.e @ s.contact @ s.e.T @ a
                            +s.e @ s.mixed.T @ u-s.charge.T @ velocity)
        return self.pack(drho, velocity, force, pi/s.mf, field_force)


class OriginalDensityResponse:
    """Original full rectangular X equation about independently evolved rho."""
    def __init__(self, system, response, mean):
        self.s, self.r, self.mean = system, response, mean
        self.packing = DenseTangent(system, response)

    def rhs(self, time, vector):
        s, r = self.s, self.r
        blocks, u, p, a, pi = self.packing.unpack(vector)
        rho, _, _, _, _, hs = self.mean.values(self.mean.solution.sol(time))
        density = sum(s.w*np.vdot(e.coupling.v, x) for e, x in zip(r.edges, blocks))
        gradient = sum(s.w*np.einsum('aij,ij->a', e.coupling.g.conj(), x)
                       for e, x in zip(r.edges, blocks))
        derivatives = []
        for edge, x in zip(r.edges, blocks):
            drive = s.k*density*edge.coupling.v+np.einsum('a,aij->ij', u, edge.coupling.g)
            if r.gamma:
                drive += np.einsum('a,aij->ij', np.sqrt(s.n)*s.e.T @ a, s.ell[edge.source])
            derivatives.append(-1j*(hs[edge.target] @ x-x @ hs[edge.source]
                                    +drive @ rho[edge.source]-rho[edge.target] @ drive))
        curvature = s.c+s.dc
        if not r.gamma:
            curvature = curvature+np.array([[0, .08j, 0], [-.08j, 0, 0], [0, 0, 0]])
        velocity, force = p/s.mass, -curvature @ u-gradient
        fdot, field_force = np.zeros(0, complex), np.zeros(0, complex)
        if r.gamma:
            ell = s.w*sum(np.einsum('aij,ji->a', l, x) for l, x in zip(s.ell, blocks))
            velocity = (p-np.sqrt(s.n)*s.charge @ a)/s.mass
            force -= np.sqrt(s.n)*s.mixed @ s.e.T @ a
            fdot = pi/s.mf
            field_force = -np.sqrt(s.n)*(s.e @ ell+s.e @ s.mixed.T @ u-s.charge.T @ (p/s.mass))
            field_force -= s.n*(s.e @ s.contact @ s.e.T+s.charge.T @ (s.charge/s.mass[:, None])) @ a
        return self.packing.pack(derivatives, velocity, force, fdot, field_force)


def fixture(end):
    s = FiniteSystem()
    # A nonstationary physical H0 with a scale well separated from the feedback.
    # Preserve full spaces and the original fixture's nonorthogonal frames/D.
    s.h0 = [h+np.diag(np.linspace(-30., 45., n)) for h, n in zip(s.h0, s.sizes)]
    s.mean = solve_ivp(s.mean_rhs, (0, end), s.mean0, method='DOP853',
                       rtol=3e-13, atol=3e-14, dense_output=True)
    require('moving_frame_mean_completed', s.mean.success,
            message=s.mean.message, rhs_evaluations=s.mean.nfev)
    return s


def split_provider(system, *, freeze=False):
    def supply(time):
        sample = 0. if freeze else time
        frames, _, _, _, _, _, _ = system.mean_values(system.mean.sol(sample))
        _, u, _, a, _ = system.unpack_mean(system.mean.sol(sample))
        rho = [ps @ d @ ps.conj().T for ps, d in zip(frames, system.ds)]
        density = system.w*sum(np.trace(v @ x) for v, x in zip(system.v, rho)).real
        # Deliberately construct residual directly, never subtract H0 actions.
        residual = [system.k*(density-.04)*v+np.einsum('a,aij->ij', u.real, g)
                    +np.einsum('a,aij->ij', (system.e.T @ a).real, ell)
                    for v, g, ell in zip(system.v, system.g, system.ell)]
        hs = [h+b for h, b in zip(system.h0, residual)]
        stage = JointStage(tuple(frames), tuple((lambda x, h=h: h @ x) for h in hs),
                           tuple((lambda x, h=h: h.conj().T @ x) for h in hs),
                           float(time), 'independent-current-moving-finite-stage' if not freeze else 'deliberately-wrong-frozen-mean')
        return SplitJointStage(stage, tuple((lambda x, b=b: b @ x) for b in residual),
                               tuple((lambda x, b=b: b.conj().T @ x) for b in residual))
    return supply


def free_flow(system, response):
    # Only this finite oracle diagonalizes dense H0. Native requires full-PW
    # action exponentials with an independently assessed action tolerance.
    cache = {}

    def action(index, vectors, tau, adjoint=False):
        key = index, float(tau), adjoint
        if key not in cache:
            h = system.h0[index].conj().T if adjoint else system.h0[index]
            eigenvalues, modes = np.linalg.eigh(h)
            cache[key] = (modes*np.exp(-1j*tau*eigenvalues)) @ modes.conj().T
        return cache[key] @ vectors

    return JointH0Flow(response, lambda i, x, t: action(i, x, t),
                        lambda i, x, t: action(i, x, t, True))


def tangent_vector(tangent):
    return np.concatenate([*(np.asarray(x).ravel() for x in tangent.ket),
                           *(np.asarray(x).ravel() for x in tangent.bra),
                           *(getattr(tangent, name) for name in NAMES)])


def actual_split_rhs(system, response, model, time, tangent):
    context = model.context(time, tangent)
    ket, bra = [], []
    for i, edge in enumerate(response.edges):
        l, b = model.edge_rhs(context, i, (tangent.ket[i], tangent.bra[i]), tangent)
        ket.append(l-1j*system.h0[edge.target] @ tangent.ket[i])
        bra.append(b+1j*tangent.bra[i] @ system.h0[edge.source])
    return JointTangent(tuple(ket), tuple(bra), *(getattr(context.small_rhs, name) for name in NAMES))


def density_rate(system, response, mapper, time, tangent, derivative):
    frames, _, _, _, _, _, hs = system.mean_values(system.mean.sol(time))
    blocks = []
    for edge, l, b, dl, db in zip(response.edges, tangent.ket, tangent.bra, derivative.ket, derivative.bra):
        ps, pt = frames[edge.source], frames[edge.target]
        blocks.append(dl @ ps.conj().T+l @ (1j*ps.conj().T @ hs[edge.source])
                      -1j*hs[edge.target] @ pt @ b+pt @ db)
    return mapper.pack(blocks, *(getattr(derivative, name) for name in NAMES))


def local_exposing_blocks(system, response, mean, initial, model, label):
    time = .137
    stage = model.stage(time).current
    actual = actual_split_rhs(system, response, model, time, initial)
    original = response.rhs(stage, initial)
    error = max_error(tangent_vector(actual), tangent_vector(original))
    require(label+'_split_is_actual_joint_rhs', error < 3e-12, max_error=error)
    reference = OriginalDensityResponse(system, response, mean)
    physical = reference.packing.physical(stage, initial)
    expected = reference.rhs(time, physical)
    physical_error = max_error(density_rate(system, response, reference.packing, time, initial, actual), expected)
    require(label+'_moving_original_density_rhs', physical_error < 3e-9, max_error=physical_error)

    density = sum(system.w*np.vdot(e.coupling.v, x)
                  for e, x in zip(response.edges, reference.packing.unpack(physical)[0]))
    imaginary = float(abs(density.imag))
    require(label+'_complex_density_is_exposed', imaginary > 1e-5, imaginary_density=imaginary)
    bad_complex = copy.copy(response)
    bad_complex.kernel = lambda n, **kw: system.k*n.real
    wrong = actual_split_rhs(system, bad_complex,
              JointResidual(bad_complex, split_provider(system), model.free_pair), time, initial)
    separation = max_error(tangent_vector(wrong), tangent_vector(original))
    require(label+'_reject_real_density_projection', separation > 1e-5, derivative_separation=separation)

    bad_normalization = copy.copy(response)
    if response.gamma:
        bad_normalization.sqrt_n = 1.
        mutation = 'omit orthonormal sqrt(N) matter/global-field coupling'
    else:
        bad_normalization.sectors = tuple(replace(s, weight=1.) for s in response.sectors)
        mutation = 'omit original explicit-slot 1/N density/force weight'
    wrong = actual_split_rhs(system, bad_normalization,
              JointResidual(bad_normalization, split_provider(system), model.free_pair), time, initial)
    separation = max_error(tangent_vector(wrong), tangent_vector(original))
    require(label+'_reject_wrong_normalization', separation > 1e-5,
            mutation=mutation, derivative_separation=separation)

    bad_d = copy.copy(response)
    bad_d.sectors = tuple(replace(s, d=s.d.real) for s in response.sectors)
    wrong = actual_split_rhs(system, bad_d,
              JointResidual(bad_d, split_provider(system), model.free_pair), time, initial)
    separation = max_error(tangent_vector(wrong), tangent_vector(original))
    require(label+'_reject_real_raw_D', separation > 1e-5, derivative_separation=separation)


def free_stiff_block(system, response, initial, flow, label):
    tau = .8
    stage = system.stage(0.)
    mapper = DenseTangent(system, response)
    xs = mapper.unpack(mapper.physical(stage, initial))[0]
    errors = []
    for i, edge in enumerate(response.edges):
        l, b = flow(i, (initial.ket[i], initial.bra[i]), tau)
        ps = expm(-1j*tau*system.h0[edge.source]) @ stage.frames[edge.source]
        pt = expm(-1j*tau*system.h0[edge.target]) @ stage.frames[edge.target]
        expected = expm(-1j*tau*system.h0[edge.target]) @ xs[i] @ expm(1j*tau*system.h0[edge.source])
        errors.append(max_error(l @ ps.conj().T+pt @ b, expected))
    max_phase = max(float(np.max(np.abs(np.linalg.eigvalsh(h))))*tau for h in system.h0)
    z = -1j*max_phase
    rk4_modulus = float(abs(1+z+z*z/2+z**3/6+z**4/24))
    require(label+'_stiff_H0_exact_full_row_flow', max(errors) < 1e-11 and rk4_modulus > 10,
            density_flow_error=max(errors), largest_free_phase=max_phase,
            rk4_free_mode_amplification=rk4_modulus, physical_free_mode_amplification=1.)


def exercise(system, mean, gamma, end, steps):
    label = 'Gamma' if gamma else 'nonGamma'
    response = system.response(gamma)
    rng = np.random.default_rng(281004 if gamma else 281005)

    def complex_random(shape):
        return .2*(rng.normal(size=shape)+1j*rng.normal(size=shape))

    zero = response.zero_tangent()
    initial = response.tangent(tuple(complex_random(x.shape) for x in zero.ket),
              tuple(complex_random(x.shape) for x in zero.bra), complex_random((3,)),
              complex_random((3,)), complex_random((response.nfield,)), complex_random((response.nfield,)))
    flow = free_flow(system, response)
    model = JointResidual(response, split_provider(system), flow)
    local_exposing_blocks(system, response, mean, initial, model, label)
    free_stiff_block(system, response, initial, flow, label)
    reference = OriginalDensityResponse(system, response, mean)
    y0 = reference.packing.physical(system.stage(0.), initial)
    loose = solve_ivp(reference.rhs, (0, end), y0, method='DOP853', rtol=3e-11, atol=3e-12)
    tight = solve_ivp(reference.rhs, (0, end), y0, method='DOP853', rtol=3e-13, atol=3e-14)
    require(label+'_independent_density_response_completed', loose.success and tight.success,
            loose_rhs_evaluations=loose.nfev, tight_rhs_evaluations=tight.nfev)
    expected = tight.y[:, -1]
    reference_change = max_error(loose.y[:, -1], expected)
    require(label+'_independent_reference_tolerance_refinement', reference_change < 3e-9,
            max_change=reference_change)
    scale = max(1., float(np.linalg.norm(expected)))
    nonhermitian = max(float(np.linalg.norm(x-x.conj().T)) for x in reference.packing.unpack(y0)[0]) if gamma else None
    if gamma:
        require(label+'_nonhermitian_complex_source', nonhermitian > .01, norm=nonhermitian)
    else:
        require(label+'_no_nonGamma_global_photon', response.nfield == 0 and len(initial.coordinate) == 0,
                nfield=response.nfield, ket_shapes=[list(x.shape) for x in initial.ket],
                bra_shapes=[list(x.shape) for x in initial.bra])

    def observable(solution):
        return reference.packing.physical(system.stage(solution.time), solution.state)

    rungs, solutions, values = [], [], []
    for count in steps:
        started = monotonic()
        solution = SourceLocalMidpoint(model).evolve(initial, np.linspace(0, end, count+1))
        value = observable(solution)
        error = float(np.linalg.norm(value-expected))/scale
        rungs.append(dict(steps=count, dt=end/count, relative_scaled_density_canonical_error=error,
                          context_bytes=solution.context_bytes, operation_counts=solution.operation_counts,
                          finite_oracle_elapsed_seconds=monotonic()-started))
        solutions.append(solution)
        values.append(value)
    ratios = [rungs[i]['relative_scaled_density_canonical_error']/rungs[i+1]['relative_scaled_density_canonical_error']
              for i in range(len(rungs)-1)]
    orders = [float(np.log2(ratio)) for ratio in ratios]
    require(label+'_moving_temporal_refinement', all(np.isfinite(ratios)) and
            all(x > 0 for x in ratios) and 2.5 < ratios[-1] < 5.8 and
            rungs[-1]['relative_scaled_density_canonical_error'] < 3e-4,
            rungs=rungs, error_ratios=ratios, observed_orders=orders)
    differences = [float(np.linalg.norm(values[i+1]-values[i]))/scale for i in range(len(values)-1)]
    require(label+'_self_refinement_agrees_with_reference_order',
            2.5 < differences[-2]/differences[-1] < 5.8,
            observable_differences=differences, last_difference_ratio=differences[-2]/differences[-1])
    estimate = refinement_difference(solutions[-2], solutions[-1], observable, np.linalg.norm)
    estimate['scope'] = 'empirical second-order estimate in finite packed density/canonical norm; not a physical-unit error bound'

    frozen_provider = split_provider(system, freeze=True)
    frozen = SourceLocalMidpoint(JointResidual(response, frozen_provider, flow)).evolve(
                initial, np.linspace(0, end, steps[-2]+1))
    wrong = reference.packing.physical(frozen_provider(end).current, frozen.state)
    frozen_error = float(np.linalg.norm(wrong-expected))/scale
    require(label+'_reject_wrong_frozen_mean', frozen_error > 20*rungs[-1]['relative_scaled_density_canonical_error']
            and frozen_error > 1e-3, frozen_relative_error=frozen_error,
            moving_fine_relative_error=rungs[-1]['relative_scaled_density_canonical_error'])
    return dict(q_fractional=response.q_fractional.tolist(), nfield=response.nfield,
                rungs=rungs, observed_orders=orders, refinement_estimate=estimate,
                wrong_frozen_mean_relative_error=frozen_error)


def binding(path):
    data = path.read_bytes()
    return dict(path=str(path.relative_to(ROOT)), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--steps', type=int, nargs=4, default=[12, 24, 48, 96],
                        help='four doubling refinement levels; default 12 24 48 96')
    args = parser.parse_args()
    if args.output is not None and args.output.exists():
        raise FileExistsError('preserve existing result: '+str(args.output))
    if args.steps[0] < 1 or any(b != 2*a for a, b in zip(args.steps[:-1], args.steps[1:])):
        parser.error('--steps must be four positive doubling levels')
    started, end = monotonic(), .27
    symbolic_blocks()
    system = fixture(end)
    mean = OriginalDensityMean(system, end)
    defects = []
    for time in (0., .137, end):
        frame_rho = system.mean_values(system.mean.sol(time))[5]
        density_rho = mean.values(mean.solution.sol(time))[0]
        defects.append(max(max_error(a, b) for a, b in zip(frame_rho, density_rho)))
    require('independent_frame_density_mean_agreement', max(defects) < 3e-9, max_errors=defects)
    moving = float(np.linalg.norm(system.mean.sol(end)-system.mean0))
    require('genuinely_moving_mean_not_stationary', moving > .1, mean_state_change=moving)
    require('original_complex_singular_raw_D_and_unequal_ranks',
            np.max(np.abs(system.ds[0].imag)) > .01 and abs(np.linalg.det(system.ds[1])) < 1e-12
            and len(set(system.ranks)) > 1, ranks=list(system.ranks), full_pw_dimensions=list(system.sizes),
            raw_D_imaginary_max=float(np.max(np.abs(system.ds[0].imag))),
            singular_raw_D_determinant_abs=float(abs(np.linalg.det(system.ds[1]))))
    results = [exercise(system, mean, gamma, end, args.steps) for gamma in (True, False)]
    result = dict(schema='P254-moving-exponential-response-oracle-v1', status='PASS',
        scope='Finite moving coupled response only; no Native allocation, material finite-time success, covariance adjoint or speedup claim',
        method='source-local H0 interaction-picture explicit midpoint, independently complex full-row ket/bra',
        source_bindings=[binding(HERE/'moving_exponential_response.py'), binding(Path(__file__).resolve()),
                         binding(HERE/'run_joint_response_oracle.py'),
                         binding(ROOT/'src/substrate_framework/joint_bloch_response.py')],
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                         sympy=sp.__version__, platform=platform.platform(),
                         OPENBLAS_NUM_THREADS=os.environ.get('OPENBLAS_NUM_THREADS'),
                         OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS')),
        elapsed_seconds=monotonic()-started, finite_process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        end_time_ry_inverse=end, results=results, checks=CHECKS)
    if args.output is None:
        print(json.dumps(result, indent=2, allow_nan=False))
    else:
        with args.output.open('x') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print(json.dumps(dict(status=result['status'], checks=len(CHECKS), output=str(args.output),
                              elapsed_seconds=result['elapsed_seconds'])))


if __name__ == '__main__':
    main()
