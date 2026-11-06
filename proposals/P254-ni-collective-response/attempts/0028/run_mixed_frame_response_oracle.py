#!/usr/bin/env python3
"""Run finite mixed-frame mathematics, emitting actual structured evidence.

Usage: python run_mixed_frame_response_oracle.py --output /path/report.json
Requires NumPy/SciPy.  This is a finite-model mathematical oracle, NOT native
material propagation, an energy-transfer result, or a finite-amplitude bound.
No report verdict exists until this program is actually executed.
"""
import argparse
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import traceback

import numpy as np
from scipy.linalg import expm

import mixed_frame_response as production
from mixed_frame_oracle_finite_model import (
    analytic_density, dense_density_ode, embed_frames, evolve_production,
    exact_frame_jets, extract_frames, finite_model, materialize_descriptors,
    physical_actions,
)


class Evidence:
    """Record numerical defects and actual mutation sensitivities separately."""
    def __init__(self):
        self.checks = []
        self.mutations = []

    def close(self, name, actual, expected, *, atol=2e-9, rtol=2e-9):
        actual, expected = np.asarray(actual), np.asarray(expected)
        defect = float(np.max(np.abs(actual - expected), initial=0))
        scale = float(np.max(np.abs(expected), initial=0))
        threshold = atol + rtol * scale
        self.checks.append(dict(name=name, passed=bool(np.isfinite(defect) and defect <= threshold),
                                max_abs_defect=defect, threshold=threshold))

    def require(self, name, condition, **measurements):
        self.checks.append(dict(name=name, passed=bool(condition), **measurements))

    def reject(self, name, actual, expected, *, threshold=1e-7):
        defect = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected)), initial=0))
        self.mutations.append(dict(name=name, rejected=bool(np.isfinite(defect) and defect > threshold),
                                   max_abs_defect=defect, rejection_threshold=threshold))

    def rejects_input(self, name, function):
        try:
            function()
        except (ValueError, TypeError) as exc:
            self.mutations.append(dict(name=name, rejected=True, exception=type(exc).__name__, message=str(exc)))
        else:
            self.mutations.append(dict(name=name, rejected=False, exception=None))


def descriptors(model, stages):
    return production.density_dyads(*stages, model.D, target_plus=model.target_plus,
                                   target_minus=model.target_minus, source_weights=model.weights)


def check_model_contract(evidence, model, name):
    rho0 = model.initial_frame @ model.weighted_D @ model.initial_frame.conj().T
    commutator_norm = float(np.linalg.norm(model.h0 @ rho0 - rho0 @ model.h0))
    gram = model.initial_frame.conj().T @ model.initial_frame
    raw_gram_defect = float(np.linalg.norm(gram - np.eye(gram.shape[0])))
    evidence.close(name + '/Hminus_is_adjoint', model.hm, model.hp.conj().T, atol=1e-14)
    evidence.close(name + '/H0_Hermitian', model.h0, model.h0.conj().T, atol=1e-14)
    evidence.close(name + '/Hmixed_Hermitian', model.hx, model.hx.conj().T, atol=1e-14)
    ranks = [int(np.linalg.matrix_rank(d)) for d in model.D]
    evidence.require(name + '/singular_complex_D', all(rank == 1 for rank in ranks)
                     and max(float(np.max(np.abs(d.imag))) for d in model.D) > 1e-4,
                     numerical_D_ranks=ranks, original_column_ranks=list(model.ranks))
    evidence.require(name + '/nonstationary_raw_frame', commutator_norm > 1e-4 and raw_gram_defect > 1e-4,
                     commutator_norm=commutator_norm, raw_gram_defect=raw_gram_defect)


def check_rhs(evidence, model, stages, prefix):
    f0, fp, fm, fx = embed_frames(model, stages)
    dense = (-1j * (model.h0 @ fp + model.hp @ f0),
             -1j * (model.h0 @ fm + model.hm @ f0),
             -1j * (model.h0 @ fx + model.hp @ fm + model.hm @ fp + model.hx @ f0))
    result = production.frame_rhs(*stages, target_plus=model.target_plus,
                                  target_minus=model.target_minus, **physical_actions(model))
    full = embed_frames(model, (stages[0],) + result)
    for name, actual, expected in zip(('plus', 'minus', 'mixed'), full[1:], dense):
        evidence.close(prefix + '/RHS_' + name, actual, expected, atol=1e-12, rtol=1e-12)
    evidence.reject(prefix + '/omit_Hmixed_drive', dense[2] + 1j * model.hx @ f0, dense[2])
    evidence.reject(prefix + '/omit_reference_connection', dense[2] + 1j * model.h0 @ fx, dense[2])
    return result


def check_borrowed_storage(evidence, model, stages, records):
    originals = tuple(a for stage in stages for a in stage)
    okay = True
    contribution_count = 0
    for blocks in records.values():
        for block in blocks:
            for term in block.dyads:
                contribution_count += 1
                okay = (okay and any(term.left is a for a in originals)
                        and any(term.right is a for a in originals)
                        and term.D is model.D[term.source_index]
                        and term.weight == model.weights[term.source_index])
    evidence.require('density/copy_free_original_source_storage', okay,
                     contribution_count=contribution_count)
    wrong = {name: tuple(replace(block, dyads=tuple(replace(term, weight=model.weights[0])
                                                   for term in block.dyads)) for block in blocks)
             for name, blocks in records.items()}
    evidence.reject('density/replaced_original_source_weights',
                    materialize_descriptors(model, wrong)[3], materialize_descriptors(model, records)[3])


def check_invariants(evidence, model, frame_jets, densities, prefix):
    f0, fp, fm, fx = frame_jets
    gram0 = model.initial_frame.conj().T @ model.initial_frame
    gp = f0.conj().T @ fp + fm.conj().T @ f0
    gm = f0.conj().T @ fm + fp.conj().T @ f0
    gx = (f0.conj().T @ fx + fx.conj().T @ f0
          + fp.conj().T @ fp + fm.conj().T @ fm)
    evidence.close(prefix + '/Gram_mean', f0.conj().T @ f0, gram0)
    for name, value in (('plus', gp), ('minus', gm), ('mixed', gx)):
        evidence.close(prefix + '/Gram_' + name, value, np.zeros_like(value))
    initial_trace = np.trace(model.initial_frame @ model.weighted_D @ model.initial_frame.conj().T)
    evidence.close(prefix + '/trace_mean', np.trace(densities[0]), initial_trace)
    for name, value in zip(('plus', 'minus', 'mixed'), densities[1:]):
        evidence.close(prefix + '/trace_' + name, np.trace(value), 0)
    bad_gx = f0.conj().T @ fx + fx.conj().T @ f0
    evidence.reject(prefix + '/Gram_omit_cross_terms', bad_gx, np.zeros_like(gx))
    evidence.reject(prefix + '/trace_omit_cross_terms',
                    np.trace(densities[3] - fp @ model.weighted_D @ fp.conj().T
                             - fm @ model.weighted_D @ fm.conj().T), 0)


def check_density_mutations(evidence, model, jets, density):
    f0, fp, fm, fx = jets
    d = model.weighted_D
    mixed_mean = fx @ d @ f0.conj().T + f0 @ d @ fx.conj().T
    cross = fp @ d @ fp.conj().T + fm @ d @ fm.conj().T
    evidence.reject('density/omit_mixed_electronic_mean', cross, density[3])
    evidence.reject('density/omit_product_cross_terms', mixed_mean, density[3])
    evidence.reject('density/half_mixed_derivative', density[3] / 2, density[3])
    wrong_bra = fp @ d @ f0.conj().T + f0 @ d @ fp.conj().T
    evidence.reject('density/unswapped_bra_first_direction', wrong_bra, density[1])
    wrong_cross = fp @ d @ fm.conj().T + fm @ d @ fp.conj().T
    evidence.reject('density/unswapped_bra_cross_products', mixed_mean + wrong_cross, density[3])


def check_finite_continuation(evidence, model, time, density):
    steps = (0.04, 0.02, 0.01, 0.005)
    mixed_errors, plus_errors, minus_errors = [], [], []
    for step in steps:
        pp = analytic_density(model, time, step, step)
        pm = analytic_density(model, time, step, -step)
        mp = analytic_density(model, time, -step, step)
        mm = analytic_density(model, time, -step, -step)
        mixed = (pp - pm - mp + mm) / (4 * step * step)
        plus = (analytic_density(model, time, step, 0)
                - analytic_density(model, time, -step, 0)) / (2 * step)
        minus = (analytic_density(model, time, 0, step)
                 - analytic_density(model, time, 0, -step)) / (2 * step)
        mixed_errors.append(float(np.max(np.abs(mixed - density[3]))))
        plus_errors.append(float(np.max(np.abs(plus - density[1]))))
        minus_errors.append(float(np.max(np.abs(minus - density[2]))))
    for name, errors in (('mixed', mixed_errors), ('plus', plus_errors), ('minus', minus_errors)):
        ratios = [a / b for a, b in zip(errors[:-1], errors[1:])]
        evidence.require('continuation/refinement_' + name,
                         all(3.5 < ratio < 4.5 for ratio in ratios) and errors[-1] < 2e-4,
                         steps=list(steps), max_abs_errors=errors, refinement_ratios=ratios)
    lam = 0.08 + 0.11j
    eta = lam.conjugate()
    h = model.h0 + lam * model.hp + eta * model.hm + lam * eta * model.hx
    u = expm(-1j * time * h)
    rho0 = model.initial_frame @ model.weighted_D @ model.initial_frame.conj().T
    evidence.close('continuation/physical_slice_unitary_density',
                   analytic_density(model, time, lam, eta), u @ rho0 @ u.conj().T, atol=1e-12)
    lam, eta = 0.13, -0.07
    h = model.h0 + lam * model.hp + eta * model.hm + lam * eta * model.hx
    u = expm(-1j * time * h)
    evidence.reject('continuation/naive_dagger_off_physical_slice',
                    u @ rho0 @ u.conj().T, analytic_density(model, time, lam, eta))
    step = steps[0]
    finite = (analytic_density(model, time, step, step)
              - analytic_density(model, time, step, -step)
              - analytic_density(model, time, -step, step)
              + analytic_density(model, time, -step, -step)) / (4 * step * step)
    evidence.reject('continuation/half_factor_against_finite_exponentials', density[3] / 2, finite)


def check_target_and_input_mutations(evidence):
    model = finite_model((3, 3, 3), seed=25429)
    full = exact_frame_jets(model, 0.37)
    stages = extract_frames(model, full)
    actions = physical_actions(model)
    p, m = model.target_plus, model.target_minus
    expected = production.frame_rhs(*stages, target_plus=p, target_minus=m, **actions)
    psi, plus, minus, mixed = stages
    wrong_plus = tuple(-1j * (actions['h0'](m[i], plus[i]) + actions['hplus'](i, psi[i]))
                       for i in range(3))
    wrong_mixed = tuple(-1j * (actions['h0'](i, mixed[i])
                                + actions['hplus'](p[i], minus[i])
                                + actions['hminus'](m[i], plus[i])
                                + actions['hmixed'](i, psi[i])) for i in range(3))
    evidence.reject('targets/H0_at_wrong_endpoint', np.concatenate([a.ravel() for a in wrong_plus]),
                    np.concatenate([a.ravel() for a in expected[0]]))
    evidence.reject('targets/swapped_incoming_cross_targets',
                    np.concatenate([a.ravel() for a in wrong_mixed]),
                    np.concatenate([a.ravel() for a in expected[2]]))
    evidence.rejects_input('validation/noninverse_targets', lambda: production.frame_rhs(
        *stages, target_plus=p, target_minus=p, **actions))
    malformed = list(plus)
    malformed[0] = np.zeros((3, model.ranks[0] + 1), complex)
    evidence.rejects_input('validation/lost_original_column_rank', lambda: production.frame_rhs(
        psi, malformed, minus, mixed, target_plus=p, target_minus=m, **actions))
    bad_d = list(model.D)
    bad_d[1] = bad_d[1].copy()
    bad_d[1][0, 1] += 0.2j
    evidence.rejects_input('validation/nonhermitian_D_not_repaired', lambda: production.density_dyads(
        *stages, bad_d, target_plus=p, target_minus=m, source_weights=model.weights))
    evidence.rejects_input('validation/missing_physical_action', lambda: production.frame_rhs(
        *stages, target_plus=p, target_minus=m, **{**actions, 'hmixed': None}))


def run_oracle(evidence):
    time = 0.43
    metadata = []
    for label, model in (('rectangular_q', finite_model()), ('gamma', finite_model(gamma=True, seed=25430))):
        check_model_contract(evidence, model, label)
        stages, integration = evolve_production(model, time, production.frame_rhs)
        actual_jets = embed_frames(model, stages)
        exact = exact_frame_jets(model, time)
        for name, actual, expected in zip(('mean', 'plus', 'minus', 'mixed'), actual_jets, exact):
            evidence.close(label + '/evolved_frame_vs_full_exponential_' + name, actual, expected)
        check_rhs(evidence, model, stages, label)
        records = descriptors(model, stages)
        actual_density = materialize_descriptors(model, records)
        independent_density = dense_density_ode(model, time)
        for name, actual, expected in zip(('mean', 'plus', 'minus', 'mixed'), actual_density, independent_density):
            evidence.close(label + '/density_vs_dense_Liouville_' + name, actual, expected)
        check_invariants(evidence, model, actual_jets, actual_density, label)
        if label == 'rectangular_q':
            check_borrowed_storage(evidence, model, stages, records)
            check_density_mutations(evidence, model, actual_jets, actual_density)
            check_finite_continuation(evidence, model, time, independent_density)
        metadata.append(dict(name=label, dimensions=list(model.dimensions), ranks=list(model.ranks),
                             target_plus=list(model.target_plus), target_minus=list(model.target_minus),
                             weights=list(model.weights), evolution_time=time, integration=integration))
    check_target_and_input_mutations(evidence)
    return metadata


def raw_bindings():
    root = Path(__file__).resolve().parent
    names = ('mixed_frame_response.py', 'mixed_frame_oracle_finite_model.py',
             'run_mixed_frame_response_oracle.py')
    result = []
    for name in names:
        path = root / name
        raw = path.read_bytes()
        result.append(dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path, help='JSON evidence output path')
    args = parser.parse_args()
    evidence = Evidence()
    report = dict(schema='P254-mixed-frame-response-oracle-v1',
                  started_at=datetime.now(timezone.utc).isoformat(),
                  scope='Finite mathematical models only; no native propagation or material energy result.',
                  convention='fixed original D; mixed=d2/(dlambda deta); eta=lambda* physically; bra directions swap',
                  implementation_raw_bytes=raw_bindings(),
                  numpy_version=np.__version__)
    try:
        report['models'] = run_oracle(evidence)
    except Exception as exc:
        report['execution_error'] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
    report['checks'] = evidence.checks
    report['mutations'] = evidence.mutations
    report['passed'] = (not report.get('execution_error') and bool(evidence.checks)
                        and bool(evidence.mutations) and all(c['passed'] for c in evidence.checks)
                        and all(m['rejected'] for m in evidence.mutations))
    report['finished_at'] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(passed=report['passed'], checks=len(evidence.checks),
                          mutations=len(evidence.mutations), output=str(args.output))))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
