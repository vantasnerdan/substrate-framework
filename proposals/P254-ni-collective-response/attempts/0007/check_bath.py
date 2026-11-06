"""Expose two-point oscillator license and ensemble-sensitive static subtraction."""
import numpy as np
import sympy as sp

from substrate_framework.finite_fermion_state import fermi_dirac
from substrate_framework.finite_fermion_response import noise_lines, retarded, isolated_static
from substrate_framework.finite_fermion_bath import oscillator_map, oscillator_noise, oscillator_retarded
from substrate_framework.finite_fermion_thermodynamics import grandcanonical_curvature, calibrated_inverse
from substrate_framework.material_matter_field import static_matched_bath_hessian, eliminate_bath
from fixtures import OPERATORS, ENERGIES, BETA, HBAR, Z, zero, serial
from free_energy_oracle import free_energy_functions, negative_hessian


def run(ledger):
    outputs = []
    calibrated = sp.Matrix([[17, 2], [2, 19]])
    for energies in ENERGIES:
        f = fermi_dirac(energies, BETA)
        modes = oscillator_map(energies, f, OPERATORS)
        noise = oscillator_noise(modes, BETA, 2)
        electronic = noise_lines(energies, f, OPERATORS)
        ledger.check('complex Gaussian oscillator noise matches all nonzero fermion lines',
                     set(noise) == set(electronic)-{0} and all(
                         zero(noise[g]-electronic[g]) for g in noise))
        chi = retarded(energies, f, OPERATORS, Z, HBAR)
        bath_chi = oscillator_retarded(modes, Z, HBAR, 2)
        ledger.check('both oscillator poles reproduce full retarded matrix', zero(chi-bath_chi))
        real_modes = [(g, coupling.applyfunc(sp.re)) for g, coupling in modes]
        ledger.check('discarded complex coupling phases mutation rejected',
                     not zero(oscillator_retarded(real_modes, Z, HBAR, 2)-chi))
        vacuum_noise = {gap: coupling*coupling.H for gap, coupling in modes}
        ledger.check('zero Bose population mutation rejected at finite temperature',
                     any(not zero(vacuum_noise[g]-noise[g]) for g in vacuum_noise))
        ledger.check('positive-gap oscillator noise misses zero occupation sector',
                     0 not in noise and not zero(electronic[0]))
        static = isolated_static(energies, f, OPERATORS)
        gc = grandcanonical_curvature(energies, f, OPERATORS, BETA)
        qpcoupling = sp.Matrix.hstack(*(sp.Matrix.hstack(
            sp.sqrt(2)*g.applyfunc(sp.re), -sp.sqrt(2)*g.applyfunc(sp.im)) for _, g in modes))
        bath_stiffness = sp.diag(*(v for gap, _ in modes for v in (gap, gap)))
        ledger.check('complex oscillator canonical quadrature static Schur term',
                     zero(qpcoupling*bath_stiffness.inv()*qpcoupling.T-static))
        full = static_matched_bath_hessian(calibrated, bath_stiffness, qpcoupling)
        ledger.check('matched isolated quadratic action counterterm returns calibrated stiffness',
                     zero(eliminate_bath(full, 2, 0)-calibrated))
        valid_gc = calibrated_inverse(calibrated, gc, gc, sp.eye(2), 0)
        ledger.check('same-ensemble GC action subtraction returns calibrated stiffness',
                     zero(valid_gc-calibrated))
        invalid = calibrated_inverse(calibrated, gc, static, sp.eye(2), 0)
        defect = invalid-calibrated
        ledger.check('isolated counterterm used in GC action mutation exposes population defect',
                     not zero(defect) and zero(defect+BETA*electronic[0]))
        # Actual free-energy action: 1/2 x (K_cal+C_gc) x + Omega(x)-Omega(0).
        grand = free_energy_functions([float(e) for e in energies])[0]
        ledger.check('action curvatures are exactly real', zero(gc-gc.conjugate())
                     and zero(static-static.conjugate()))
        bare = np.array((calibrated+gc).applyfunc(sp.simplify).evalf(), float)
        wrong_bare = np.array((calibrated+static).applyfunc(sp.simplify).evalf(), float)
        def action(x):
            return x@bare@x/2+grand(x)
        def wrong_action(x):
            return x@wrong_bare@x/2+grand(x)
        observed = -negative_hessian(action, 0.0005)
        wrong_observed = -negative_hessian(wrong_action, 0.0005)
        ledger.check('GC matched action verified by independent full Fock free-energy derivatives',
                     np.max(np.abs(observed-np.array(calibrated, float))) < 0.0002)
        ledger.check('invalid subtraction independently exposed in actual free-energy action',
                     np.max(np.abs(wrong_observed-np.array(calibrated, float))) > 1e-3)
        outputs.append({'energies': energies, 'modes_gap_complex_g': modes,
                        'oscillator_noise': noise, 'missing_zero_noise': electronic[0],
                        'isolated_static_counterterm': static, 'GC_static_counterterm': gc,
                        'invalid_GC_action_minus_calibrated': defect,
                        'independent_GC_action_hessian': observed.tolist(),
                        'independent_wrong_action_hessian': wrong_observed.tolist()})
    gap = sp.symbols('Delta', positive=True)
    fi = sp.Rational(1, 2)
    fj = 1/(1+sp.exp(gap))
    strength = fi-fj
    positive_noise_limit = sp.limit(strength/(1-sp.exp(-gap)), gap, 0, dir='+')
    static_limit = sp.limit(2*strength/gap, gap, 0, dir='+')
    dynamic_limit = sp.limit(strength*(1/(gap-HBAR*Z)+1/(gap+HBAR*Z)), gap, 0, dir='+')
    ledger.check('coalescing-gap finite noise and noncommuting static/dynamic limits',
                 positive_noise_limit == sp.Rational(1, 4)
                 and static_limit == sp.Rational(1, 2) and dynamic_limit == 0)
    return serial({'fixtures': outputs, 'zero_gap_limit': {
        'one_positive_line_noise': positive_noise_limit, 'static_first_then_gap': static_limit,
        'gap_first_at_nonzero_z': dynamic_limit,
        'exception': 'Exactly zero gaps have finite covariance but no retarded commutator; not a normalized Delta=0 thermal oscillator.'}})
