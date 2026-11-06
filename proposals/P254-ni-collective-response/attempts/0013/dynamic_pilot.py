"""Compare one finite functional's direct dynamics with two screening choices."""
from __future__ import annotations

import json
import numpy as np
import sympy as sp

from fixture import make_fixture, symbolic_reference
from substrate_framework.finite_fermion_response import retarded
from substrate_framework.finite_fermion_thermodynamics import grandcanonical_curvature
from substrate_framework.screened_force_response import (
    linearized_density_response, screened_response, static_screened_mixed)


def joint_response(fixture, z):
    operators = tuple(sp.Matrix(a) for a in (*fixture['B'], *fixture['D']))
    energies, occupations = symbolic_reference()
    return np.array(retarded(energies, occupations, operators, z)).astype(complex)


def constrained_static(fixture, fixed_m=False):
    operators = tuple(sp.Matrix(a) for a in (*fixture['B'], *fixture['D']))
    constraints = (sp.eye(4), sp.Matrix(fixture['Z'])) if fixed_m else (sp.eye(4),)
    full = np.array(grandcanonical_curvature(
        *symbolic_reference(), operators + constraints, sp.S.One)).astype(complex)
    return (full[:4, :4] - full[:4, 4:]
            @ np.linalg.solve(full[4:, 4:], full[4:, :4]))


def serialize(matrix):
    matrix = np.asarray(matrix)
    return {'real': matrix.real.tolist(), 'imag': matrix.imag.tolist()}


def exact_frequency_freezing_witness():
    """Scalar channel algebra at arbitrary independent chi(z), chi_static."""
    c, c0, kernel = sp.symbols('c c0 kernel', real=True)
    b = sp.symbols('b', real=True)
    full = b**2 * c / (1 + kernel * c)
    frozen = b**2 * c / (1 + kernel * c0)
    difference = sp.factor(full - frozen)
    expected = b**2 * c * kernel * (c0 - c) / (
        (1 + kernel*c) * (1 + kernel*c0))
    assert sp.cancel(difference - expected) == 0
    witness = sp.factor(difference.subs({b: 2, c: sp.Rational(1, 3),
                                       c0: sp.Rational(1, 2), kernel: 1}))
    assert witness != 0
    assert sp.simplify(difference.subs(c, c0)) == 0
    return {'difference': str(difference), 'exact_witness': str(witness),
            'scope': 'Conditional finite scalar identity, not a material error estimate.'}


def probe(fixture, z):
    joint = joint_response(fixture, z)
    stat = constrained_static(fixture)
    full = screened_response(joint, 2, fixture['U'])
    frozen = static_screened_mixed(joint, stat, 2, fixture['U'])
    direct = linearized_density_response(
        np.diag(fixture['h0']), np.diag(fixture['gamma0']),
        fixture['B'], fixture['D'], fixture['U'], z)
    error = float(np.max(np.abs(full - direct)))
    assert error < 2e-12, ('direct dynamics disagrees', error)
    discrepancy = float(np.max(np.abs(frozen - direct)))
    assert discrepancy > 1e-6, ('frequency freezing witness insensitive', discrepancy)
    wrong_sign = screened_response(joint, 2, -fixture['U'])
    sign_error = float(np.max(np.abs(wrong_sign - direct)))
    assert sign_error > 1e-6, ('interaction-sign mutation insensitive', sign_error)
    # Complex same-eigenbasis phases test contraction order, not native finite-q.
    phases = np.exp(1j * np.array([0.0, .31, -.23, .59]))
    gauge = np.diag(phases)
    changed = dict(fixture)
    for key in ('B', 'D'):
        changed[key] = np.array([gauge.conj().T @ a @ gauge for a in fixture[key]])
    gauged = screened_response(joint_response(changed, z), 2, fixture['U'])
    gauge_error = float(np.max(np.abs(gauged - full)))
    assert gauge_error < 2e-12, ('gauge contraction error', gauge_error)
    return {'z': {'real': z.real, 'imag': z.imag}, 'selfconsistent': serialize(full),
            'direct_density_equation': serialize(direct), 'frozen_static_vertex': serialize(frozen),
            'full_vs_direct_max_error': error,
            'frozen_vs_direct_max_discrepancy': discrepancy,
            'wrong_U_sign_discrepancy': sign_error, 'complex_gauge_error': gauge_error,
            'frozen_cross_channel_asymmetry': float(np.max(np.abs(frozen-frozen.T)))}


def main():
    fixture = make_fixture()
    finite = probe(fixture, 1/3 + 1j/5)
    joint0 = joint_response(fixture, 0)
    isolated0 = screened_response(joint0, 2, fixture['U'])
    isothermal = screened_response(constrained_static(fixture), 2, fixture['U'])
    ensemble_gap = float(np.max(np.abs(isothermal - isolated0)))
    assert ensemble_gap > 1e-6, 'population/static ensemble witness insensitive'
    result = {'kind': 'finite-functional-screening-construction',
              'exact_frequency_freezing': exact_frequency_freezing_witness(),
              'finite_frequency': finite,
              'isolated_zero': serialize(isolated0),
              'fixed_mean_N_isothermal': serialize(isothermal),
              'isothermal_minus_isolated_max': ensemble_gap,
              'decision_scope': 'Full same-functional induced density solves the direct dynamics. '
                  'Frozen static screening is assessed as an approximation, not refuted as an approximation. '
                  'Toy success is not native SG15/PBE transfer, q admission, physical loss or production authority.'}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
