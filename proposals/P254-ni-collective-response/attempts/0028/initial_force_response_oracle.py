"""Sensitive initial-jet oracle: symbolic moving coefficients and full coupled RHS."""
import json
import sys
from substrate_framework.streamed_material_state import file_binding
import numpy as np
import sympy as sp

from initial_force_response import force_rate, full_grid
from run_joint_response_oracle import FiniteSystem
from substrate_framework.joint_response_observables import local_momentum_kick, total_force_probe, readout
from substrate_framework.joint_spatial_readout import coherent_local


def symbolic_block():
    t = sp.symbols('t', real=True)
    c0, c1, g0, g1, b0, b1, p, mass, rho2, u2, a2 = sp.symbols('C0 C1 G0 G1 B0 B1 P M rho2 u2 a2', nonzero=True)
    u = p*t/mass+u2*t**2/2
    rho, a = rho2*t**2/2, a2*t**2/2
    force = -(c0+c1*t)*u-(g0+g1*t)*rho-(b0+b1*t)*a
    assert force.subs(t, 0) == 0
    assert sp.simplify(sp.diff(force, t).subs(t, 0)+c0*p/mass) == 0
    n, phase = sp.symbols('N phase', positive=True)
    assert sp.simplify((-c0*p/mass).subs(p, phase/sp.sqrt(n))+c0*phase/(mass*sp.sqrt(n))) == 0
    return 'Exact first derivative with independent time-dependent C/G/field coefficients; no stationary-reference assumption.'


def coupled_block():
    system = FiniteSystem()
    checks = []
    for gamma in (True, False):
        response = system.response(gamma)
        cell = np.array([0, 0, 1])
        direction = np.array([0., 1., 0.])
        initial = local_momentum_kick(response, cell, direction, 1.)
        expected = force_rate(response.curvature, response.masses, response.q_fractional,
            ncells=response.ncells, cell=cell, amplitude=1., coordinate=1)
        for time in (0., .13):
            stage = system.stage(time)
            derivative = response.rhs(stage, initial)
            actual = []
            for axis in range(response.ndof):
                target = np.eye(response.ndof)[axis]
                probe = total_force_probe(response, target, identity='finite initial-force jet')
                np.testing.assert_allclose(readout(response, stage, initial, probe)['total'], 0, atol=1e-14)
                actual.append(readout(response, stage, derivative, probe)['total'])
            np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
        checks.append(dict(gamma=gamma, actual_coupled_rhs='matched at two distinct moving mean stages'))
        if not gamma:
            for wrong in (expected*np.sqrt(response.ncells), expected.conj(),
                    -response.curvature.T[:, 1]*initial.p[1]/response.masses[1]):
                if np.allclose(wrong, expected, rtol=1e-8, atol=1e-10):
                    raise AssertionError('normalization, opposite phase, or transpose mutation was insensitive')
    return checks


def spatial_block():
    rng = np.random.default_rng(225386)
    mesh = np.array([3, 2, 1])
    c = rng.normal(size=(*mesh, 3, 3))+1j*rng.normal(size=(*mesh, 3, 3))
    cell = np.array([1, 1, 0])
    result = full_grid(c, [2., 3., 5.], mesh=mesh, cell=cell, amplitude=1.7, coordinate=2)
    direct = np.column_stack([coherent_local(result['q_fractional'], values, result['cells'])
                             for values in result['force_rate_q'].T])
    np.testing.assert_allclose(result['force_rate_cells'], direct, rtol=2e-13, atol=2e-13)
    # Independent real-space convolution, not a round-trip of the same helper.
    lattice = np.fft.ifftn(c, axes=(0, 1, 2))
    expected = np.array([-lattice[tuple((r-cell) % mesh)][:, 2]*1.7/5 for r in result['cells']])
    np.testing.assert_allclose(result['force_rate_cells'], expected, rtol=2e-13, atol=2e-13)
    bad = np.fft.ifftn(result['force_rate_q'].reshape(*mesh, 3), axes=(0, 1, 2)).reshape(-1, 3)
    if np.allclose(bad, expected, rtol=1e-8, atol=1e-10):
        raise AssertionError('extra 1/sqrtN spatial normalization was insensitive')
    return 'Complete unequal3x2x1 character grid, nonzero seed cell, complex nonreciprocal raw C and independent convolution; no repair.'


if __name__ == '__main__':
    sources = [file_binding(__file__)]+[file_binding(sys.modules[name].__file__) for name in (
        'initial_force_response', 'run_joint_response_oracle',
        'substrate_framework.joint_bloch_response',
        'substrate_framework.joint_response_observables',
        'substrate_framework.joint_spatial_readout')]
    result = dict(schema='P254-exact-initial-force-response-oracle-v1',
        symbolic=symbolic_block(), coupled=coupled_block(), spatial=spatial_block(),
        scope='Exact initial-jet algebra and original probe/RHS conventions only; not finite-time or material accuracy.')
    if sources != [file_binding(source['path']) for source in sources]:
        raise ValueError('oracle source changed during verification')
    result['source_bindings'] = sources
    result['status'] = 'passed'
    print(json.dumps(result, indent=2))
