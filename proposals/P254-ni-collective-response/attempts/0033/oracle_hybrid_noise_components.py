"""New saved-coordinate/source-attribution and finite-channel influence oracle.

Only these uncertain consumer boundaries are exercised. Existing A29/A32
Hamiltonian/dynamics proofs and all material trajectories are carried forward.
"""
import json
import numpy as np
from scipy.linalg import expm

from oracle_hybrid_gamma_model import fixture, fock_covariance, close
from hybrid_noise_components import decompose_initial_noise
from finite_bath_influence import frozen_force_influence


def mechanical_fixture():
    # Three oscillator vacuum modes, represented in a two-level tensor space.
    # Its first/second moments are exact; no high occupation is truncated.
    a = np.asarray([[0., 1.], [0., 0.]])
    lowering = []
    for mode in range(3):
        operator = np.ones((1, 1))
        for axis in range(3):
            operator = np.kron(operator, a if axis == mode else np.eye(2))
        lowering.append(operator)
    q = [(operator+operator.T)/np.sqrt(2) for operator in lowering]
    p = [(operator-operator.T)/(1j*np.sqrt(2)) for operator in lowering]
    original = np.asarray([q[0], q[1], p[0], p[1], 2*q[2], p[2]/2])
    transform = np.eye(6)
    cosine, sine = np.cos(.37), np.sin(.37)
    for first, second in ((0, 4), (2, 5)):
        transform[first, first] = transform[second, second] = cosine
        transform[first, second] = sine
        transform[second, first] = -sine
    operators = np.einsum('ij,jab->iab', transform, original)
    covariance = np.asarray([[(left@right)[0, 0] for right in operators]
                             for left in operators])
    return operators, covariance


def source_attribution():
    model, c, bank, directions, rng = fixture()
    mechanical, w0 = mechanical_fixture()
    duals = rng.normal(size=(model.size, 3))+1j*rng.normal(size=(model.size, 3))
    coordinates = dict(right_eigenvectors=model.v, inverse_right_eigenvectors=model.vi)
    def calculate(values, saved_coordinates):
        return decompose_initial_noise(values, w0, c.d, c.rho0, bank,
            saved_coordinates, source_slot=0, factor_identity='synthetic-original-frame',
            basis_identity='synthetic-finite-bath', state_authority='synthetic-original-D')
    actual = calculate(duals, coordinates)
    nx, m, b = model.nx, c.m, model.b
    expected = np.zeros((3, 3), complex)
    physical = []
    for probe in range(3):
        z = duals[-b*b:, probe].reshape(b, b)
        physical.append(model.v@z.T@model.vi)
    mechanical_error = 0.
    for left in range(3):
        for right in range(3):
            left_operator = np.einsum('i,iab->ab', duals[:nx, left], mechanical)
            right_operator = np.einsum('i,iab->ab', duals[:nx, right], mechanical)
            expected[left, right] = (left_operator@right_operator)[0, 0]
            for a, ia in [('ionic', range(4)), ('field', range(4, 6))]:
                for bname, ib in [('ionic', range(4)), ('field', range(4, 6))]:
                    la = sum(duals[i, left]*mechanical[i] for i in ia)
                    rb = sum(duals[i, right]*mechanical[i] for i in ib)
                    mechanical_error = max(mechanical_error, close(
                        actual['mechanical'][a+'_'+bname][left, right],
                        (la@rb)[0, 0], 'literal oscillator cross-order'))
            expected[left, right] += fock_covariance(bank['projected_density'],
                physical[left], physical[right])
            for slot in range(1, c.ns):
                start = nx+(slot-1)*m*m
                ol = duals[start:start+m*m, left].reshape(m, m).T
                ort = duals[start:start+m*m, right].reshape(m, m).T
                entry = fock_covariance(c.rho0[slot], ol, ort)
                close(actual['electronic_by_slot'][slot, left, right], entry,
                      'literal per-sector Fock covariance')
                expected[left, right] += entry
    total_error = close(actual['total'], expected, 'independent product Fock/oscillator budget')
    q = bank['projected_active_coordinates']; projector = q@q.conj().T
    retained = [projector@operator@projector for operator in physical]
    omitted = [operator-part for operator, part in zip(physical, retained)]
    for a, lefts in [('retained', retained), ('omitted', omitted)]:
        for bname, rights in [('retained', retained), ('omitted', omitted)]:
            direct = np.asarray([[fock_covariance(bank['projected_density'], l, r)
                                  for r in rights] for l in lefts])
            close(actual['source0'][a+'_'+bname], direct, 'direct split Fock cross-order')
    gauge = np.diag([.7, 1.3, np.exp(.4j), 1.6])
    inverse_gauge = np.linalg.inv(gauge)
    changed = duals.copy()
    z = np.moveaxis(changed[-model.b**2:].reshape(model.b, model.b, 3), -1, 0)
    z = gauge.T@z@inverse_gauge.T
    changed[-model.b**2:] = np.moveaxis(z, 0, -1).reshape(model.b**2, 3)
    same = calculate(changed, dict(right_eigenvectors=model.v@gauge,
        inverse_right_eigenvectors=inverse_gauge@model.vi))
    gauge_error = close(same['total'], expected, 'actual persisted eigen-gauge invariance')
    cross = actual['mechanical']['ionic_field']+actual['mechanical']['field_ionic']
    split_cross = actual['source0']['retained_omitted']+actual['source0']['omitted_retained']
    if np.max(np.abs(cross)) < 1e-6 or np.max(np.abs(split_cross)) < 1e-6:
        raise AssertionError('omitted mechanical/electronic cross-order mutations escaped')
    wrong = calculate(changed, coordinates)
    if np.max(np.abs(wrong['total']-expected)) < 1e-6:
        raise AssertionError('wrong saved eigencoordinate bank mutation escaped')
    return dict(total_error=total_error, mechanical_cross_error=mechanical_error,
        eigen_gauge_error=gauge_error, missing_mechanical_cross_separation=float(np.max(np.abs(cross))),
        missing_source0_cross_separation=float(np.max(np.abs(split_cross))))


def force_influence():
    model, c, original_bank, directions, rng = fixture()
    bank = dict(original_bank)
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    h = np.diag([.1, .4, .8, 1.3])+.05*(raw+raw.conj().T)
    bank['hamiltonian'] = h
    for mode in range(2):
        raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        force = (raw+raw.conj().T)/7
        bank['force_mode'+str(mode)] = force
        bank['adjoint_force_mode'+str(mode)] = force.conj().T
    energies, vectors = np.linalg.eig(h)
    coordinates = dict(eigenvalues=energies, right_eigenvectors=vectors,
        inverse_right_eigenvectors=np.linalg.inv(vectors))
    times = np.asarray([0., .4, 1.1]); scale = .5
    actual = frozen_force_influence(bank, coordinates, times, source_factor=scale,
        factor_identity='synthetic-original-frame', basis_identity='synthetic-finite-bath',
        state_authority='synthetic-original-D')
    expected = np.empty((3, 2, 3, 2), complex)
    for ti, time in enumerate(times):
        for mode in range(2):
            left = expm(1j*h*time)@(scale*bank['force_mode'+str(mode)])@expm(-1j*h*time)
            for tj, second_time in enumerate(times):
                for second_mode in range(2):
                    right = expm(1j*h*second_time)@(scale*bank['force_mode'+str(second_mode)])@expm(-1j*h*second_time)
                    expected[ti, mode, tj, second_mode] = fock_covariance(
                        bank['projected_density'], left, right)
    error = close(actual['ordered'], expected, 'literal nonstationary Fock influence')
    reverse = expected.transpose(2, 3, 0, 1)
    theta = np.where(times[:, None] > times[None, :], 1.,
                     np.where(times[:, None] == times[None, :], .5, 0.))
    response = 1j*theta[:, None, :, None]*(expected-reverse)
    response_error = close(actual['retarded'], response, 'physical -uF Kubo sign/causality')
    if np.max(np.abs(2*response)) < 1e-7:
        raise AssertionError('retarded-sign mutation escaped')
    nonstationarity = np.max(np.abs(expected[0, :, 0, :]-expected[-1, :, -1, :]))
    if nonstationarity < 1e-7:
        raise AssertionError('stationary time-difference substitution mutation escaped')
    return dict(ordered_error=error, retarded_error=response_error,
                equal_elapsed_time_nonstationarity_separation=float(nonstationarity))


if __name__ == '__main__':
    print(json.dumps(dict(schema='P254-A33-new-consumer-owning-oracle-v1',
        source_budget=source_attribution(), force_influence=force_influence(),
        scope='new source/cross-order/saved-coordinate and nonstationary influence consumers only; no material calibration'), sort_keys=True))
