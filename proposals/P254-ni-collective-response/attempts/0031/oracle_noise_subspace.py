"""Small owning A31 oracle; importing performs no execution.

Run explicitly with PYTHONPATH=src python <this file>. Exact complex physical
operators and a literal occupation-bit Fock ensemble supply the answer, not an
admitted material output or another raw-factor contraction. Full-complement
projection is compared to physical CAR/Fock covariances including both mixed
retained/omitted orders. Rank-zero/one retain the complete original state but
have a nonzero, directly measured omitted-action selection proxy. No material,
late-time, stationarity, loss or physical-tail-bound claim follows.
"""
from __future__ import annotations

import copy
import json
import numpy as np
import sympy as sp

from noise_subspace import build_noise_subspace
from substrate_framework.physical_factor_preparation import ActionFactors

ATOL = 3e-11
PARTS = ('retained', 'omitted')


def canonical(matrix):
    """Canonicalize each small symbolic block before composing Fock products."""
    return matrix.applyfunc(lambda entry: sp.simplify(sp.expand_complex(entry)))


def annihilators(modes):
    matrices = []
    for mode in range(modes):
        action = sp.zeros(2**modes)
        for bits in range(2**modes):
            if bits & (1 << mode):
                parity = (bits & ((1 << mode)-1)).bit_count()
                action[bits ^ (1 << mode), bits] = (-1)**parity
        matrices.append(action)
    return matrices


def lift(operator, car):
    value = sp.zeros(car[0].rows)
    for a in range(operator.rows):
        for b in range(operator.cols):
            value += operator[a, b]*car[a].H*car[b]
    return canonical(value)


def exact_case():
    i, q = sp.I, sp.Rational
    unitary = sp.Matrix([[1, 0, 1, 0], [i, 0, -i, 0],
                         [0, 1, 0, 1], [0, i, 0, -i]]) / sp.sqrt(2)
    raw_map = sp.Matrix([[1, i/3], [q(1, 4), q(6, 5)]])
    prepared = sp.Matrix([[q(2, 5), i/10], [-i/10, q(1, 2)]])
    psi = canonical(unitary[:, :2]*raw_map)
    inverse = canonical(raw_map.inv())
    density = canonical(inverse*prepared*inverse.H)
    coordinates = inverse[:, :1]
    forces = {
        'mode0': sp.Matrix([[1, 1+i, 2*i, q(1, 3)], [2-i, -1, i, 2],
                            [1, 3, 2, i], [2*i, -1, 1-i, 3]]),
        'mode1': sp.Matrix([[i, 2, 1-i, -1], [-1, i, 2*i, 3],
                            [2+i, 1, 3, -i], [1, 1+i, 2, -2]])}
    seed = sp.Matrix([[1, i, 2, 1-i], [2, -1, i, 3],
                      [-i, 1, 2, 1], [3, i, 2-i, -2]])
    hamiltonian = canonical((seed+seed.H)/2)
    gram = canonical(psi.H*psi)
    rho = canonical(psi*density*psi.H)
    active = canonical(psi*coordinates)
    projector = canonical(active*active.H)
    assert canonical(rho*projector-projector*rho) != sp.zeros(4)
    assert canonical(coordinates.H*gram*coordinates) == sp.eye(1)
    assert canonical(unitary.H*unitary) == sp.eye(4)
    # Exact two-orbital quasifree ensemble, with the two other modes vacuum.
    # Its Fock weights/coherences are constructed directly, without Wick's rule.
    ensemble = sp.zeros(16)
    double = sp.simplify(prepared.det())
    ensemble[0, 0] = sp.simplify((sp.eye(2)-prepared).det())
    ensemble[3, 3] = double
    one_particle = canonical(prepared-double*sp.eye(2))
    for a in range(2):
        for b in range(2):
            ensemble[1 << a, 1 << b] = one_particle[a, b]
    assert sp.simplify(sp.trace(ensemble)-1) == 0
    car = annihilators(4)
    pieces, lifted, means = {}, {}, {}
    for name, operator in forces.items():
        retained = canonical(projector*operator*projector)
        for part, value in (('full', operator), ('retained', retained),
                            ('omitted', canonical(operator-retained))):
            key = name, part
            pieces[key] = value
            lifted[key] = lift(canonical(unitary.H*value*unitary), car)
            means[key] = sp.simplify(sp.trace(ensemble*lifted[key]))
            assert sp.simplify(means[key]-sp.trace(rho*value)) == 0
    covariances = {}
    for a in forces:
        for b in forces:
            for pa, pb in [('full', 'full')] + [(x, y) for x in PARTS for y in PARTS]:
                key = a, pa, b, pb
                left, right = (a, pa), (b, pb)
                exact = sp.simplify(sp.trace(ensemble*lifted[left]*lifted[right])
                                    - means[left]*means[right])
                physical = sp.trace(rho*pieces[left]*(sp.eye(4)-rho)*pieces[right])
                assert sp.simplify(exact-physical) == 0
                covariances[key] = exact
    cross = [sp.simplify(covariances[a, 'retained', b, 'omitted']
                        + covariances[a, 'omitted', b, 'retained'])
             for a in forces for b in forces]
    assert any(value != 0 for value in cross)
    numeric = lambda value: np.array(value.evalf(), dtype=np.complex128)
    case = dict(psi=numeric(psi), gram=numeric(gram), d=numeric(density),
                r=numeric(coordinates), rho=numeric(rho), projector=numeric(projector),
                h=numeric(hamiltonian), forces={k: numeric(v) for k, v in forces.items()})
    reference = {key: complex(value) for key, value in covariances.items()}
    return case, reference


def build(case, rank, *, false_adjoint=False, omit_extra=False):
    operators = {'hamiltonian': case['h'], **case['forces']}
    factors = {}
    for name, operator in operators.items():
        adjoint = operator if false_adjoint and name != 'hamiltonian' else operator.conj().T
        factors[name] = ActionFactors(operator@case['psi'], adjoint@case['psi'],
            'exact-fixture:'+name, 'exact-fixture:'+name+':independent-adjoint',
            'fixture-energy', 'four-physical-orbitals', 'literal exact finite fixture')
    calls = []

    def apply(name, vectors, adjoint=False):
        # Actual callback input is tested geometrically, not by a copied key.
        assert vectors.shape == (4, rank)
        np.testing.assert_allclose(case['psi'].conj().T@vectors, 0, atol=ATOL)
        calls.append((name, adjoint, vectors.shape[1]))
        operator = operators[name]
        if adjoint and not (false_adjoint and name != 'hamiltonian'):
            operator = operator.conj().T
        return np.zeros_like(vectors) if omit_extra else operator@vectors

    result = build_noise_subspace(case['psi'], case['gram'], case['d'], case['r'],
        factors['hamiltonian'], {k: factors[k] for k in case['forces']}, apply,
        extra_rank=rank)
    # A rank-zero invocation must not ask a physical operator to act on Psi.
    assert len(calls) == (2*len(operators) if rank else 0)
    return result


def close(actual, expected):
    np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=ATOL)


def physical_basis(case, result):
    # Only the tiny oracle forms Q0, the physical projector and dense operators.
    return np.column_stack((case['psi']@result['inverse_sqrt_gram'], result['extra_basis']))


def projected_covariances(result):
    rho = result['projected_density']
    projector = result['active_projector']
    matrices, adjoints = {}, {}
    for name, operator in result['forces'].items():
        independent_adjoint = result['force_adjoints'][name]
        retained = projector@operator@projector
        retained_adjoint = projector@independent_adjoint@projector
        for part, value, adjoint in (
                ('full', operator, independent_adjoint),
                ('retained', retained, retained_adjoint),
                ('omitted', operator-retained, independent_adjoint-retained_adjoint)):
            matrices[name, part] = value
            adjoints[name, part] = adjoint
    values = {}
    for a in result['forces']:
        for b in result['forces']:
            for pa, pb in [('full', 'full')] + [(x, y) for x in PARTS for y in PARTS]:
                left, right = (a, pa), (b, pb)
                # The independently supplied true-adjoint leg must participate.
                values[a, pa, b, pb] = np.trace(rho@adjoints[left].conj().T@matrices[right]) \
                    - np.trace(rho@matrices[left]@rho@matrices[right])
    return values


def check_geometry_and_operators(case, result):
    basis = physical_basis(case, result)
    close(basis.conj().T@basis, np.eye(basis.shape[1]))
    close(basis@result['projected_frame'], case['psi'])
    close(basis@result['projected_density']@basis.conj().T, case['rho'])
    close(basis@result['active_projector']@basis.conj().T, case['projector'])
    operators = {'hamiltonian': case['h'], **case['forces']}
    for name, physical in operators.items():
        forward = result['hamiltonian'] if name == 'hamiltonian' else result['forces'][name]
        adjoint = result['hamiltonian_adjoint'] if name == 'hamiltonian' else result['force_adjoints'][name]
        close(forward, basis.conj().T@physical@basis)
        close(adjoint, basis.conj().T@physical.conj().T@basis)
    return basis


def direct_proxy(case):
    # Independent physical-row construction uses sqrt(|rho|) rather than raw
    # coordinate transport. It has the same singular values (additional zero
    # columns only), and does not share the production Gram weighting formula.
    original, _ = np.linalg.qr(case['psi'], mode='reduced')
    # Work on the exactly declared occupied support, not spurious roundoff
    # eigenvalues of its two known vacuum modes. No state eigenvalue is clipped.
    values, vectors = np.linalg.eigh(original.conj().T@case['rho']@original)
    orbitals = original@vectors
    root = (orbitals*np.sqrt(np.abs(values)))@orbitals.conj().T
    complement = np.eye(4)-original@original.conj().T
    return np.column_stack([complement@operator@root
        for force in case['forces'].values() for operator in (force, force.conj().T)])


def proxy_checks(case, result, basis):
    proxy = direct_proxy(case)
    singular_values = np.linalg.svd(proxy, compute_uv=False)
    close(result['singular_values'], singular_values)
    residual = proxy-basis@(basis.conj().T@proxy)
    close(result['diagnostics']['discarded_proxy_frobenius'], np.linalg.norm(residual))
    close(result['diagnostics']['discarded_proxy_spectral'], np.linalg.norm(residual, ord=2))
    expected = np.linalg.eigvalsh(case['rho'])
    # The physical empty complement contributes exactly P-r additional zeros.
    actual = np.sort(np.r_[result['state_eigenvalues'], np.zeros(2)])
    close(actual, expected)
    return float(np.linalg.norm(residual))


def must_expose(name, operation):
    try:
        operation()
    except AssertionError:
        return name
    raise AssertionError('oracle failed to expose '+name)


def run_oracle():
    case, reference = exact_case()
    results, omitted = {}, {}
    for rank in (0, 1, 2):
        result = build(case, rank)
        basis = check_geometry_and_operators(case, result)
        omitted[rank] = proxy_checks(case, result, basis)
        results[rank] = result
    assert omitted[0] > omitted[1] > 1e-3
    assert omitted[2] < ATOL
    actual = projected_covariances(results[2])
    for key, value in reference.items():
        close(actual[key], value)
    # A raw Hermitian but indefinite state is still transported unchanged: its
    # negative eigenvalue is a selection weight, not clipped into a new state.
    signed = copy.deepcopy(case)
    vector = np.array([1, 1j], dtype=complex)
    signed['d'] = case['d']-np.outer(vector, vector.conj())
    signed['rho'] = signed['psi']@signed['d']@signed['psi'].conj().T
    signed_result = build(signed, 1)
    signed_basis = check_geometry_and_operators(signed, signed_result)
    proxy_checks(signed, signed_result, signed_basis)
    assert signed_result['state_eigenvalues'].min() < -0.1
    assert not signed_result['diagnostics']['state_within_pauli_roundoff']

    mutations = []
    wrong_gram = copy.deepcopy(case)
    r = case['r'][:, 0]
    perpendicular = np.array([-r[1].conjugate(), r[0].conjugate()])
    wrong_gram['gram'] += .2*np.outer(perpendicular, perpendicular.conj())
    mutations.append(must_expose('wrong supplied original Gram', lambda:
        check_geometry_and_operators(case, build(wrong_gram, 0))))
    wrong_d = copy.deepcopy(case)
    wrong_d['d'] = np.linalg.inv(case['gram'])@case['d']@np.linalg.inv(case['gram'])
    mutations.append(must_expose('wrong raw-D metric transport', lambda:
        check_geometry_and_operators(case, build(wrong_d, 2))))
    mutations.append(must_expose('forward action substituted for true adjoint', lambda:
        check_geometry_and_operators(case, build(case, 2, false_adjoint=True))))
    mutations.append(must_expose('extra-column physical callback omitted', lambda:
        check_geometry_and_operators(case, build(case, 2, omit_extra=True))))
    wrong_mixed = copy.deepcopy(results[2])
    original_rank = case['psi'].shape[1]
    for matrix in wrong_mixed['forces'].values():
        matrix[:original_rank, original_rank:] = matrix[original_rank:, :original_rank].conj().T
    mutations.append(must_expose('false Hermitian mixed operator block', lambda:
        check_geometry_and_operators(case, wrong_mixed)))

    def missing_covariance_cross():
        for a in case['forces']:
            for b in case['forces']:
                diagonal_only = actual[a, 'retained', b, 'retained']+actual[a, 'omitted', b, 'omitted']
                close(diagonal_only, reference[a, 'full', b, 'full'])
    mutations.append(must_expose('retained/omitted covariance cross terms dropped', missing_covariance_cross))
    return dict(scope='exact four-orbital initial finite-subspace CAR/Fock oracle only',
                exact_fock_directed_covariances=len(reference),
                ranks=[0, 1, 2], omitted_action_proxy=omitted,
                full_complement_max_covariance_error=float(max(
                    abs(actual[key]-value) for key, value in reference.items())),
                exposed_mutations=mutations,
                signed_state='negative eigenvalue retained; absolute value used only for selection')


if __name__ == '__main__':
    print(json.dumps(run_oracle(), indent=2, allow_nan=False))
