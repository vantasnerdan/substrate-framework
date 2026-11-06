"""Independent finite SymPy, dense-PW and Fock oracles for the t0 split.

Run explicitly with PYTHONPATH=src python <this file>; importing runs nothing.
The production path uses small full-row moments and FactorPreparation. This
oracle instead constructs tiny PW-square projectors and actual fermionic Fock
operators/ensembles. Exact rational complex data test the encoded CAR identity;
independent seeded complex frames test the numerical public interface and raw-D
transport. Neither is material-source authentication or full128/64q, late-time
loss, damping, bath fitting, FDT or lifetime evidence.

Convention sources: physical_factor_preparation.py (ordered AB, raw D and empty
complement); attempts/0029/reduced_active_basis.py (Q=Psi R, not raw-D slicing).
No admitted oracle, implementation moment helper or saved material output is
used as the independent answer. Missing-cross, wrong-D/metric, discarded empty
complement and false-adjoint mutations are evaluated numerical counterexamples.
"""
from __future__ import annotations

import itertools
import json

import numpy as np
import sympy as sp

from substrate_framework.material_noise_projection import projected_noise_decomposition
from substrate_framework.physical_factor_preparation import ActionFactors


PARTS = ('retained', 'omitted')
ATOL = 2e-12
RTOL = 2e-11


def annihilators(modes):
    """Literal occupation-bit CAR matrices, independent of Wick contraction."""
    out = []
    for mode in range(modes):
        matrix = np.zeros((2**modes, 2**modes), dtype=int)
        for bits in range(2**modes):
            if bits & (1 << mode):
                parity = (bits & ((1 << mode)-1)).bit_count()
                matrix[bits ^ (1 << mode), bits] = (-1)**parity
        out.append(matrix)
    return out


def symbolic_lift(operator, car):
    out = sp.zeros(car[0].rows)
    for i in range(operator.rows):
        for j in range(operator.cols):
            out += operator[i, j]*car[i].H*car[j]
    return out.applyfunc(sp.expand_complex)


def symbolic_fixture():
    """Complex nonorthogonal raw frame, coherent state and empty third mode."""
    i, rat = sp.I, sp.Rational
    unitary = sp.Matrix([[1/sp.sqrt(2), 0, 1/sp.sqrt(2)],
                         [i/sp.sqrt(2), 0, -i/sp.sqrt(2)], [0, 1, 0]])
    raw_map = sp.Matrix([[1, i/5], [rat(1, 7), 1]])
    prepared = sp.Matrix([[rat(2, 5), i/10], [-i/10, rat(1, 2)]])
    psi = (unitary[:, :2]*raw_map).applyfunc(sp.expand_complex)
    inverse = raw_map.inv().applyfunc(sp.expand_complex)
    density = (inverse*prepared*inverse.H).applyfunc(sp.expand_complex)
    coordinates = inverse[:, :1]
    operators = {
        'A': sp.Matrix([[1, 1+i, 2*i], [2-i, -1, i], [1, 3, 2]]),
        'B': sp.Matrix([[i, 2, 1-i], [-1, i, 2*i], [2+i, 1, 3]])}
    # Explicit two-orbital gauge-invariant Gaussian ensemble, third mode vacuum.
    # Its one-particle block is rho_2-det(rho_2) I, with independent vacuum and
    # double-occupancy probabilities. This does not use the covariance formula.
    fock = sp.zeros(8)
    double = prepared.det()
    fock[0, 0] = (sp.eye(2)-prepared).det()
    fock[3, 3] = double
    one_particle = prepared-double*sp.eye(2)
    for a in range(2):
        for b in range(2):
            fock[1 << a, 1 << b] = one_particle[a, b]
    return psi, density, coordinates, operators, unitary, prepared, fock


def symbolic_reference():
    """Exact PW versus actual Fock block identity, then public numerical check."""
    psi, density, coordinates, operators, unitary, _, fock = symbolic_fixture()
    # Canonicalize each exact block before composition, not one huge expression
    # after repeated matrix products. This preserves all exact Fock comparisons.
    gram = (psi.H*psi).applyfunc(sp.expand_complex)
    q = (psi*coordinates).applyfunc(sp.expand_complex)
    projector = (q*q.H).applyfunc(sp.expand_complex)
    rho = (psi*density*psi.H).applyfunc(sp.expand_complex)
    assert sp.simplify((rho*projector-projector*rho).norm()) != 0
    assert sp.simplify(sp.trace(fock)-1) == 0
    assert sp.simplify(coordinates.H*gram*coordinates) == sp.eye(1)
    car = [sp.Matrix(a) for a in annihilators(3)]
    pieces, lifted, means = {}, {}, {}
    for name, operator in operators.items():
        retained = (projector*operator*projector).applyfunc(sp.expand_complex)
        for part, value in (('full', operator), ('retained', retained),
                            ('omitted', operator-retained)):
            key = name, part
            pieces[key] = value
            lifted[key] = symbolic_lift((unitary.H*value*unitary).applyfunc(sp.expand_complex), car)
            means[key] = sp.simplify(sp.trace(fock*lifted[key]))
            assert sp.simplify(means[key]-sp.trace(rho*value)) == 0
    covariance = {}
    pairs = [(a, 'full', b, 'full') for a in operators for b in operators]
    pairs += [(a, pa, b, pb) for a in operators for b in operators
              for pa in PARTS for pb in PARTS]
    for a, pa, b, pb in pairs:
        left, right = (a, pa), (b, pb)
        exact = sp.simplify(sp.trace(fock*lifted[left]*lifted[right])
                            - means[left]*means[right])
        dense = sp.trace(rho*pieces[left]*(sp.eye(3)-rho)*pieces[right])
        assert sp.simplify(exact-dense) == 0
        covariance[a, pa, b, pb] = exact
    cross = [sp.simplify(covariance[a, 'retained', b, 'omitted']
                        + covariance[a, 'omitted', b, 'retained'])
             for a in operators for b in operators]
    assert any(value != 0 for value in cross)
    for a in operators:
        for b in operators:
            assert sp.simplify(covariance[a, 'full', b, 'full']-sum(
                covariance[a, pa, b, pb] for pa in PARTS for pb in PARTS)) == 0
    numeric = lambda matrix: np.array(matrix.evalf(), dtype=np.complex128)
    case = dict(psi=numeric(psi), gram=numeric(gram), d=numeric(density),
                coordinates=numeric(coordinates),
                operators={name: numeric(value) for name, value in operators.items()})
    reference = dict(
        full={(a, b): complex(value) for (a, pa, b, pb), value in covariance.items()
              if pa == pb == 'full'},
        blocks={key: complex(value) for key, value in covariance.items()
                if key[1] in PARTS and key[3] in PARTS},
        means={key: complex(value) for key, value in means.items()})
    error = compare(evaluate(case), reference)
    return {'exact_directed_pairs': len(covariance), 'public_interface_max_error': error}


def random_unitary(rng, dimension):
    value = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
    return np.linalg.qr(value)[0]


def numerical_fixture():
    """Independent complex raw coordinates; state and selected P do not commute."""
    rng = np.random.default_rng(2540030)
    unitary = random_unitary(rng, 5)
    rotation = random_unitary(rng, 3)
    selected = random_unitary(rng, 3)[:, :2]
    raw_map = np.array([[1.05, .07+.11j, 0], [.02j, .95, .09],
                        [.03, -.02j, 1.02]], dtype=complex)
    occupations = np.array([.15, .37, .61])
    inverse = np.linalg.solve(raw_map, np.eye(3))
    prepared = (rotation*occupations)@rotation.conj().T
    psi = unitary[:, :3]@raw_map
    operators = {}
    for name in ('A', 'B', 'H'):
        value = rng.normal(size=(5, 5))+1j*rng.normal(size=(5, 5))
        operators[name] = value+value.conj().T if name == 'H' else value
    return dict(psi=psi, gram=psi.conj().T@psi,
                d=inverse@prepared@inverse.conj().T,
                coordinates=inverse@selected, operators=operators,
                orthogonal_frame=unitary[:, :3], raw_map=raw_map,
                prepared=prepared, orbitals=unitary[:, :3]@rotation,
                occupations=occupations)


def evaluate(case, *, false_adjoint=False):
    psi = case['psi']
    probes = {}
    for name, operator in case['operators'].items():
        action = operator@psi
        adjoint = action if false_adjoint else operator.conj().T@psi
        probes[name] = ActionFactors(action, adjoint, 'finite:'+name,
            'finite:'+name+':independent-adjoint', 'fixture-energy',
            'finite-PW-basis', 'literal finite operator, not a material-source claim')
    return projected_noise_decomposition(psi, case['gram'], case['d'],
        case['coordinates'], probes, factor_identity='original-finite-frame',
        pw_basis_identity='finite-PW-basis',
        state_authority='declared finite quasifree state with empty complement')


def split_operators(case):
    q = case['psi']@case['coordinates']
    projector = q@q.conj().T
    out = {}
    for name, operator in case['operators'].items():
        retained = projector@operator@projector
        out[name, 'full'] = operator
        out[name, 'retained'] = retained
        out[name, 'omitted'] = operator-retained
    return out


def collect_reference(names, pieces, mean, covariance):
    return dict(
        full={(a, b): covariance((a, 'full'), (b, 'full')) for a in names for b in names},
        blocks={(a, pa, b, pb): covariance((a, pa), (b, pb))
                for a in names for b in names for pa in PARTS for pb in PARTS},
        means={key: mean(key) for key in pieces})


def dense_reference(case):
    pieces = split_operators(case)
    rho = case['psi']@case['d']@case['psi'].conj().T
    hole = np.eye(len(rho))-rho
    return collect_reference(case['operators'], pieces,
        lambda a: np.trace(rho@pieces[a]),
        lambda a, b: np.trace(rho@pieces[a]@hole@pieces[b]))


def fock_ensemble(orbitals, occupations, car):
    """Mixture of actual Slater vectors; unrepresented orbitals remain vacuum."""
    dimension = len(car[0])
    creators = [sum((orbitals[i, a]*car[i].conj().T for i in range(len(car))),
                    np.zeros((dimension, dimension), complex))
                for a in range(len(occupations))]
    ensemble = np.zeros((dimension, dimension), complex)
    for occupied in itertools.product((0, 1), repeat=len(occupations)):
        probability = np.prod([p if bit else 1-p for p, bit in zip(occupations, occupied)])
        vector = np.zeros(dimension, complex)
        vector[0] = 1
        for a in reversed(range(len(occupied))):
            if occupied[a]:
                vector = creators[a]@vector
        ensemble += probability*np.outer(vector, vector.conj())
    np.testing.assert_allclose(np.trace(ensemble), 1, atol=ATOL, rtol=RTOL)
    return ensemble


def fock_reference(case):
    pieces = split_operators(case)
    car = [a.astype(complex) for a in annihilators(len(case['psi']))]
    ensemble = fock_ensemble(case['orbitals'], case['occupations'], car)
    bilinears = {(i, j): a.conj().T@b for i, a in enumerate(car)
                for j, b in enumerate(car)}
    lifted = {key: sum((value[i, j]*term for (i, j), term in bilinears.items()),
                       np.zeros_like(ensemble)) for key, value in pieces.items()}
    means = {key: np.trace(ensemble@value) for key, value in lifted.items()}
    return collect_reference(case['operators'], pieces, lambda a: means[a],
        lambda a, b: np.trace(ensemble@lifted[a]@lifted[b])-means[a]*means[b])


def compare(actual, reference):
    """Compare every requested numerical quantity, not just reconstruction."""
    error = 0.
    for field in ('full', 'blocks', 'means'):
        assert actual[field].keys() == reference[field].keys(), field
        for key, expected in reference[field].items():
            np.testing.assert_allclose(actual[field][key], expected, atol=ATOL, rtol=RTOL,
                                       err_msg=str((field, key)))
            error = max(error, float(abs(actual[field][key]-expected)))
    if 'reconstruction_defect' in actual:
        for value in actual['reconstruction_defect'].values():
            np.testing.assert_allclose(value, 0, atol=ATOL, rtol=RTOL)
    return error


def changed_by(actual, mutated, field):
    difference = max(abs(actual[field][key]-mutated[field][key]) for key in actual[field])
    assert difference > 1e-5, ('insensitive mutation', field, difference)
    return float(difference)


def mutation_checks(case, actual):
    """Expose physically plausible omissions by evaluating altered inputs/results."""
    missing_cross = max(abs(actual['full'][a, b]
        - actual['blocks'][a, 'retained', b, 'retained']
        - actual['blocks'][a, 'omitted', b, 'omitted'])
        for a in case['operators'] for b in case['operators'])
    assert missing_cross > 1e-5
    wrong_d = evaluate(dict(case, d=case['prepared']))
    wrong_metric = evaluate(dict(case, gram=np.eye(3),
                                coordinates=case['raw_map']@case['coordinates']))
    span = case['orthogonal_frame']@case['orthogonal_frame'].conj().T
    truncated = {name: span@operator@span for name, operator in case['operators'].items()}
    lost_complement = evaluate(dict(case, operators=truncated))
    fake_adjoint = evaluate(case, false_adjoint=True)
    return dict(missing_cross=float(missing_cross),
                wrong_raw_D=changed_by(actual, wrong_d, 'full'),
                wrong_Gram_metric=changed_by(actual, wrong_metric, 'blocks'),
                discarded_initially_empty_complement=changed_by(actual, lost_complement, 'full'),
                false_nonhermitian_adjoint=changed_by(actual, fake_adjoint, 'full'))


def reframing_check(case, actual):
    change = np.array([[1, .2j, -.1], [.12, .9, .07j], [.03j, -.1, 1.1]], complex)
    inverse = np.linalg.solve(change, np.eye(3))
    transformed = dict(case, psi=case['psi']@change,
        gram=change.conj().T@case['gram']@change,
        d=inverse@case['d']@inverse.conj().T,
        coordinates=inverse@case['coordinates'])
    # The independent frame changes while the physical state and projector do not.
    return compare(evaluate(transformed), actual)


def empty_selection_check(case, actual):
    empty = evaluate(dict(case, coordinates=np.empty((3, 0), complex)))
    for a in case['operators']:
        for b in case['operators']:
            np.testing.assert_allclose(empty['blocks'][a, 'omitted', b, 'omitted'],
                                       actual['full'][a, b], atol=ATOL, rtol=RTOL)
            for pa, pb in itertools.product(PARTS, repeat=2):
                if 'retained' in (pa, pb):
                    np.testing.assert_allclose(empty['blocks'][a, pa, b, pb], 0,
                                               atol=ATOL, rtol=RTOL)
    return max(abs(value) for value in empty['reconstruction_defect'].values())


def main():
    symbolic = symbolic_reference()
    case = numerical_fixture()
    rho = case['psi']@case['d']@case['psi'].conj().T
    q = case['psi']@case['coordinates']
    projector = q@q.conj().T
    commutator = float(np.linalg.norm(rho@projector-projector@rho))
    assert commutator > 1e-3
    complement = np.eye(len(rho))-case['orthogonal_frame']@case['orthogonal_frame'].conj().T
    np.testing.assert_allclose(complement@rho, 0, atol=ATOL, rtol=RTOL)
    actual = evaluate(case)
    dense_error = compare(actual, dense_reference(case))
    fock_error = compare(actual, fock_reference(case))
    report = dict(symbolic=symbolic, numerical_dense_max_error=dense_error,
        numerical_Fock_max_error=fock_error, state_projector_commutator_frobenius=commutator,
        raw_frame_transport_max_error=reframing_check(case, actual),
        empty_selection_reconstruction_defect=float(empty_selection_check(case, actual)),
        mutation_separations=mutation_checks(case, actual),
        scope='finite t0 ordered electronic decomposition only; no material or late-time inference')
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
