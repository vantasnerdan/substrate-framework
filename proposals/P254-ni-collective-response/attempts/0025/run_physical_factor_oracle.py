"""New factor/Fock/Fourier oracle; no material ensemble or old oracle replay.

Run from repository root with PYTHONPATH=src. All numeric comparisons use
512*eps*max(1, number of scalar entries)*(1+max(abs(reference))). No fitting,
Hermitization, clipping, material execution, or timing claim occurs here.
"""
from itertools import combinations
import json
import numpy as np
import sympy as sp
from scipy.linalg import expm, null_space

from substrate_framework.physical_factor_preparation import (
    ActionFactors, FactorFunctional, FactorPreparation,
    coefficients_from_orthonormal_populations,
)
from substrate_framework.material_observable_factors import (
    NumberProbeGrid, gamma_probe_actions, CompleteGridNormalization,
    gamma_correlations,
)


def second_quantized(matrix):
    """Literal occupation-bit CAR, independent of factor implementation."""
    n = len(matrix)
    result = np.zeros((2**n, 2**n), complex)
    for ket in range(2**n):
        for j in range(n):
            if not (ket >> j) & 1:
                continue
            removed = ket ^ (1 << j)
            sj = (-1)**((ket & ((1 << j)-1)).bit_count())
            for i in range(n):
                if (removed >> i) & 1:
                    continue
                si = (-1)**((removed & ((1 << i)-1)).bit_count())
                result[removed | (1 << i), ket] += si*sj*matrix[i, j]
    return result


def symbolic_block():
    """Exact nonorthogonal coefficients and empty-mode Fock trace."""
    p, q = sp.symbols('p q', real=True)
    r = sp.Matrix([[1, 1+sp.I], [0, 2]])
    psi = sp.Matrix([[1, 1+sp.I], [0, 2], [0, 0]])
    d = r.inv()*sp.diag(p, q)*r.inv().H
    a = sp.Matrix([[1, sp.I, 2], [-sp.I, 0, 1], [2, 1, -1]])
    b = sp.Matrix([[0, 2, sp.I], [1, 1, 3], [2*sp.I, 0, 1]])
    # Build exact CAR without passing symbolic entries through floating point.
    car = []
    for j in range(3):
        c = sp.zeros(8)
        for ket in range(8):
            if ket & (1 << j):
                c[ket ^ (1 << j), ket] = (-1)**((ket & ((1 << j)-1)).bit_count())
        car.append(c)
    def lift(m):
        return sum((m[i,j]*car[i].H*car[j] for i in range(3)
                    for j in range(3)), sp.zeros(8))
    weights = [((p if k & 1 else 1-p)*(q if k & 2 else 1-q)
                if not k & 4 else 0) for k in range(8)]
    state = sp.diag(*weights)
    ah, bh = lift(a), lift(b)
    exact = sp.trace(state*ah*bh)-sp.trace(state*ah)*sp.trace(state*bh)
    ma, mb = psi.H*a*psi, psi.H*b*psi
    tab = (a.H*psi).H*(b*psi)
    candidate = sp.trace(d*tab)-sp.trace(d*ma*d*mb)
    assert sp.simplify(candidate-exact) == 0
    assert sp.simplify(sp.trace(d*(psi.H*psi))-(p+q)) == 0
    complement_lost = sp.trace(d*ma*(psi.H*psi).inv()*mb)-sp.trace(d*ma*d*mb)
    assert sp.simplify(candidate-complement_lost) != 0
    return {'exact_connected_identity': str(sp.factor(exact)),
            'empty_complement_difference': str(sp.factor(candidate-complement_lost))}


def fock_state(unitary, populations):
    """Exterior minors construct the full density independently of Wick."""
    n = len(populations)
    state = np.zeros((2**n, 2**n), complex)
    for occupied_bits in range(2**n):
        occupied = [i for i in range(n) if occupied_bits & (1 << i)]
        probability = np.prod([populations[i] if i in occupied else 1-populations[i]
                               for i in range(n)])
        vector = np.zeros(2**n, complex)
        for rows in combinations(range(n), len(occupied)):
            vector[sum(1 << i for i in rows)] = np.linalg.det(unitary[np.ix_(rows, occupied)])
        state += probability*np.outer(vector, vector.conj())
    return state


def close(actual, reference):
    actual, reference = np.asarray(actual), np.asarray(reference)
    error = float(np.max(np.abs(actual-reference), initial=0))
    allowance = 512*np.finfo(float).eps*max(1, actual.size)*(1+float(np.max(np.abs(reference), initial=0)))
    assert error <= allowance, (error, allowance)
    return error


def different(actual, mutation):
    scale = 1+float(np.max(np.abs(actual)))
    assert np.max(np.abs(np.asarray(actual)-mutation)) > 1e-8*scale


def functional(psi, operators):
    probes = {name: ActionFactors(a@psi, a.conj().T@psi, name, name+'^dagger',
                                 'oracle', 'oracle-PW', 'literal finite operator')
              for name, a in operators.items()}
    return FactorFunctional.from_actions(psi, probes, factor_identity='oracle factors',
                                         pw_basis_identity='oracle-PW')


def finite_fock_block():
    seed = np.array([[0, 1+2j, 2, -1j], [-1+2j, 0, 1j, 1],
                     [-2, 1j, 0, 2j], [-1j, -1, 2j, 0]], complex)
    u = expm(seed/5)
    r = np.array([[1.4, .3+.6j], [0, .8]], complex)
    f = np.array([.23, .71])
    psi = u[:, :2]@r
    rinv = np.linalg.inv(r)
    d = rinv@np.diag(f)@rinv.conj().T
    a = np.array([[1, 1j, 2, 0], [-1j, -2, .3, 1j], [2, .3, 1, 2], [0, -1j, 2, 0]])
    b = np.array([[0, 2j, 1, 1], [1, 2, -1j, 0], [2j, 1, 0, 3], [0, 1j, 1, -1]])
    operators = {'A': a, 'B': b, 'Bd': b.conj().T, 'diagonal': u@np.diag([2, -1, 4, 4])@u.conj().T,
                 'degenerate': u@np.array([[0, 1j, 0, 0], [-1j, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])@u.conj().T}
    h0 = u@np.diag([1, 1, 3, 4])@u.conj().T
    close(h0@operators['degenerate']-operators['degenerate']@h0, np.zeros((4,4)))
    kernel = functional(psi, operators)
    prep = FactorPreparation(kernel, d, 'synthetic oracle independent quasifree state')
    state = fock_state(u, np.r_[f, 0, 0])
    close(np.trace(state), 1)
    errors = []
    for left, al in operators.items():
        ah = second_quantized(al)
        errors.append(close(prep.mean(left), np.trace(state@ah)))
        for right, br in operators.items():
            bh = second_quantized(br)
            exact = np.trace(state@ah@bh)-np.trace(state@ah)*np.trace(state@bh)
            errors.append(close(prep.ordered(left, right), exact))
    positive = prep.ordered('Bd', 'B')
    assert abs(positive.imag) < 1e-11 and positive.real >= 0
    all_f = np.r_[f, 0, 0]
    b_natural = u.conj().T@b@u
    positive_sum = np.sum((1-all_f)[:,None]*abs(b_natural)**2*all_f[None,:])
    errors.append(close(positive, positive_sum))
    close(prep.ordered('diagonal', 'diagonal'), np.sum(f*(1-f)*np.array([4, 1])))
    close(prep.ordered('degenerate', 'degenerate'), f[0]*(1-f[1])+f[1]*(1-f[0]))
    m, s = kernel.overlaps, kernel.gram
    dropped = np.trace(d@m['A']@np.linalg.solve(s, m['B']))-np.trace(d@m['A']@d@m['B'])
    different(prep.ordered('A', 'B'), dropped)
    different(prep.ordered('A', 'B'), prep.ordered('A', 'Bd'))
    different(prep.ordered('A', 'B'), prep.ordered('B', 'A'))
    wrong_t = (a.T@psi).conj().T@(b@psi)
    wrong_conjugation = np.trace(d@wrong_t)-np.trace(d@m['A']@d@m['B'])
    different(prep.ordered('A', 'B'), wrong_conjugation)
    # Arbitrary invertible factor-coordinate change and physical unitary basis change.
    v = np.array([[1, .2j], [.1, 1.2]], complex)
    vinv = np.linalg.inv(v)
    rotated = FactorPreparation(functional(psi@v, operators), vinv@d@vinv.conj().T, 'same state coordinates')
    basis = expm(seed/7)
    changed = FactorPreparation(functional(basis@psi, {k: basis@x@basis.conj().T for k,x in operators.items()}), d, 'same physical state')
    for other in (rotated, changed):
        errors.append(close(other.ordered('A', 'B'), prep.ordered('A', 'B')))
    # Explicit Gram mapping means populations in the symmetric orthonormal frame.
    mapped = coefficients_from_orthonormal_populations(s, f, mapping_authority='synthetic symmetric frame')
    values, vectors = np.linalg.eigh(s)
    frame = psi@((vectors/np.sqrt(values))@vectors.conj().T)
    expected = frame@np.diag(f)@frame.conj().T
    errors.append(close(psi@mapped@psi.conj().T, expected))
    FactorPreparation(kernel, mapped, 'synthetic mapped state')
    for bad in (np.eye(2)*10, -np.eye(2), d+np.array([[0, 1j], [0, 0]])):
        try:
            FactorPreparation(kernel, bad, 'invalid oracle state')
        except ValueError:
            pass
        else:
            raise AssertionError('inadmissible coefficient accepted')
    return {'max_error': max(errors), 'positive_Bdagger_B': positive.real,
            'nonorthogonal_gram_defect': float(np.max(np.abs(s-np.eye(2))))}


def literal_fourier(shape, miller, kind):
    """Independent real-grid sum, not the FFT action or its index helper."""
    matrix = np.zeros((miller.shape[1], miller.shape[1]), complex)
    for index in np.ndindex(shape):
        fraction = np.array(index)/shape
        field = (np.cos if kind == 'cos' else np.sin)(2*np.pi*fraction[0])
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                matrix[i,j] += field*np.exp(2j*np.pi*np.dot(miller[:,j]-miller[:,i], fraction))/np.prod(shape)
    return matrix


def material_probe_block():
    shape = (7, 5, 3)
    miller = np.array([[0, 1, -1, 2], [0, 0, 0, 1], [0, 0, 0, 0]])
    cell = np.diag([4., 5., 6.])
    grid = NumberProbeGrid(shape, miller, cell, 2*np.pi*np.linalg.inv(cell).T, 'oracle-PW')
    psi = np.array([[1, .2j], [.1j, .8], [.3, .2], [0, .1j]], complex)
    g = np.array([[1, 1j, 2, 0], [-1j, 2, 0, 1], [2, 0, 0, 0], [0, 1, 0, 1]])
    preparations = {}
    sector_expected = {}
    errors = []
    for spin in (1, 2):
        probes = gamma_probe_actions(psi, spin_one_based=spin, pw_basis_identity='oracle-PW',
            source_authority='synthetic saved action', gpsi_h32_x=g@psi,
            gpsi_ni1_x=2*g@psi, number_grid=grid)
        for kind in ('cos', 'sin'):
            matrix = literal_fourier(shape, miller, kind)
            errors.append(close(probes['number_'+kind+'_b1'].action, matrix@psi))
            errors.append(close(probes['number_'+kind+'_b1'].adjoint_action, matrix.conj().T@psi))
        errors.append(close(probes['force_H32_x'].action, -g@psi))
        kernel = FactorFunctional.from_actions(psi, probes, factor_identity='synthetic spin '+str(spin), pw_basis_identity='oracle-PW')
        d = coefficients_from_orthonormal_populations(kernel.gram, [.15*spin, .25*spin], mapping_authority='oracle only')
        preparations[spin] = FactorPreparation(kernel, d, 'independent synthetic spin '+str(spin))
        arrays, metadata = kernel.archive_payload()
        assert metadata['state'] is None and metadata['density_coefficients'] is None
        close(arrays['gram'], psi.conj().T@psi)
        for entry in metadata['array_layout']['ordered_products']:
            close(arrays[entry['array']], kernel.ordered_products[entry['left'], entry['right']])
        values, vectors = np.linalg.eigh(kernel.gram)
        frame = psi@((vectors/np.sqrt(values))@vectors.conj().T)
        complete_frame = np.column_stack((frame, null_space(frame.conj().T)))
        state = fock_state(complete_frame, np.array([.15*spin, .25*spin, 0, 0]))
        force_fock = second_quantized(-g)
        mean = np.trace(state@force_fock)
        sector_expected[spin] = np.trace(state@force_fock@force_fock)-mean**2
        errors.append(close(preparations[spin].mean('force_H32_x'), mean))
        different(mean, -mean)
        for kind in ('cos', 'sin'):
            number_fock = second_quantized(literal_fourier(shape, miller, kind))
            expected = np.trace(state@force_fock@number_fock)-mean*np.trace(state@number_fock)
            errors.append(close(preparations[spin].ordered('force_H32_x', 'number_'+kind+'_b1'), expected))
    norm = CompleteGridNormalization((4,4,4), 'chosen oracle collective representation')
    result = gamma_correlations(preparations, 'force_H32_x', 'force_H32_x', norm, independence_authority='synthetic spin product')
    raw = sum(sector_expected.values())
    errors.append(close(result['electronic']['ordered_AB'], raw/64))
    different(raw/64, raw/64**2)
    different(raw/64, 2*sector_expected[1]/64)
    assert result['total'] is None and result['ionic'] is None and result['cross'] is None
    assert result['complete'] is False
    return {'max_error': max(errors), 'normalization': norm.metadata()}


def main():
    symbolic = symbolic_block()
    finite = finite_fock_block()
    material = material_probe_block()
    print(json.dumps({'scope': 'new finite algebra and literal Fourier oracle only',
                      'symbolic': symbolic, 'finite_fock': finite,
                      'material_adapter': material}, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
