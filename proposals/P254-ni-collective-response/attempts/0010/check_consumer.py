"""Exercise the frozen production gauge API against independent exact fixtures."""
import importlib.util
import numpy as np
import sympy as sp

from oracle_evidence import require


def array(matrix):
    return np.array([[complex(sp.N(value, 30)) for value in row] for row in matrix.tolist()], dtype=complex)


def run(ledger, consumer_path):
    spec = importlib.util.spec_from_file_location('p254_frozen_gauge_consumer', consumer_path)
    consumer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(consumer)
    shifts = [sp.Rational(2, 7), sp.Integer(0), -sp.Rational(3, 11) + sp.I / 13]
    header = {'fermi_shift_applied': True, 'fermi_shift_ry_per_bohr': [complex(sp.N(value, 30)) for value in shifts]}
    fixtures, expected_matrices, tolerances, errors = [], [], [], []
    for index, shift in enumerate(shifts):
        number = sp.Integer(index + 1)
        absolute = sp.Matrix([[number / 3, sp.Rational(2, 7) + sp.I / 5, -sp.Rational(1, 11)],
                              [sp.Rational(2, 7) - sp.I / 5, -number / 2, sp.Rational(5, 13) + sp.I / 17],
                              [-sp.Rational(1, 11), sp.Rational(5, 13) - sp.I / 17, number + sp.Rational(1, 19)]])
        native = absolute - shift * sp.eye(3)
        record = {'kind': 'SCREENED_NATIVE', 'perturbation': index + 1, 'matrix': array(native)}
        expected = array(absolute)
        tolerance = 64 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(expected))))
        first = consumer.screened_hamiltonian_matrix(header, record)
        first_correct = np.allclose(first, expected, rtol=0, atol=tolerance)
        repeated = consumer.screened_hamiltonian_matrix(header, record)
        require(ledger, f'consumer.absolute_local_perturbation_{index + 1}_repeated_query',
                first_correct and np.allclose(repeated, expected, rtol=0, atol=tolerance))
        require(ledger, f'consumer.diagonal_local_perturbation_{index + 1}',
                np.allclose(consumer.screened_hamiltonian_diagonal(header, record), np.diag(expected), rtol=0, atol=tolerance))
        fixtures.append(record)
        expected_matrices.append(expected)
        tolerances.append(tolerance)
        errors.append(float(np.max(np.abs(repeated - expected))))

    no_shift_header = {'fermi_shift_applied': False, 'fermi_shift_ry_per_bohr': [0j, 0j, 0j]}
    no_shift_record = fixtures[1]
    require(ledger, 'consumer.branch_false_has_no_scalar_correction',
            np.array_equal(consumer.screened_hamiltonian_matrix(no_shift_header, no_shift_record), expected_matrices[1])
            and np.array_equal(consumer.screened_hamiltonian_diagonal(no_shift_header, no_shift_record), np.diag(expected_matrices[1])))
    invalid_header = {'fermi_shift_applied': False, 'fermi_shift_ry_per_bohr': [1j, 0j, 0j]}
    try:
        consumer.screened_hamiltonian_matrix(invalid_header, fixtures[0])
    except ValueError:
        branch_refused = True
    else:
        branch_refused = False
    require(ledger, 'consumer.refuses_nonzero_correction_outside_applied_branch', branch_refused)
    bare = dict(fixtures[0], kind='BARE')
    try:
        consumer.screened_hamiltonian_matrix(header, bare)
    except ValueError:
        bare_refused = True
    else:
        bare_refused = False
    require(ledger, 'consumer.refuses_screened_gauge_on_bare_operator', bare_refused)
    return {'cases': ['first local perturbation', 'applied branch with zero shift', 'last local perturbation with complex shift',
                      'branch false with zero removed-shift coefficient', 'branch false with nonzero coefficient refused', 'bare misuse refused'],
            'arithmetic_allowance': '64*float64_epsilon*max(1,max_abs(exact_expected_matrix)); not a physical convergence tolerance',
            'absolute_matrix_max_errors': errors, 'absolute_arithmetic_allowances': tolerances,
            'scope': 'Production consumer is executed from its immutable exact-byte archive, not a duplicated implementation.'}
