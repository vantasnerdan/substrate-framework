"""Exact complex perturbation projection and orthonormal band-basis covariance."""
import sympy as sp

from oracle_evidence import require


def zero(matrix):
    return all(sp.simplify(value) == 0 for value in matrix)


def combine(coefficients, values):
    return sum((coefficient * value for coefficient, value in zip(coefficients, values)), sp.zeros(2))


def run(ledger):
    root = sp.sqrt(2)
    v12 = sp.Matrix([[1, sp.I, 0], [sp.I, 1, 0], [0, 0, root]]) / root
    v23 = sp.Matrix([[root, 0, 0], [0, 1, sp.I], [0, sp.I, 1]]) / root
    patterns = v12 * v23
    require(ledger, 'projection.unitary_complex_pattern_basis', zero(patterns.conjugate().T * patterns - sp.eye(3)))
    target = sp.Matrix([1, 0, 0])
    coefficients = patterns.conjugate().T * target
    require(ledger, 'projection.solved_cartesian_coefficients', zero(patterns * coefficients - target))
    cartesian = [sp.Matrix([[1, sp.I / 3], [-sp.I / 3, 2]]),
                 sp.Matrix([[-2, (2 + 3 * sp.I) / 7], [(2 - 3 * sp.I) / 7, 4]]),
                 sp.Matrix([[3, -sp.Rational(1, 5)], [-sp.Rational(1, 5), -1]])]
    shifts = [sp.Rational(2, 7), -sp.Rational(3, 11), sp.Rational(5, 13)]
    patterned_absolute = [combine(patterns[:, j], cartesian) for j in range(3)]
    patterned_shift = [sum(patterns[i, j] * shifts[i] for i in range(3)) for j in range(3)]
    patterned_native = [matrix - shift * sp.eye(2) for matrix, shift in zip(patterned_absolute, patterned_shift)]
    projected_native = combine(coefficients, patterned_native)
    projected_shift = sp.simplify(sum(c * value for c, value in zip(coefficients, patterned_shift)))
    projected_absolute = combine(coefficients, patterned_absolute)
    require(ledger, 'projection.cartesian_absolute_recovered', zero(projected_absolute - cartesian[0]))
    require(ledger, 'projection.cartesian_native_shift_recovered', sp.simplify(projected_shift - shifts[0]) == 0)
    require(ledger, 'projection.reconstruct_before_or_after_projection',
            zero(projected_native + projected_shift * sp.eye(2) - projected_absolute))

    wrong_conjugate_shift = sum(sp.conjugate(c) * value for c, value in zip(coefficients, patterned_shift))
    missing_mode_shift = sum(coefficients[i] * patterned_shift[i] for i in range(2))
    require(ledger, 'projection.exposes_extra_conjugation', sp.simplify(wrong_conjugate_shift - projected_shift) != 0)
    require(ledger, 'projection.exposes_missing_shift_mode', sp.simplify(missing_mode_shift - projected_shift) != 0)
    require(ledger, 'projection.exposes_conjugated_operator_coefficients',
            not zero(combine([sp.conjugate(c) for c in coefficients], patterned_absolute) - projected_absolute))

    unitary = sp.Matrix([[1, sp.I], [sp.I, 1]]) / root
    require(ledger, 'projection.unitary_band_basis', zero(unitary.conjugate().T * unitary - sp.eye(2)))
    transformed_native = unitary.conjugate().T * projected_native * unitary
    transformed_absolute = unitary.conjugate().T * projected_absolute * unitary
    require(ledger, 'projection.scalar_correction_unitary_covariance',
            zero(transformed_native + projected_shift * sp.eye(2) - transformed_absolute))
    nonunitary = sp.diag(2, 1)
    wrong_metric = nonunitary.T * projected_native * nonunitary + projected_shift * sp.eye(2)
    require(ledger, 'projection.nonorthogonal_basis_counterexample',
            not zero(wrong_metric - nonunitary.T * projected_absolute * nonunitary))
    return {'cartesian_shift': str(projected_shift), 'coefficients': [str(sp.simplify(c)) for c in coefficients],
            'mutants_refused': ['extra_conjugation', 'missing_shift_mode', 'conjugated_operator_coefficients'],
            'domain': 'Ordinary NC orthonormal same-k/same-spin basis; nonorthogonal overlap metrics remain outside the contract.'}
