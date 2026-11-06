"""Exact operator identities and isolated two-level energy derivatives."""
import sympy as sp

from oracle_evidence import require


def zero(matrix):
    return all(sp.simplify(value) == 0 for value in matrix)


def isolated_slope(matrix, parameter, eigenvalue, energy):
    polynomial = (matrix - energy * sp.eye(2)).det()
    numerator = -sp.diff(polynomial, parameter).subs({parameter: 0, energy: eigenvalue})
    denominator = sp.diff(polynomial, energy).subs({parameter: 0, energy: eigenvalue})
    return sp.cancel(numerator / denominator)


def run(ledger):
    a, b, zr, zi, shift, mu0, t, energy = sp.symbols('a b zr zi shift mu0 t energy', real=True)
    gap = sp.Symbol('gap', positive=True)
    absolute = sp.Matrix([[a, zr + sp.I * zi], [zr - sp.I * zi, b]])
    native = absolute - shift * sp.eye(2)
    require(ledger, 'algebra.reconstruction_all_entries', zero(native + shift * sp.eye(2) - absolute))
    require(ledger, 'algebra.offdiagonals_not_shifted',
            native[0, 1] == absolute[0, 1] and native[1, 0] == absolute[1, 0])

    h_absolute = sp.diag(0, gap) + t * absolute
    h_native = h_absolute - (mu0 + t * shift) * sp.eye(2)
    slopes = []
    for index, root in enumerate((0, gap)):
        raw_slope = isolated_slope(h_absolute, t, root, energy)
        relative_slope = isolated_slope(h_native, t, root - mu0, energy)
        require(ledger, f'algebra.isolated_absolute_branch_{index}', sp.simplify(raw_slope - absolute[index, index]) == 0)
        require(ledger, f'algebra.isolated_relative_branch_{index}', sp.simplify(relative_slope - native[index, index]) == 0)
        require(ledger, f'algebra.isolated_reconstruction_branch_{index}', sp.simplify(relative_slope + shift - raw_slope) == 0)
        slopes.append({'absolute': str(raw_slope), 'relative': str(relative_slope)})

    concrete_shift = sp.Rational(2, 7)
    witness = sp.Matrix([[1, (2 + 3 * sp.I) / 5], [(2 - 3 * sp.I) / 5, -2]])
    relative = witness - concrete_shift * sp.eye(2)
    mutants = {
        'omitted': relative,
        'opposite_sign': relative - concrete_shift * sp.eye(2),
        'doubled': relative + 2 * concrete_shift * sp.eye(2),
        'offdiagonal_shift': relative + concrete_shift * sp.ones(2),
        'rounded_or_fitted_substitute': relative + sp.Rational(3, 10) * sp.eye(2),
    }
    for name, candidate in mutants.items():
        require(ledger, 'algebra.exposes_' + name, not zero(candidate - witness))

    degenerate_perturbation = sp.Matrix([[0, 1], [1, 0]])
    require(ledger, 'algebra.degenerate_single_band_counterexample',
            set(degenerate_perturbation.eigenvals()) == {-1, 1}
            and degenerate_perturbation[0, 0] == 0)
    native_choices = [-sp.Rational(4, 5), sp.Rational(3, 4)]
    absolute_choices = [value + sp.Rational(1, 10) for value in native_choices]
    require(ledger, 'algebra.absolute_gauge_requires_fresh_selection',
            max(range(2), key=lambda i: abs(native_choices[i])) !=
            max(range(2), key=lambda i: abs(absolute_choices[i])))
    return {'isolated_branch_slopes': slopes, 'mutants_refused': list(mutants),
            'domain': 'Isolated slopes use a strictly positive gap; degeneracy counterexample is not licensed.',
            'scope': 'Exact finite algebra, not a physical material convergence certificate.'}
