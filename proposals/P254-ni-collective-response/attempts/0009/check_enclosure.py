"""Bounded correction of attempt0008's force-print enclosure verifier."""
import numpy as np
import sympy as sp
from scipy.linalg import norm

from substrate_framework.force_print_precision import (
    barycentric_box_identity, component_box_max_squared, corner_squared_norms,
)
from substrate_framework.material_primary_boundary import force_rounding_upper


def check(ledger, name, condition):
    ledger.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def theory(ledger):
    identity = barycentric_box_identity()
    check(ledger, 'box interpolation weights sum to one', identity['weight_sum'] == 1)
    check(ledger, 'box interpolation reproduces all three coordinates',
          identity['position_gaps'] == (0, 0, 0))
    variance = identity['variance']
    check(ledger, 'exact squared-norm interpolation variance identity',
          sp.expand(variance - identity['expected_variance']) == 0)
    t = identity['parameters'][0]
    check(ledger, 't(1-t) cannot be negative on the declared closed unit interval',
          sp.reduce_inequalities([t >= 0, t <= 1, t * (1 - t) < 0], t) == sp.false)
    check(ledger, 'radius squared is nonnegative',
          (identity['radius'] ** 2).is_nonnegative is True)
    midpoint = {value: sp.Rational(1, 2) for value in identity['parameters']}
    midpoint[identity['radius']] = 1
    check(ledger, 'omitted variance term exposed by an exact interior witness',
          variance != 0 and variance.subs(midpoint) == 3)
    return {'weight_sum': str(identity['weight_sum']),
            'variance': str(sp.factor(variance)),
            'domain': 'h>=0 and each0<=t_i<=1; nonnegative product weights and variance bound every interior point by corner squared norms',
            'omitted_variance_witness': 'h=1,t_i=1/2 gives exact gap3'}


def boxes(ledger, rows, radius, label):
    result = []
    for index, row in enumerate(rows):
        exact = component_box_max_squared(row, radius)
        enumerated = max(corner_squared_norms(row, radius))
        check(ledger, f'{label} row{index}: independent eight-corner maximum', exact == enumerated)
        result.append(exact)
    return result


def encloses(candidate, rows, radius, exact_maxima, slack):
    upper = np.asarray(candidate(rows, radius), float)
    if upper.shape != (len(rows),) or not np.isfinite(upper).all() or (upper < 0).any():
        return False
    allowance = sp.Rational(str(slack))
    return all((sp.Rational(str(bound)) + allowance) ** 2 >= exact
               for bound, exact in zip(upper, exact_maxima))


def native(ledger, rows, radius, gate, original_maximum):
    exact = boxes(ledger, rows, radius, 'native')
    slack = 64 * np.finfo(float).eps * max(1, float(np.max(np.abs(rows))))
    check(ledger, 'actual triangle bounds enclose every native component box',
          encloses(force_rounding_upper, rows, radius, exact, slack))
    upper = force_rounding_upper(rows, radius)
    check(ledger, 'native maximum exactly preserves original receipt',
          float(max(upper)) == original_maximum)
    check(ledger, 'source-print enclosure still satisfies frozen primary gate',
          sp.Rational(str(max(upper))) < sp.Rational(str(gate)))
    return {'rows': len(rows), 'exact_box_squared_maxima': list(map(str, exact)),
            'max_triangle_upper_ev_per_angstrom': float(max(upper)),
            'gate_ev_per_angstrom': gate, 'arithmetic_allowance_ev_per_angstrom': slack,
            'precision_scope': 'Exact decimal embeddings of the actual float function inputs; physical CODATA conversion and DFT uncertainty retain original separate budgets.'}


def mutations(ledger, radius, gate):
    h, g = sp.Rational(str(radius)), sp.Rational(str(gate))
    rows = np.array([[0, 0, 0], [-float(h), 2 * float(h), -3 * float(h)],
                     [float(g - h / 2), 0, 0]], float)
    exact = boxes(ledger, rows, radius, 'witness')
    slack = 64 * np.finfo(float).eps
    check(ledger, 'correct function encloses zero, mixed-sign and near-gate boxes',
          encloses(force_rounding_upper, rows, radius, exact, slack))
    center_squared = sum(sp.Rational(str(value)) ** 2 for value in rows[-1])
    check(ledger, 'near-gate printed center is below gate but an allowed corner exceeds it',
          center_squared < g ** 2 < exact[-1])
    mutants = {
        'omitted_or_zero_radius': lambda f, r: np.array([norm(row) for row in f]),
        'halved_radius': lambda f, r: force_rounding_upper(f, r / 2),
        'missing_sqrt3': lambda f, r: np.array([norm(row) for row in f]) + r,
    }
    refusals = {}
    for name, candidate in mutants.items():
        refusals[name] = not encloses(candidate, rows, radius, exact, slack)
        check(ledger, name + ' fails independent exact component-box enclosure', refusals[name])
    correct_accepts = bool(force_rounding_upper(rows[-1:], radius)[0] <= gate)
    omitted_accepts = bool(mutants['omitted_or_zero_radius'](rows[-1:], radius)[0] <= gate)
    check(ledger, 'near-threshold omitted radius produces an exposed false acceptance',
          not correct_accepts and omitted_accepts)
    return {'witness_rows': rows.tolist(), 'exact_box_squared_maxima': list(map(str, exact)),
            'mutant_refusals': refusals,
            'equivalent_mutations': 'Omission and zeroing are the same mathematical function and form one nonduplicated class.',
            'near_threshold': {'printed_squared_norm': str(center_squared),
                               'exact_max_squared_norm': str(exact[-1]),
                               'gate_squared': str(g ** 2),
                               'correct_accepts': correct_accepts, 'omitted_radius_accepts': omitted_accepts}}
