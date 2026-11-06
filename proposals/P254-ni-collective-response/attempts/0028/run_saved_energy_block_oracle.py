#!/usr/bin/env python3
"""Independent exact complex trace and original-metadata action envelope."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
import numpy as np
import sympy as sp
from mixed_frame_response import DensityDyad
from mixed_material_moments import _one_body
from substrate_framework.streamed_material_state import authenticate, file_binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source_paths = [Path(__file__), HERE/'mixed_energy_block_work.py',
                    HERE/'mixed_material_moments.py', HERE/'mixed_frame_response.py',
                    HERE/'run_saved_energy_continuation.py']
    sources = [file_binding(path) for path in source_paths]
    prepared_binding = file_binding(args.prepared)
    prepared = json.loads(authenticate(prepared_binding).read_text())
    # Exact rectangular complex factors; no Hermiticity/conjugacy shortcut.
    left = sp.Matrix(4, 3, lambda i, j: sp.Rational(i+2*j+1, 7)+sp.I*sp.Rational(i-j, 11))
    right = sp.Matrix(4, 2, lambda i, j: sp.Rational(2*i-j+2, 9)+sp.I*sp.Rational(i+j+1, 13))
    coefficient = sp.Matrix([[1+sp.I, 2], [0, 0], [sp.I/3, sp.Rational(2, 5)]])
    h = sp.Matrix(4, 4, lambda i, j: sp.Rational(i+2*j+1, 17)+sp.I*sp.Rational(i-j+2, 19))
    weight = sp.Rational(2, 7)
    expected = sp.simplify(weight*sp.trace(h*left*coefficient*right.conjugate().T))
    convert = lambda value: np.asarray(value.evalf(), dtype=complex)
    dyad = DensityDyad(convert(left), convert(coefficient), convert(right), 0, float(weight))
    h_numeric = convert(h)
    errors = {}
    for block in (1, 2, 64):
        calls = []
        def action(vectors):
            calls.append(vectors.shape[1])
            return h_numeric@vectors
        actual = _one_body(dyad, action, block)
        if len(calls) != int(sp.ceiling(sp.Rational(left.cols, block))) or sum(calls) != left.cols:
            raise AssertionError('actual action/column count differs from exact block envelope')
        error = abs(actual-complex(expected))
        tolerance = 2e-12*max(1., abs(complex(expected)))
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError('complex rectangular block trace differs from exact SymPy')
        errors[str(block)] = dict(absolute_error=float(error), tolerance=float(tolerance),
                                 actual_h0_actions=str(len(calls)), actual_action_columns=str(sum(calls)))
    mutants = dict(missing_weight=expected/weight,
        missing_bra_conjugation=sp.simplify(weight*sp.trace(h*left*coefficient*right.T)),
        conjugated_coefficient=sp.simplify(weight*sp.trace(h*left*coefficient.conjugate()*right.conjugate().T)))
    if any(sp.simplify(value-expected) == 0 for value in mutants.values()):
        raise AssertionError('complex density trace mutation not exposed')
    columns = [slot['shape'][1] for slot in prepared['identity']['source']['slots']]
    if len(columns) != 128 or len(set(columns)) != 1:
        raise AssertionError('original uniform128-source raw-column coverage required')
    symbolic = sum(4*sp.ceiling(sp.Rational(width, 64)) for width in columns)
    direct = sum(4*((width+63)//64) for width in columns)
    cumulative = sp.Integer(1024)+symbolic-sum(4*sp.ceiling(sp.Rational(columns[i], 64)) for i in (0, 1, 2))
    direct_cumulative = 1024+sum(4*((width+63)//64) for i, width in enumerate(columns) if i not in (0, 1, 2))
    if symbolic != direct or cumulative != direct_cumulative or not cumulative <= 4096:
        raise AssertionError('original charged-cost/donor action envelope differs')
    action_mutants = dict(omitted_historical_cost=cumulative-1024,
        omitted_source=cumulative-4*sp.ceiling(sp.Rational(columns[-1], 64)),
        discarded_column_remainder=1024+sum(4*(width//64) for width in columns[3:]))
    if any(sp.simplify(value-cumulative) == 0 for value in action_mutants.values()):
        raise AssertionError('original charged-cost/source/remainder mutation not exposed')
    if sources != [file_binding(path) for path in source_paths] or prepared_binding != file_binding(args.prepared):
        raise ValueError('source/input changed during exact calculation')
    result = dict(schema='P254-saved-energy-block-oracle-v1', status='passed',
        source_bindings=sources, original_prepared=prepared_binding,
        exact_rectangular_complex_trace=str(expected), block_comparisons=errors,
        rejected_mutations=list(mutants)+list(action_mutants), original_source_count=str(len(columns)),
        original_columns_per_source=str(columns[0]), wide64_full_bare_action_upper_bound=str(symbolic),
        original_admitted_actions='1024', retained_completed_source_indices=['0', '1', '2'],
        cumulative_consumer_action_upper_bound_with_donors=str(cumulative),
        cumulative_action_ceiling='4096',
        scope='Exact algebra, grouping and original-metadata action bound; original diagnostic actions remain separately charged. No whole-run time/memory/physical-convergence claim.')
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
