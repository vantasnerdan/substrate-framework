#!/usr/bin/env python3
"""Check actual matched-source final-RK costs, without inventing an ETA or cause."""
import argparse
import json
import hashlib
from pathlib import Path

import sympy as sp
import run_stage_sensitive_completion_assessment as saved


def original(binding):
    path = Path(binding['path'])
    raw = path.read_bytes()
    if len(raw) != int(binding['bytes']) or hashlib.sha256(raw).hexdigest() != binding['sha256']:
        raise ValueError('original cost evidence bytes differ')
    return json.loads(raw)


def refutes_proportional(previous, measured):
    return bool(measured > previous*sp.Rational(8, 7))


def observe(source):
    rows = []
    for expected, pair in enumerate(source['pairs']):
        prior, actual = original(pair['prior']), original(pair['actual'])
        if (prior['index'] != expected or actual['index'] != expected
                or prior['stage'] != 2 or actual['stage'] != 3
                or prior['complete'] is not True or actual['complete'] is not True):
            raise ValueError('matched complete source/stage evidence differs')
        left, right = prior['observed_execution'], actual['observed_execution']
        cpu = 6+expected % 5
        if (left['cpu'] != cpu or right['cpu'] != cpu
                or left['affinity'] != [cpu] or right['affinity'] != [cpu]
                or left['hamiltonian_actions'] != 7 or right['hamiltonian_actions'] != 8):
            raise ValueError('original matched source/action/affinity differs')
        if sp.Integer(left['action_columns'])*8 != sp.Integer(right['action_columns'])*7:
            raise ValueError('original full-column action scaling differs')
        previous = saved.rational(left['elapsed_seconds'])
        measured = saved.rational(right['elapsed_seconds'])
        predicted = previous*sp.Rational(8, 7)
        rows.append(dict(index=expected, cpu=cpu, prior=pair['prior'], actual=pair['actual'],
            prior_seven_action_elapsed_seconds=str(previous),
            actual_eight_action_elapsed_seconds=str(measured),
            proportional_eight_over_seven_elapsed_seconds=str(predicted),
            actual_to_proportional_ratio=str(measured/predicted),
            actual_exceeds_proportional=refutes_proportional(previous, measured),
            original_prior_measurements=left, original_actual_measurements=right))
    if len(rows) != 5:
        raise ValueError('one actual complete original item per allocated CPU required')
    observed_verdict = all(row['actual_exceeds_proportional'] for row in rows)
    # Replace actual observations with the previously assumed proportional
    # durations: the same consumer must no longer reject that counterfactual.
    counterfactual_verdict = all(refutes_proportional(
        saved.rational(row['prior_seven_action_elapsed_seconds']),
        saved.rational(row['proportional_eight_over_seven_elapsed_seconds'])) for row in rows)
    if not observed_verdict or observed_verdict == counterfactual_verdict:
        raise AssertionError('actual runtime counterexample/mutation is not sensitive')
    return dict(schema='P254-actual-matched-final-RK-cost-observation-v1', status='PASS',
        samples=rows, checks=[
            dict(name='matched_original_stages_sources_actions_and_affinities', passed=True),
            dict(name='exact_original_full_column_eight_over_seven_scaling', passed=True),
            dict(name='actual_elapsed_refutes_proportional_callback_time_model', passed=observed_verdict),
            dict(name='counterfactual_proportional_observations_change_verdict', passed=True)],
        scope='First five actual final-source transitions only; no full physical endpoint, runtime bound, ETA, optimization cause or energy success inferred. Action-count proportionality is not elapsed-time proportionality.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, plan = json.loads(args.source.read_bytes()), json.loads(args.plan.read_bytes())
    if plan['source'] != saved.binding(args.source):
        raise ValueError('preregistration source bytes differ')
    result = observe(source)
    result.update(plan=plan, sources=[saved.binding(path) for path in
                  (args.source, args.plan, Path(__file__), Path(saved.__file__))])
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'], samples=len(result['samples']), checks=len(result['checks']))))


if __name__ == '__main__':
    main()
