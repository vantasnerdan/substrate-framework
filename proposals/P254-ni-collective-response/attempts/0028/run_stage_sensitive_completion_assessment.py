#!/usr/bin/env python3
"""Favorable empirical scenarios from actual Native stage receipts, not an ETA."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import sympy as sp


def rational(value):
    return sp.Rational(str(value))


def median(values):
    values = sorted(values)
    mid = len(values)//2
    return values[mid] if len(values) % 2 else (values[mid-1]+values[mid])/2


def assigned(total, workers):
    return [len(range(index, total, workers)) for index in range(workers)]


def makespan(counts, costs):
    return max(count*cost for count, cost in zip(counts, costs))


def timestamp(unix):
    return datetime.fromtimestamp(float(unix), timezone.utc).isoformat(timespec='seconds')


def binding(path):
    raw = path.read_bytes()
    return dict(path=str(path), bytes=str(len(raw)), sha256=hashlib.sha256(raw).hexdigest())


def run(source, plan):
    records = [item['native'] for item in source['current_completed_prefix']]
    if [item['index'] for item in records] != list(range(20)):
        raise ValueError('worksheet requires the admitted20-source prefix')
    workers = list(range(6, 11))
    groups = {cpu: [] for cpu in workers}
    for item in records:
        if item['cpu'] != 6+item['index'] % 5 or item['hamiltonian_actions'] != 7:
            raise ValueError('actual pinned stage/assignment differs')
        groups[item['cpu']].append(rational(item['elapsed_seconds']))
    total_counts = assigned(128, len(workers))
    # Deliberately favorable: count the current in-flight item on each CPU paid.
    remaining_counts = [total-len(groups[cpu])-1 for total, cpu in zip(total_counts, workers)]
    if min(remaining_counts) < 0:
        raise ValueError('invalid admitted stage prefix')
    moment = (rational(source['moment_context_original_budget']['elapsed_seconds'])
              -rational(source['first_global_stage_original_budget']['elapsed_seconds']))
    if moment <= 0:
        raise ValueError('positive measured complete moment span required')
    anchor = rational(datetime.fromisoformat(plan['anchor_utc']).timestamp())
    deadline = rational(plan['original_deadline_unix'])
    scenarios = []
    for name, statistic, ratio in (
            ('observed_cpu_minimum_no_extra_eighth_action', min, sp.S.One),
            ('observed_cpu_median_no_extra_eighth_action', median, sp.S.One),
            ('observed_cpu_median_action_scaled8over7', median, sp.Rational(8,7)),
            ('observed_cpu_maximum_action_scaled8over7', max, sp.Rational(8,7))):
        costs = [statistic(groups[cpu]) for cpu in workers]
        current = makespan(remaining_counts, costs)
        final = makespan(total_counts, costs)*ratio
        physical = anchor+current+moment+final
        endpoint_prepared = physical+moment
        scenarios.append(dict(name=name,
            per_cpu_source_seconds=[str(value) for value in costs],
            current_stage_remaining_seconds=str(current),
            next_final_transition_seconds=str(final),
            one_future_moment_seconds=str(moment),
            physical_endpoint_scenario_utc=timestamp(physical),
            final_endpoint_moments_before_bare_H0_scenario_utc=timestamp(endpoint_prepared),
            physical_endpoint_before_deadline=bool(physical < deadline),
            endpoint_moments_before_deadline=bool(endpoint_prepared < deadline),
            seconds_past_deadline_before_ANY_bare_energy=str(endpoint_prepared-deadline),
            final_transition_action_time_scale=str(ratio)))
    toy = makespan([1,0,0], [sp.S.One, sp.Integer(2), sp.Integer(3)])
    if toy != 1 or assigned(7,3) != [3,2,2]:
        raise AssertionError('exact cyclic schedule block failed')
    if makespan([4,4], [sp.S.One, sp.Integer(3)]) == 4+12:
        raise AssertionError('serial-schedule mutation was not distinguished')
    changed_verdict = any(item['physical_endpoint_before_deadline']
                          != item['endpoint_moments_before_deadline'] for item in scenarios)
    if not changed_verdict:
        raise AssertionError('omitted-required-final-context mutation not sensitive on this data')
    return dict(schema='P254-stage-sensitive-completion-scenarios-v1',
        status='scenario-model-exercised-not-runtime-guarantee',
        anchor_utc=plan['anchor_utc'], original_deadline_utc=timestamp(deadline),
        observed_sources='20/128', actual_stage='2', actual_counted_actions_per_source='7',
        future_structural_final_actions_per_source='8', workers=workers,
        total_assignment=total_counts,
        favorable_remaining_assignment_ignoring_one_inflight_per_worker=remaining_counts,
        observed_per_cpu_elapsed_seconds={str(cpu): [str(value) for value in values]
                                         for cpu, values in groups.items()},
        measured_full_nonzero_moment_span_seconds=str(moment),
        future_moment_barriers='2; final-evaluation input AND distinct physical endpoint',
        omitted_time='ALL future bare-H0 energy time, all future global stage aggregation and startup; no extra cost for nonzero mixed-jet moment terms. One in-flight source per worker treated fully paid.',
        assumptions='Observed per-CPU timing range extrapolated, not a hard bound or confidence interval. Future data/work/cache may differ. Action-scaled8/7 is a separate proportional-cost assumption; unscaled rows deliberately ignore the extra callback.',
        scenarios=scenarios,
        checks=[dict(name='exact_cyclic_schedule_example', passed=True),
                dict(name='serial_instead_of_parallel_sensitive', passed=True, mutation=True),
                dict(name='omitting_required_physical_endpoint_prepare_changes_deadline_verdict',
                     passed=changed_verdict, mutation=True)],
        conclusion='Physical endpoint may fit some favorable scenarios. Full11-energy completion is not a credible promise under this owner when even final-state moment preparation is late in these scenarios before any bare-energy time. Preserve original deadline; evaluate standalone no-replay post-endpoint readout, not internal-state adoption or clock reset.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, plan = json.loads(args.source.read_bytes()), json.loads(args.plan.read_bytes())
    result = run(source, plan)
    result['sources'] = [binding(path) for path in (args.source, args.plan, Path(__file__))]
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
