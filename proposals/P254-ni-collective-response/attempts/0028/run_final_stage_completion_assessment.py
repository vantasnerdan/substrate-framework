#!/usr/bin/env python3
"""Empirical completion scenarios from the actual nonzero mixed-jet work.

Uses original Native costs, exact rational schedule arithmetic and the unchanged
producer deadline. Scenarios are not runtime bounds or permission to extend it.
"""
import argparse
from datetime import datetime
import json
from pathlib import Path

import sympy as sp

import run_stage_sensitive_completion_assessment as schedule


def checked_groups(entries, actions, expected_count):
    records = [entry['native'] for entry in entries]
    if [record['index'] for record in records] != list(range(expected_count)):
        raise ValueError('original ordered source coverage differs')
    groups = {cpu: [] for cpu in range(6, 11)}
    for record in records:
        if (record['cpu'] != 6 + record['index'] % 5
                or record['hamiltonian_actions'] != actions
                or record['blas_threads'] != 1):
            raise ValueError('original pinned source/action/thread assignment differs')
        duration = schedule.rational(record['elapsed_seconds'])
        if duration <= 0:
            raise ValueError('positive observed work-unit duration required')
        groups[record['cpu']].append(duration)
    return groups


def run(source, plan):
    moments = checked_groups(source['current_moments'], 0, len(source['current_moments']))
    transitions = checked_groups(source['completed_stage2_transitions'], 7, 128)
    cpus = list(range(6, 11))
    totals = schedule.assigned(128, len(cpus))
    # Deliberately favorable: an in-flight item is treated as entirely paid.
    remaining = [total - len(moments[cpu]) - 1 for total, cpu in zip(totals, cpus)]
    if min(remaining) < 0:
        raise ValueError('remaining source assignment differs')
    anchor = schedule.rational(datetime.fromisoformat(source['snapshot_utc']).timestamp())
    deadline = schedule.rational(source['original_deadline_unix'])
    scenarios = []
    for name, statistic, action_ratio in (
            ('observed_cpu_minimum_ignore_extra_eighth_action', min, sp.S.One),
            ('observed_cpu_median_ignore_extra_eighth_action', schedule.median, sp.S.One),
            ('observed_cpu_median_action_scaled8over7', schedule.median, sp.Rational(8, 7)),
            ('observed_cpu_maximum_action_scaled8over7', max, sp.Rational(8, 7))):
        moment_costs = [statistic(moments[cpu]) for cpu in cpus]
        transition_costs = [statistic(transitions[cpu]) for cpu in cpus]
        next_context = schedule.makespan(remaining, moment_costs)
        final_transition = schedule.makespan(totals, transition_costs) * action_ratio
        physical = anchor + next_context + final_transition
        scenarios.append(dict(name=name,
            per_cpu_moment_seconds=list(map(str, moment_costs)),
            per_cpu_seven_action_transition_seconds=list(map(str, transition_costs)),
            favorable_remaining_context_seconds=str(next_context),
            future_final_transition_seconds=str(final_transition),
            final_transition_action_time_scale=str(action_ratio),
            next_context_scenario_utc=schedule.timestamp(anchor + next_context),
            physical_endpoint_scenario_utc=schedule.timestamp(physical),
            physical_before_original_deadline=bool(physical < deadline),
            physical_seconds_past_original_deadline=str(physical - deadline),
            omitting_required_final_transition_before_deadline=bool(anchor + next_context < deadline)))
    # Single-purpose verified blocks remain distinct from the runtime assumption.
    if schedule.assigned(7, 3) != [3, 2, 2]:
        raise AssertionError('cyclic assignment example failed')
    if schedule.makespan([2, 1], [sp.Integer(3), sp.Integer(5)]) != 6:
        raise AssertionError('parallel makespan example failed')
    if schedule.makespan([2, 1], [sp.Integer(3), sp.Integer(5)]) == 11:
        raise AssertionError('serial instead of parallel mutation not sensitive')
    sensitive = any(row['physical_before_original_deadline'] !=
                    row['omitting_required_final_transition_before_deadline'] for row in scenarios)
    if not sensitive:
        raise AssertionError('omitted final RK transition mutation not sensitive')
    return dict(schema='P254-actual-final-stage-completion-scenarios-v1',
        status='actual-cost-scenarios-exercised-not-runtime-guarantee',
        anchor_utc=source['snapshot_utc'], original_deadline_utc=schedule.timestamp(deadline),
        observed_final_stage_moment_sources=str(len(source['current_moments']))+'/128',
        observed_previous_seven_action_transition_sources='128/128',
        total_cyclic_assignment=totals,
        favorable_remaining_moment_assignment_ignoring_one_inflight_per_worker=remaining,
        scenarios=scenarios,
        omitted_time='ALL global aggregation/startup/retirement overhead and ALL post-physical-endpoint moment preparation/bare-energy work. One in-flight moment per CPU is treated completely paid.',
        assumptions='Per-CPU observed minimum/median/maximum extrapolation, not a bound or confidence interval. Future inputs/caches may differ. Eight/seven proportional action cost is an explicit separate assumption; unscaled rows deliberately omit the extra callback.',
        checks=[dict(name='exact_cyclic_assignment_example', passed=True),
                dict(name='exact_parallel_makespan_example', passed=True),
                dict(name='serial_schedule_mutation_sensitive', passed=True, mutation=True),
                dict(name='omitted_required_final_RK_transition_changes_deadline_verdict',
                     passed=sensitive, mutation=True)],
        scope='Planning evidence only; original producer maximum, physical criteria, code and cumulative costs are not changed.',
        plan=plan)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source.read_bytes())
    plan = json.loads(args.plan.read_bytes())
    if plan['source'] != schedule.binding(args.source):
        raise ValueError('plan does not bind the exact original cost snapshot')
    result = run(source, plan)
    result['sources'] = [schedule.binding(path) for path in
        (args.source, args.plan, Path(__file__), Path(schedule.__file__))]
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
