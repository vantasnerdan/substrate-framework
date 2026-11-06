#!/usr/bin/env python3
"""Copy compact scientific observations from completed, byte-bound publications.

No propagation, numerical readout, comparison, fit or uncertainty is recomputed.
The original command, publication and failed orchestration receipts are retained.
"""
import argparse
import json
from pathlib import Path

from material_run_contract import binding, bound_json
from matched_reference_io import save_json
from run_matched_clamped_batch import load_result


def scalar_fields(value):
    """Copy dictionaries and scalar leaves; omit potentially large arrays."""
    if isinstance(value, dict):
        return {key: scalar_fields(item) for key, item in value.items()
                if not isinstance(item, list)}
    return value


def compact_readout(readout, readout_index):
    result = {key: readout[key] for key in (
        'case', 'control', 'mass', 'original_result', 'original_input',
        'field_authority', 'kick', 'limits', 'whole_coupled_lattice_ranking')}
    result['observations'] = []
    for index, observation in enumerate(readout['observations']):
        row = dict(
            source_json_pointer=f'/readouts/{readout_index}/observations/{index}',
            time_inverse_ry=observation['time_inverse_ry'],
            reference_time_inverse_ry=observation['reference_time_inverse_ry'],
            producer_diagnostics=observation['producer_diagnostics'], states={})
        for state_name, state in observation['states'].items():
            compact = dict(maxwell=scalar_fields(state['maxwell']))
            for motion in ('absolute_displacement', 'kick_delta_displacement'):
                if motion in state:
                    compact[motion] = {
                        species: scalar_fields(state[motion][species])
                        for species in ('all_ions', 'Ni', 'H')}
            row['states'][state_name] = compact
        result['observations'].append(row)
    return result


def consume_publication(job, command_directory):
    intent_binding = binding(command_directory/'intent.json')
    intent = bound_json(intent_binding)
    if intent['argv'] != job['argv']:
        raise ValueError('original command differs from planned producer')
    outcome_binding = binding(command_directory/'outcome.json')
    outcome = bound_json(outcome_binding)
    if outcome['returncode'] != 0 or outcome['timed_out']:
        raise ValueError('original scientific producer did not succeed')
    publication = bound_json(outcome['stdout'])
    if publication['status'] != 'complete':
        raise ValueError('producer did not publish a complete result')
    original = publication['output']
    if Path(original['path']).resolve() != Path(job['receipt']).resolve():
        raise ValueError('publication is not the planned scientific result')
    value = load_result(original)
    if value['status'] != 'complete':
        raise ValueError('saved scientific result is incomplete')
    return dict(
        schema='P254-clamped-scientific-compact-summary-v1', status='complete',
        original_scientific_output=original, original_command=intent_binding,
        original_command_outcome=outcome_binding, original_stdout=outcome['stdout'],
        projection='Exact copied scalar/dictionary fields; arrays and per-atom '
                   'tables omitted. No new calculation or uncertainty estimate.',
        oracle=value['oracle'], source_bindings=value['source_bindings'],
        numerical_accuracy_certified=value['numerical_accuracy_certified'],
        bound_pair_comparison_pointer='/bound_pair_comparison',
        readouts=[compact_readout(item, index)
                  for index, item in enumerate(value['readouts'])])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan_binding = binding(args.plan)
    plan = bound_json(plan_binding)
    prior_binding = binding(args.batch/'result.json')
    prior = bound_json(prior_binding)
    if prior['plan'] != plan_binding:
        raise ValueError('original orchestration does not bind this plan')
    args.output.mkdir(parents=True, exist_ok=False)
    completed = {}
    try:
        for index, job in enumerate(plan['jobs']):
            summary = consume_publication(job, args.batch/f'command-{index:03d}')
            path = args.output/(job['key']+'.json')
            save_json(path, summary)
            completed[job['key']] = binding(path)
            print(json.dumps(dict(event='summary-complete', key=job['key'],
                                  output=completed[job['key']])), flush=True)
    finally:
        complete = len(completed) == len(plan['jobs'])
        save_json(args.output/'result.json', dict(
            schema='P254-clamped-scientific-summary-index-v1',
            status='complete' if complete else 'incomplete',
            original_plan=plan_binding, original_orchestration=prior_binding,
            original_orchestration_status=prior['status'],
            implementation_sources=[binding(__file__),
                                    binding(Path(__file__).with_name(
                                        'run_matched_clamped_batch.py'))],
            completed=completed, recomputed_scientific_results=False,
            scope='Clamped-electron observations only; whole coupled E/D/C open'))
        print('[dft-done] A28 scientific summaries' if complete else
              '[dft-fail] A28 scientific summaries', flush=True)


if __name__ == '__main__':
    main()
