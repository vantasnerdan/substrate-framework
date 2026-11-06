#!/usr/bin/env python3
"""Carry the settled field oracle; exercise actual source-authority publication."""
import argparse
import json
from pathlib import Path

from run_full_field_moments_material import (
    HERE, binding, bound_json, implementation_bindings, load_inputs, pack, preload_parsers)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prior_binding = binding(HERE/'full-field-moments-oracle-001.json')
    prior = bound_json(prior_binding)
    assert prior['status'] == 'passed'
    current = implementation_bindings()
    current_by_path = {entry['path']: entry for entry in current}
    preserved = []
    for entry in prior['implementation_sources']:
        if Path(entry['path']).name == 'run_full_field_moments_material.py':
            original = Path(entry['path'])
            saved = binding(HERE/'source-generations'/(original.stem+'-'+entry['sha256']+original.suffix))
            assert saved['sha256'] == entry['sha256'] and saved['bytes'] == entry['bytes']
            preserved.append(saved)
        else:
            assert current_by_path[entry['path']] == entry
    assert prior['registration'] == binding(HERE/'full-field-moments-registration.json')
    preload_parsers()
    _, _, _, _, _, grid = load_inputs(HERE/'full-common-reference-004/full-common-reference.json')
    actual = grid.slots[0].source_binding
    publication = args.output.with_name(args.output.stem+'-source-binding.json')
    with publication.open('x') as stream:
        json.dump(pack({'source_binding': actual}), stream, indent=2, allow_nan=False)
        stream.write('\n')
    restored = json.loads(publication.read_bytes())
    assert restored['source_binding'] == dict(actual)
    result = dict(prior, scope='Original independent numerical oracle carried forward; actual saved-source binding publication adapter exercised',
        implementation_sources=current, prior_scientific_oracle=prior_binding,
        original_consumer_sources=preserved,
        publication_adapter=dict(source_type=type(actual).__name__, source_binding=pack(actual),
            saved_output=binding(publication), actual_saved_output_matches=True,
            scientific_oracle_reexecuted=False, adapter_source=binding(Path(__file__))))
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed', result=binding(args.output))))


if __name__ == '__main__':
    main()
