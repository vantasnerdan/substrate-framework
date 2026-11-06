#!/usr/bin/env python3
"""Consume original completed matched-clamped output; optionally bind a control pair."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3]/'src'), str(HERE)]
from clamped_scientific_readout import (
    REGISTRATION, binding, bound_json, compare_pair, encode, load_readout, sources,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', type=Path, required=True)
    parser.add_argument('--result-sha256', required=True)
    parser.add_argument('--independent-result', type=Path)
    parser.add_argument('--independent-result-sha256')
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if bool(args.independent_result) != bool(args.independent_result_sha256):
        parser.error('independent result path and SHA256 are required together')
    current = sources()
    registration = json.loads(REGISTRATION.read_text())
    oracle_binding = binding(args.oracle)
    oracle = bound_json(oracle_binding)
    if (oracle.get('schema') != 'P254-clamped-scientific-readout-oracle-v1'
            or oracle.get('status') != 'passed' or oracle.get('source_bindings') != current
            or sorted(oracle.get('checks', {})) != sorted(registration['checks'])
            or sorted(oracle.get('mutation_witnesses', {})) != sorted(registration['mutation_witnesses'])):
        raise ValueError('focused passed oracle for exact current source inventory required')
    readouts = [load_readout(args.result, args.result_sha256)]
    comparison = None
    if args.independent_result:
        readouts.append(load_readout(args.independent_result, args.independent_result_sha256))
        comparison = compare_pair(*readouts)
    output = dict(schema='P254-clamped-scientific-readout-bundle-v1', status='complete',
                  source_bindings=current, oracle=oracle_binding, readouts=readouts,
                  bound_pair_comparison=comparison, numerical_accuracy_certified=False)
    with args.output.open('x') as stream:
        json.dump(encode(output), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='complete', output=binding(args.output))))


if __name__ == '__main__':
    main()
