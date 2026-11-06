"""Normalize the actual matched clamped controls beside completed active curves.

Consumes the frozen complete A29 publication only. No new trajectory, electronic
projection, noise fit, decay fit or physical-loss calibration is performed.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0029'), str(HERE.parent/'0028')]

import numpy as np
from clamped_scientific_readout import encode
from reduced_coherence_readout import normalized_ordered
from reduced_comparison_inputs import decode
from substrate_framework.streamed_material_state import file_binding


def unique_index(rows, key):
    out = {}
    for row in rows:
        identity = key(row)
        if identity in out:
            raise ValueError('duplicate physical observation identity: '+str(identity))
        out[identity] = row
    return out


def contrast_curve(curve, control):
    observations = unique_index(control['observations'], lambda row: row['time_inverse_ry'])
    samples = []
    for active_sample in curve['coherence']['samples']:
        time = active_sample['time_inverse_ry']
        sample = dict(time_inverse_ry=time)
        if not active_sample['connected_available'] or time not in observations:
            samples.append(dict(sample, supported=False, reason='missing actual matched-time covariance'))
            continue
        correlation = observations[time]['connected_correlation']
        reference = active_sample['reference_time_inverse_ry']
        if (not correlation['available'] or correlation['reference_time_inverse_ry'] != reference
                or reference not in observations
                or not observations[reference]['connected_correlation']['available']):
            samples.append(dict(sample, supported=False, reason='missing exact common covariance reference'))
            continue
        clamped = normalized_ordered(decode(correlation['two_time']['ordered_ab']),
            decode(correlation['equal_time']['ordered_ab'][0]),
            decode(observations[reference]['connected_correlation']['equal_time']['ordered_ab'][0]))
        active = decode(active_sample['ordered_coherence'])
        if not active['defined'] or not clamped['defined']:
            samples.append(dict(sample, supported=False, reason='an actual normalization is undefined',
                                active=active, clamped=clamped))
            continue
        active_value = np.asarray(active['complex_degree'], dtype=complex)
        clamped_value = np.asarray(clamped['complex_degree'], dtype=complex)
        if active_value.shape != clamped_value.shape:
            raise ValueError('matched controls must retain the same original spatial cells')
        samples.append(dict(sample, supported=True, reference_time_inverse_ry=reference,
            active=active, clamped=clamped,
            normalized_complex_difference=active_value-clamped_value,
            magnitude_difference=np.abs(active_value)-np.abs(clamped_value),
            active_noise_scope=active_sample['noise_scope'],
            clamped_noise_scope=correlation['noise_scope']))
    return dict(case=curve['case'], bands=curve['bands'], features=curve['features'],
        active_source=curve['source'], clamped_source=control['source'],
        clamped_identity=control['identity'], samples=samples)


def consume(comparison, expected_sha256):
    payload = Path(comparison).read_bytes()
    actual_sha = sha256(payload).hexdigest()
    if actual_sha != expected_sha256:
        raise ValueError('complete original comparison bytes differ from requested source')
    report = json.loads(payload)
    source = dict(path=str(Path(comparison).resolve()), bytes=len(payload), sha256=actual_sha)
    del payload
    curves = unique_index(report['actual_curves'], lambda row: row['source']['sha256'])
    controls = unique_index(report['canonical_clamped_controls'], lambda row: row['identity'])
    comparisons = report['matched_clamped_comparisons']
    eligible = [row for row in comparisons if row.get('supported_matched_preparation') is True]
    rows = [contrast_curve(curves[row['finite_result']['sha256']],
                           controls[row['clamped_control_identity']]) for row in eligible]
    return dict(schema='P254-saved-coupled-clamped-coherence-v1', status='complete',
        source=source, rows=rows, eligible_comparisons=len(eligible),
        explicit_unsupported_comparisons=sum(row.get('supported') is False for row in comparisons),
        implementation=[file_binding(__file__), file_binding(
            HERE.parent/'0029/reduced_coherence_readout.py')],
        units='dimensionless ordered coherence and raw coupled-minus-clamped differences',
        scope='Exact matched preparation and saved times; declared active-Gaussian and electrons-clamped noise models. Descriptive coupling/noise-model contrast, not omitted-physics error, an independent passive bath, irreversible damping, heat or a physical lifetime.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = consume(args.comparison, args.expected_sha256)
    with Path(args.output).open('x') as stream:
        json.dump(encode(result), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=file_binding(args.output), eligible=result['eligible_comparisons'],
                          unsupported=result['explicit_unsupported_comparisons'])))


if __name__ == '__main__':
    main()
