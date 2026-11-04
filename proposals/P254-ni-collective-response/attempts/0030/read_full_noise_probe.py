"""Consume one authenticated physical-force packet without repeating its actions.

The ordered four-block decomposition uses the original raw state and admitted
Gram. This is one initial electronic sector, not a full-q or late-time result.
"""
import argparse
import json
from pathlib import Path
import resource
from time import perf_counter


def complex_record(value):
    value = complex(value)
    return dict(real=value.real, imag=value.imag)


def consume_packet(packet_path):
    import numpy as np
    from substrate_framework.material_control_reference import bound_arrays, bound_json
    from substrate_framework.material_noise_projection import projected_noise_decomposition
    from substrate_framework.physical_factor_preparation import ActionFactors
    from substrate_framework.streamed_material_state import authenticate, file_binding

    started = perf_counter()
    binding = file_binding(packet_path)
    packet = bound_json(binding)
    if (packet['schema'] != 'P254-full-noise-probe-input-v1'
            or packet['status'] != 'complete' or packet['time_inverse_ry'] != 0
            or packet['quadrature_applied'] is not False):
        raise ValueError('complete original unweighted t0 physical packet required')
    arrays = bound_arrays(packet['coordinates'])
    psi = np.load(authenticate(packet['frame']), mmap_mode='r', allow_pickle=False)
    probes = {}
    for name, probe in packet['probes'].items():
        action = np.load(authenticate(probe['action']), mmap_mode='r', allow_pickle=False)
        adjoint = np.load(authenticate(probe['adjoint_action']), mmap_mode='r', allow_pickle=False)
        probes[name] = ActionFactors(action, adjoint, probe['operator_identity'],
            probe['adjoint_identity'], probe['units'], probe['pw_basis_identity'],
            probe['source_authority'])
    restored = perf_counter()
    result = projected_noise_decomposition(psi, arrays['gram'],
        arrays['density_coefficients'], arrays['active_coordinates'], probes,
        factor_identity=packet['factor_identity'], pw_basis_identity=packet['pw_basis_identity'],
        state_authority=packet['state_authority'])
    rows = []
    for (a, b), value in result['full'].items():
        rows.append(dict(probe_a=a, probe_b=b, full=complex_record(value),
            blocks=[dict(part_a=pa, part_b=pb,
                         value=complex_record(result['blocks'][a, pa, b, pb]))
                    for pa in ('retained', 'omitted') for pb in ('retained', 'omitted')],
            reconstruction_defect=complex_record(result['reconstruction_defect'][a, b])))
    return dict(schema='P254-original-sector-noise-decomposition-v1', status='complete',
        input=binding, frame=packet['frame'], source=packet['source'],
        quadrature_weight=packet['quadrature_weight'], quadrature_applied=False,
        rows=rows, means=[dict(probe=name, part=part, value=complex_record(value))
                        for (name, part), value in result['means'].items()],
        metadata=result['metadata'], work=dict(restoration_seconds=restored-started,
            contraction_seconds=perf_counter()-restored, elapsed_seconds=perf_counter()-started,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024),
        implementation=[file_binding(__file__), file_binding(
            Path(__file__).resolve().parents[4]/'src/substrate_framework/material_noise_projection.py')],
        scope='One original unweighted electronic sector at t0; original full state, retained/omitted cross terms and empty-PW complement. Not full128/64q, a late-time kernel, damping or a lifetime.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    from export_full_noise_probe import runtime
    runtime()
    result = consume_packet(args.packet)
    with Path(args.output).open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    from substrate_framework.streamed_material_state import file_binding
    print(json.dumps(dict(output=file_binding(args.output), work=result['work']), allow_nan=False))


if __name__ == '__main__':
    main()
