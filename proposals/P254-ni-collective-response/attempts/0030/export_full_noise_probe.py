"""Export ONE original t=0 sector's unweighted physical force action pair.

--output is a new directory outside the repository. packet.json is published
last; existing outputs (including interrupted work) are never overwritten.
No retained/omitted covariance or trajectory is evaluated by this CLI.
"""
import argparse
import json
import os
from pathlib import Path
import resource
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial-result', required=True)
    parser.add_argument('--coefficient-manifest', required=True)
    parser.add_argument('--original-record', type=int, required=True)
    parser.add_argument('--atom', type=int, required=True)
    parser.add_argument('--axis', choices=('x', 'y', 'z'), required=True)
    parser.add_argument('--bands', type=int, choices=(12,), required=True)
    parser.add_argument('--output', required=True, help='exclusive new packet directory outside Git')
    return parser.parse_args()


def runtime():
    if sorted(os.sched_getaffinity(0)) != [7]:
        raise ValueError('pilot requires exactly CPU7 affinity')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise ValueError(name + '=1 required before numerical imports')


def export_packet(inputs, directory):
    """Write one action pair and small coordinates; keep original Psi readonly."""
    import numpy as np
    from substrate_framework.streamed_material_state import file_binding, save_array
    started = perf_counter()
    with (directory / 'coordinates.npz').open('xb') as stream:
        np.savez(stream, gram=inputs.gram, density_coefficients=inputs.density_coefficients,
                 active_coordinates=inputs.active_coordinates, selected=inputs.selected)
    coordinates = file_binding(directory / 'coordinates.npz')
    probes = {}
    for index, (name, probe) in enumerate(inputs.probes.items()):
        probes[name] = dict(
            action=save_array(directory / ('probe-%02d-action.npy' % index), probe.action),
            adjoint_action=save_array(directory / ('probe-%02d-adjoint-action.npy' % index), probe.adjoint_action),
            operator_identity=probe.operator_identity, adjoint_identity=probe.adjoint_identity,
            units=probe.units, pw_basis_identity=probe.pw_basis_identity,
            source_authority=probe.source_authority, shape=list(probe.action.shape))
    work = dict(inputs.operator_work)
    work.update(packet_array_write_seconds=perf_counter() - started,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)
    packet = dict(schema='P254-full-noise-probe-input-v1', status='complete',
        frame=inputs.frame, coordinates=coordinates, probes=probes,
        factor_identity=inputs.factor_identity, pw_basis_identity=inputs.pw_basis_identity,
        state_authority=inputs.state_authority, source=inputs.source,
        quadrature_weight=inputs.quadrature_weight, quadrature_applied=False,
        time_inverse_ry=0, external_q=[0, 0, 0], operator_work=work,
        scope='one original source, one physical electronic force; no wholeE or late-time claim')
    with (directory / 'packet.json').open('x') as stream:
        json.dump(packet, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    return file_binding(directory / 'packet.json')


def main():
    args = arguments()
    runtime()
    directory = Path(args.output).resolve()
    if directory == ROOT or ROOT in directory.parents:
        raise ValueError('full-PW numerical packet must remain outside Git')
    # Claim exclusively before any physical action; never replace earlier work.
    directory.mkdir(parents=True, exist_ok=False)
    from full_noise_probe_inputs import prepare_full_noise_probe
    inputs = prepare_full_noise_probe(args.initial_result, args.coefficient_manifest,
        args.original_record, args.atom, args.axis, args.bands)
    binding = export_packet(inputs, directory)
    print(json.dumps(dict(packet=binding, operator_work=inputs.operator_work), allow_nan=False))


if __name__ == '__main__':
    main()
