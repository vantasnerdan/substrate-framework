"""Actual original-sector full-PW two-mode initial force-noise producer.

No H action or trajectory. Source0 uses paid legs; every other sector evaluates
only the two requested force/adjoint pairs. Complete scalar covariance outputs
are checkpointed; in-case force arrays are temporary, not a recovery promise.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src')]+[str(HERE.parent/name) for name in ('0030','0029','0028')]
import numpy as np
from cached_mode_source0 import cached_source0
from full_noise_probe_inputs import prepare_full_noise_modes
from export_output_directed_actions import covariance_rows
from clamped_scientific_readout import encode
from substrate_framework.material_noise_projection import projected_noise_decomposition
from substrate_framework.streamed_material_state import authenticate, file_binding
from checkpoint_directories import make_directory, sync_directory
from mode_noise_compatibility import load_compatibility


def save_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(encode(value), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    sync_directory(Path(path).parent)
    return file_binding(path)


def run(args):
    started = perf_counter()
    compatibility = load_compatibility()
    output = Path(args.output).resolve()
    make_directory(output, exist_ok=False)
    mode_binding = file_binding(args.channel_archive)
    with np.load(authenticate(mode_binding), allow_pickle=False) as archive:
        directions = archive['ionic_basis'][:, :2]
    if args.original_record == 0:
        inputs = cached_source0(args.initial_result, args.coefficient_manifest,
            mode_binding, directions, args.reuse_source0_action_packet,
            args.reuse_source0_mode1_packet)
    else:
        print('[physics] restoring original record '+str(args.original_record)+
              '; only two force/adjoint pairs', flush=True)
        inputs = prepare_full_noise_modes(args.initial_result,
            args.coefficient_manifest, args.original_record, directions)
    if set(inputs.probes) != {'mode0','mode1'}:
        raise ValueError('exact requested two-mode physical action set required')
    contraction_started = perf_counter()
    moments = projected_noise_decomposition(**inputs.algebra_arguments())
    contraction_seconds = perf_counter()-contraction_started
    result = dict(schema='P254-A34-original-sector-two-mode-covariance-v1',
        status='complete', original_record=args.original_record,
        scientific_compatibility_identity=compatibility['scientific_compatibility_identity'],
        quadrature_weight=inputs.quadrature_weight, quadrature_applied=False,
        source=inputs.source, frame=inputs.frame,
        factor_identity=inputs.factor_identity, pw_basis_identity=inputs.pw_basis_identity,
        state_authority=inputs.state_authority, mode_archive=mode_binding,
        mode_indices=[0,1], covariance=covariance_rows(moments),
        covariance_metadata=moments['metadata'],
        work=dict(elapsed_seconds=perf_counter()-started,
            contraction_seconds=contraction_seconds,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024,
            physical_force_forward_actions=inputs.operator_work['force_forward_actions'],
            physical_force_adjoint_actions=inputs.operator_work['force_adjoint_actions'],
            hamiltonian_actions=0, trajectory_actions=0,
            physical_input_work=inputs.operator_work, raw_action_payloads_written=False),
        implementation=[file_binding(path) for path in (Path(__file__),
            HERE/'cached_mode_source0.py', HERE.parent/'0030/original_noise_source.py',
            HERE.parent/'0030/full_noise_probe_inputs.py')],
        scope='One original full-PW/raw-D internal source at external Gamma/t0; '
            'unweighted two physical force modes with retained/omitted/cross orders. '
            'Not late Ni displacement, all64 externalq, a dynamic noise bound, '
            'irreversible loss, lifetime, heat or full campaign closure.')
    result_binding = save_json(output/'result.json', result)
    save_json(output/'done.json', dict(result=result_binding))
    print(json.dumps(dict(result=result_binding, original_record=args.original_record,
        work=result['work']), sort_keys=True), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial-result', required=True)
    parser.add_argument('--coefficient-manifest', required=True)
    parser.add_argument('--channel-archive', required=True)
    parser.add_argument('--original-record', type=int, required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reuse-source0-action-packet', required=True)
    parser.add_argument('--reuse-source0-mode1-packet', required=True)
    run(parser.parse_args())
