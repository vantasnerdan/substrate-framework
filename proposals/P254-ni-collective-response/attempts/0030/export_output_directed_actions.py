"""One original-sector literal H and output-directed ionic force pilot.

Root484 permits one <=512MiB reusable action packet, not an all-source bank.
No stationary state, energy-label H, passive-time kernel or loss is inferred.
"""
import argparse
import json
from pathlib import Path
import resource
from time import perf_counter


def restore_filter(path):
    import numpy as np
    from canonical_force_filter import CanonicalForceFilter
    from substrate_framework.streamed_material_state import authenticate, file_binding
    binding = file_binding(path)
    if binding['sha256'] != 'c19a1b0a19566c510e5557761fbe6b67b2540e05d62d4ee317cf38f6593fef80':
        raise ValueError('the actual complete canonical filter packet is required')
    packet = json.loads(Path(path).read_bytes())
    with np.load(authenticate(packet['archive']), allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    return CanonicalForceFilter(packet['filter']['final_time_inverse_ry'],
        metadata=packet['filter'], **arrays), binding, packet['archive']


def save_channels(result, path):
    import numpy as np
    from substrate_framework.streamed_material_state import file_binding
    arrays = {key: result[key] for key in ('segment_bounds', 'boundary_times')}
    for kind in ('ionic', 'field'):
        for key in ('basis', 'gram', 'eigenvalues', 'modal_chebyshev', 'modal_boundary'):
            arrays[kind+'_'+key] = result[kind][key]
        for key, value in result[kind]['diagnostics'].items():
            if key != 'max_abs':
                arrays[kind+'_diagnostic_'+key] = value
    with path.open('xb') as stream:
        np.savez(stream, **arrays)
    return file_binding(path)


def covariance_rows(result):
    from read_full_noise_probe import complex_record
    return [dict(left=a, right=b, full=complex_record(value),
        blocks=[dict(part_left=pa, part_right=pb,
                     value=complex_record(result['blocks'][a, pa, b, pb]))
                for pa in ('retained', 'omitted') for pb in ('retained', 'omitted')],
        reconstruction_defect=complex_record(result['reconstruction_defect'][a, b]))
        for (a, b), value in result['full'].items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--filter-packet', required=True)
    parser.add_argument('--original-packet', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    from export_full_noise_probe import runtime
    runtime()
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    from clamped_scientific_readout import encode
    from force_channel_representation import decompose_force_channels
    from packet_hamiltonian_actions import restore_packet_hamiltonian
    from substrate_framework.material_noise_projection import projected_noise_decomposition
    from substrate_framework.streamed_material_state import file_binding, save_array

    started = perf_counter()
    filter_, filter_binding, filter_archive = restore_filter(args.filter_packet)
    channels = decompose_force_channels(filter_)
    channel_archive = save_channels(channels, directory/'force-channels.npz')
    context = restore_packet_hamiltonian(args.original_packet)
    h, h_work = context.hamiltonian_actions()
    direction = channels['ionic']['basis'][:, 0]
    force, force_work = context.ionic_force_actions(direction)
    probes = {'literal_initial_hamiltonian': h, 'ni_output_ionic_mode0': force}
    payload_bytes = sum(probe.action.nbytes+probe.adjoint_action.nbytes for probe in probes.values())
    if payload_bytes > 512*1024**2:
        raise ValueError('single-source reusable packet exceeds Root484 allocation')
    action_bindings = {}
    for name, probe in probes.items():
        action_bindings[name] = dict(
            action=save_array(directory/(name+'-action.npy'), probe.action),
            adjoint_action=save_array(directory/(name+'-adjoint.npy'), probe.adjoint_action),
            operator_identity=probe.operator_identity, adjoint_identity=probe.adjoint_identity,
            units=probe.units, pw_basis_identity=probe.pw_basis_identity,
            source_authority=probe.source_authority)
    moments = projected_noise_decomposition(context.psi, context.gram,
        context.density_coefficients, context.active_coordinates, probes,
        factor_identity=context.factor_identity, pw_basis_identity=context.pw_basis_identity,
        state_authority=context.state_authority)
    packet = dict(schema='P254-output-directed-original-sector-actions-v1', status='complete',
        source=context.source, frame=context.frame, coordinates=context.source['coordinates'],
        factor_identity=context.factor_identity, pw_basis_identity=context.pw_basis_identity,
        state_authority=context.state_authority, quadrature_weight=context.quadrature_weight,
        quadrature_applied=False, filter_packet=filter_binding, filter_archive=filter_archive,
        channel_archive=channel_archive, channel_metadata=channels['metadata'],
        channel_spectra={kind: channels[kind]['eigenvalues'] for kind in ('ionic', 'field')},
        channel_diagnostics={kind: channels[kind]['diagnostics']['max_abs'] for kind in ('ionic', 'field')},
        selected_physical_mode=dict(kind='ionic', mode=0, direction=direction,
            convention='sum_i basis[i,0] F_i;no conjugation;all channel modes remain in archive'),
        actions=action_bindings, ordered_covariance=covariance_rows(moments),
        covariance_metadata=moments['metadata'],
        work=dict(restoration=context.operator_work, hamiltonian=h_work, force=force_work,
            elapsed_seconds=perf_counter()-started, action_payload_bytes=payload_bytes,
            peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024),
        implementation=[file_binding(__file__), file_binding(Path(__file__).with_name(
            'force_channel_representation.py')), file_binding(Path(__file__).with_name(
            'packet_hamiltonian_actions.py'))],
        scope='One unweighted original sector,t0,literal saved matched H and one output-directed ionic-force mode. Reusable working factors,not a late passive time kernel,response correction,physical loss,lattice ranking or complete campaign.')
    packet_path = directory/'packet.json'
    with packet_path.open('x') as stream:
        json.dump(encode(packet), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(packet=file_binding(packet_path), elapsed_seconds=packet['work']['elapsed_seconds'],
                          action_payload_bytes=payload_bytes), allow_nan=False))


if __name__ == '__main__':
    main()
