"""Source0 two-mode original-PW inputs from paid physical legs, no replay."""
import json
from pathlib import Path

import numpy as np
from full_noise_probe_inputs import FullNoiseProbeInputs
from original_noise_source import _same_binding
from substrate_framework.physical_factor_preparation import ActionFactors
from substrate_framework.streamed_material_state import authenticate, file_binding


def _json(path):
    binding = file_binding(path)
    with authenticate(binding).open() as stream:
        return binding, json.load(stream)


def cached_source0(initial_result, coefficient_manifest, mode_archive, directions,
                   action_packet, second_packet):
    action_binding, packet = _json(action_packet)
    mode1_binding, second = _json(second_packet)
    initial_binding, initial = _json(initial_result)
    if (initial['schema'] != 'P254-matched-joint-case-v1' or initial['status'] != 'complete'
            or initial['mode'] != 'initialize' or initial['material_trajectory'] is not False):
        raise ValueError('completed original matched initialization required')
    source = dict(packet['source'])
    if (packet['status'] != 'complete' or second['status'] != 'complete'
            or source['original_record'] != 0
            or not _same_binding(second['original_action_packet'], action_binding)
            or not _same_binding(source['initial_result'], initial_binding)
            or not _same_binding(source['coefficient_manifest'], file_binding(coefficient_manifest))
            or not _same_binding(packet['channel_archive'], mode_archive)):
        raise ValueError('paid original source0/state/mode lineage differs')
    for name in ('preparation', 'original_boundary'):
        if not _same_binding(source[name], initial['saved'][name]):
            raise ValueError('paid source0 original authority differs: ' + name)
        authenticate(initial['saved'][name])
    saved_direction = np.asarray(second['direction_real'])+1j*np.asarray(second['direction_imag'])
    if not np.array_equal(saved_direction, directions[:, 1]):
        raise ValueError('actual paid second-mode direction differs')
    with np.load(authenticate(packet['coordinates']), allow_pickle=False) as archive:
        coordinates = {name:archive[name] for name in archive.files}
    psi = np.load(authenticate(packet['frame']), mmap_mode='r', allow_pickle=False)
    raw_d = np.load(authenticate(source['original_raw_d']), mmap_mode='r', allow_pickle=False)
    if not np.array_equal(raw_d, coordinates['density_coefficients']):
        raise ValueError('unchanged original raw D differs from paid coordinates')
    records = {'mode0':packet['actions']['ni_output_ionic_mode0'],
               'mode1':second['actions']}
    probes = {}
    for name, record in records.items():
        if (record['pw_basis_identity'] != packet['pw_basis_identity']
                or record['source_authority'] != packet['state_authority']):
            raise ValueError('paid force basis/state authority differs')
        probes[name] = ActionFactors(
            np.load(authenticate(record['action']), mmap_mode='r', allow_pickle=False),
            np.load(authenticate(record['adjoint_action']), mmap_mode='r', allow_pickle=False),
            record['operator_identity'], record['adjoint_identity'], record['units'],
            record['pw_basis_identity'], record['source_authority'])
        probes[name].validate(psi.shape, packet['pw_basis_identity'])
    source.pop('operator', None)  # The old H37 probe is lineage, not this observable.
    source['operators'] = records
    source['paid_mode_action_packets'] = [action_binding, mode1_binding]
    return FullNoiseProbeInputs(psi, coordinates['gram'], raw_d,
        coordinates['active_coordinates'], coordinates['selected'], probes,
        packet['factor_identity'], packet['pw_basis_identity'], packet['state_authority'],
        packet['frame'], source, packet['quadrature_weight'],
        dict(force_forward_actions=0, force_adjoint_actions=0,
            h0_actions=0, trajectory_actions=0, native_calls=0,
            gram_recomputations=0, population_refills=0, preparation_replays=0,
            reused_paid_force_legs=4, raw_action_payloads_written=False))
