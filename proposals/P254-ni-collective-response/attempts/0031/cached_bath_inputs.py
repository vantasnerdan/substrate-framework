"""Reuse paid source0 physical actions; add only the second selected force mode."""
from pathlib import Path
import json
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0030'),
                str(HERE.parent/'0029'), str(HERE.parent/'0028')]

import numpy as np
from export_output_directed_actions import restore_filter
from packet_hamiltonian_actions import restore_packet_hamiltonian
from clamped_scientific_readout import encode
from substrate_framework.physical_factor_preparation import ActionFactors
from substrate_framework.streamed_material_state import authenticate, file_binding, save_array


def read_action(value):
    return ActionFactors(
        np.load(authenticate(value['action']), mmap_mode='r', allow_pickle=False),
        np.load(authenticate(value['adjoint_action']), mmap_mode='r', allow_pickle=False),
        value['operator_identity'], value['adjoint_identity'], value['units'],
        value['pw_basis_identity'], value['source_authority'])


def save_action(value, directory, prefix):
    return dict(action=save_array(directory/(prefix+'-action.npy'), value.action),
        adjoint_action=save_array(directory/(prefix+'-adjoint.npy'), value.adjoint_action),
        operator_identity=value.operator_identity, adjoint_identity=value.adjoint_identity,
        units=value.units, pw_basis_identity=value.pw_basis_identity,
        source_authority=value.source_authority)


class CachedBathInputs:
    """One restored context, cached raw actions, and counted NEW-column actions."""
    def __init__(self, action_packet, filter_packet, directory, extra_rank):
        started = perf_counter()
        self.directory = Path(directory)
        self.action_binding = file_binding(action_packet)
        if self.action_binding['sha256'] != 'b6dc1bc91025f198e28559080f489cd5a7b4aaecb22908a7113f6361a5ed0521':
            raise ValueError('the actual settled output-directed source0 packet is required')
        packet = json.loads(Path(action_packet).read_bytes())
        self.filter, self.filter_binding, self.filter_archive = restore_filter(filter_packet)
        with np.load(authenticate(packet['channel_archive']), allow_pickle=False) as archive:
            self.directions = archive['ionic_basis'][:, :2]
            self.modal_coefficients = archive['ionic_modal_chebyshev'][:, :, :2, :]
            self.segment_bounds = archive['segment_bounds']
        self.context = restore_packet_hamiltonian(packet['source']['packet']['path'])
        if (self.context.state_authority != packet['state_authority']
                or self.context.pw_basis_identity != packet['pw_basis_identity']):
            raise ValueError('cached actions and restored original state/basis differ')
        self.hamiltonian = read_action(packet['actions']['literal_initial_hamiltonian'])
        self.forces = {'mode0': read_action(packet['actions']['ni_output_ionic_mode0'])}
        mode1, self.mode1_work = self.context.ionic_force_actions(self.directions[:, 1])
        self.forces['mode1'] = mode1
        mode1_binding = save_action(mode1, self.directory, 'mode1')
        self.mode1_packet = dict(status='complete', original_action_packet=self.action_binding,
            source=self.context.source, direction_real=self.directions[:, 1].real.tolist(),
            direction_imag=self.directions[:, 1].imag.tolist(), actions=mode1_binding,
            work=self.mode1_work, scope='one necessary new physical force mode;cached H/mode0 not replayed')
        with (self.directory/'mode1-complete.json').open('x') as stream:
            json.dump(encode(self.mode1_packet), stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        self.prepared = {name: self.context._coupling.ionic.prepare(
                            self.directions[:, index].reshape(-1, 3))
                         for index, name in enumerate(self.forces)}
        self.extra_rank = extra_rank
        self.extra_calls = []
        self.loading_seconds = perf_counter()-started
        # The original complete configuration is the normalization authority.
        original = json.loads(authenticate(self.filter_source_result(filter_packet)).read_bytes())
        self.mesh = original['configuration']['coefficient_specification']['mesh']
        self.ncells = int(np.prod(self.mesh))

    @staticmethod
    def filter_source_result(path):
        return json.loads(Path(path).read_bytes())['source']['original_result']

    def apply_extra(self, name, vectors, adjoint=False):
        if vectors.shape[0] != self.context.psi.shape[0] or vectors.shape[1] > self.extra_rank:
            raise ValueError('callback is restricted to admitted NEW extra columns')
        started = perf_counter()
        operator = self.context._hamiltonian if name == 'hamiltonian' else self.prepared[name]
        action = operator.adjoint_apply(vectors) if adjoint else operator.apply(vectors)
        if name != 'hamiltonian':
            action *= -1
        self.extra_calls.append(dict(name=name, adjoint=bool(adjoint), columns=vectors.shape[1],
                                     elapsed_seconds=perf_counter()-started))
        return action
