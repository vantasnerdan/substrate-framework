"""Restore the original source0 context and reuse all six paid raw action legs.

Only the explicit new692-column bank reaches physical operators. This loader
never instantiates A31 CachedBathInputs or evaluates H/F on the original346.
"""
from pathlib import Path
import json
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(HERE), str(ROOT/'src')] + [
    str(HERE.parent/name) for name in ('0031', '0030', '0029', '0028', '0034')]

import numpy as np
from clamped_scientific_readout import encode
from export_output_directed_actions import restore_filter
from original_noise_source import _identity, _same_binding
from packet_hamiltonian_actions import restore_packet_hamiltonian
from substrate_framework.streamed_material_state import authenticate, file_binding

CONTRACT = dict(path=str(HERE/'implementation-contract-001.json'), bytes=9415,
    sha256='037bfd5014b6f9ce994a0658af5f0e8d4b36031fd5c3fb92238ced5a197a61f1')
SOURCE_COMPATIBILITY = dict(
    path=str(HERE/'cached-vs-restored-source-compatibility-001.json'), bytes=2860,
    sha256='9f296a8a2de87bb9b6ea8ea2dc2d03913a21f5d71c0b9a228f52a45d4cac9500')
EXTRA_RANK = 692
CALL_ORDER = tuple((name, adjoint) for name in ('hamiltonian', 'mode0', 'mode1')
                   for adjoint in (False, True))


def load_contract():
    """Authenticate the settled preregistration, not an output-derived design."""
    return json.loads(authenticate(CONTRACT).read_bytes())


class ExpandedBathInputs:
    """Authenticated paid factors and counted actions on one new selected bank."""
    def __init__(self, action_packet, mode1_packet, filter_packet, extra_rank,
                 progress=None):
        if type(extra_rank) is not int or extra_rank != EXTRA_RANK:
            raise ValueError('A37 preregistration requires exactly692 extra columns')
        started = perf_counter()
        self.extra_rank = extra_rank
        self.extra_calls = []
        self._extra_vectors = None
        self.progress = progress or (lambda event: print(
            '[physics] '+json.dumps(encode(event), sort_keys=True, allow_nan=False), flush=True))
        self.contract = load_contract()
        self.settled = {Path(item['path']).resolve(): item
                        for item in self.contract['settled_inputs']}
        # Freeze the reused algorithms before loading or building physical inputs.
        for name in ('cached_bath_inputs.py', 'noise_subspace.py',
                     'filtered_quadratic_bath.py', 'export_frozen_bath_comparator.py'):
            authenticate(self.settled[(HERE.parent/'0031'/name).resolve()])
        self.action_binding, packet = self._packet(action_packet)
        self.mode1_binding, second = self._packet(mode1_packet)
        self.filter_binding, filter_packet_data = self._packet(filter_packet)
        self._check_packets(packet, second, filter_packet_data)
        self.filter, restored_binding, self.filter_archive = restore_filter(filter_packet)
        if not _same_binding(restored_binding, self.filter_binding):
            raise ValueError('filter changed during restoration')
        self._bound(self.filter_archive)
        self._load_channels(packet, second)
        self.source_compatibility = json.loads(authenticate(SOURCE_COMPATIBILITY).read_bytes())
        authenticate(self.source_compatibility['current_restoration'])
        if not _same_binding(packet['source']['hamiltonian_implementation'],
                             self.source_compatibility['original_computation']):
            raise ValueError('cached computation implementation differs from admitted original')
        self.original_source = packet['source']
        self._emit(stage='restore_original_context', status='started')
        self.context = restore_packet_hamiltonian(packet['source']['packet']['path'])
        self._check_context(packet, second)
        self._emit(stage='restore_original_context', status='complete',
                   elapsed_seconds=perf_counter()-started)
        self._load_actions(packet, second)
        self._prepare_directions(packet, second)
        original_binding = filter_packet_data['source']['original_result']
        original = json.loads(self._bound(original_binding).read_bytes())
        self.mesh = original['configuration']['coefficient_specification']['mesh']
        if (len(self.mesh) != 3 or any(type(n) is not int or n <= 0 for n in self.mesh)
                or self.mesh != [4, 4, 4]):
            raise ValueError('actual canonical filter original4x4x4 mesh required')
        self.ncells = int(np.prod(self.mesh))
        self.original_result_binding = original_binding
        self.packet = packet
        self.mode1_packet = second
        self.loading_seconds = perf_counter()-started

    def _emit(self, **event):
        self.progress(event)

    def _check_binding(self, binding):
        expected = self.settled.get(Path(binding['path']).resolve())
        if expected is None or not _same_binding(binding, expected):
            raise ValueError('input is not an exact preregistered binding: '+binding['path'])
        return expected

    def _bound(self, binding):
        return authenticate(self._check_binding(binding))

    def _packet(self, path):
        expected = self.settled.get(Path(path).resolve())
        if expected is None:
            raise ValueError('canonical contract-bound packet path required: '+str(path))
        return dict(expected), json.loads(authenticate(expected).read_bytes())

    def _check_packets(self, packet, second, filter_data):
        if (packet['schema'] != 'P254-output-directed-original-sector-actions-v1'
                or packet['status'] != 'complete' or second['status'] != 'complete'
                or filter_data['schema'] != 'P254-canonical-force-filter-material-v1'
                or filter_data['status'] != 'complete' or packet['quadrature_applied'] is not False
                or packet['source']['original_record'] != 0
                or not _same_binding(second['original_action_packet'], self.action_binding)
                or not _same_binding(packet['filter_packet'], self.filter_binding)
                or not _same_binding(packet['filter_archive'], filter_data['archive'])
                or second['source'] != packet['source']):
            raise ValueError('complete original source0/action/mode/filter lineage differs')
        self._bound(packet['source']['packet'])
        self._bound(packet['channel_archive'])
        self._bound(filter_data['archive'])

    def _load_channels(self, packet, second):
        with np.load(authenticate(packet['channel_archive']), allow_pickle=False) as archive:
            self.directions = archive['ionic_basis'][:, :2]
            self.modal_coefficients = archive['ionic_modal_chebyshev'][:, :, :2, :]
            bounds = archive['segment_bounds']
        self.segment_bounds = self.filter.segment_bounds
        metadata = self.filter.metadata
        if (self.directions.shape != (metadata['ndof'], 2)
                or self.modal_coefficients.shape != (len(bounds), 8, 2, 1)
                or metadata['q_index'] != 0 or metadata['nprobes'] != 1
                or not np.array_equal(bounds, self.segment_bounds)
                or not all(np.isfinite(a).all() for a in
                           (self.directions, self.modal_coefficients, self.segment_bounds))):
            raise ValueError('actual canonical two-mode filter/channel layout differs')
        selected = packet['selected_physical_mode']
        direction0 = np.asarray([complex(v['real'], v['imag']) for v in selected['direction']])
        direction1 = np.asarray(second['direction_real'])+1j*np.asarray(second['direction_imag'])
        if (selected['kind'] != 'ionic' or selected['mode'] != 0
                or not np.array_equal(direction0, self.directions[:, 0])
                or not np.array_equal(direction1, self.directions[:, 1])):
            raise ValueError('paid physical force directions differ, including complex parts')

    def _check_context(self, packet, second):
        context = self.context
        self.restoration_source = encode(context.source)
        actual_source = dict(self.restoration_source)
        original_source = dict(self.original_source)
        actual_implementation = actual_source.pop('hamiltonian_implementation')
        original_source.pop('hamiltonian_implementation')
        if (not _same_binding(actual_implementation,
                              self.source_compatibility['current_restoration'])
                or context.psi.shape[1] != 346 or actual_source != original_source
                or context.state_authority != packet['state_authority']
                or context.factor_identity != packet['factor_identity']
                or context.pw_basis_identity != packet['pw_basis_identity']
                or not _same_binding(context.frame, packet['frame'])
                or context.quadrature_weight != packet['quadrature_weight']
                or encode(context.hamiltonian_specification) != packet['work']['hamiltonian']['operator_specification']):
            raise ValueError('restored original source/state/basis/literal H differs from paid actions')

    def _load_actions(self, packet, second):
        from cached_bath_inputs import read_action
        records = dict(hamiltonian=packet['actions']['literal_initial_hamiltonian'],
                       mode0=packet['actions']['ni_output_ionic_mode0'], mode1=second['actions'])
        specifications = dict(hamiltonian=packet['work']['hamiltonian']['operator_specification'],
            mode0=packet['work']['force']['operator_specification'],
            mode1=second['work']['operator_specification'])
        factors = {}
        for name, record in records.items():
            for leg in ('action', 'adjoint_action'):
                self._check_binding(record[leg])  # read_action authenticates the actual bytes.
            identity = _identity(specifications[name])
            if (record['operator_identity'] != identity
                    or record['adjoint_identity'] != identity+':actual-Hilbert-adjoint'
                    or record['pw_basis_identity'] != self.context.pw_basis_identity
                    or record['source_authority'] != self.context.state_authority
                    or record['units'] != ('Ry' if name == 'hamiltonian' else 'Ry/bohr')):
                raise ValueError('paid operator identity/state/units differs: '+name)
            factors[name] = read_action(record)
            factors[name].validate(self.context.psi.shape, self.context.pw_basis_identity)
        self.hamiltonian = factors.pop('hamiltonian')
        self.forces = factors
        self.cached_action_records = records
        self.operator_specifications = specifications

    def _prepare_directions(self, packet, second):
        self.prepared = {}
        self.direction_preparations = []
        for index, name in enumerate(self.forces):
            direction = self.directions[:, index].reshape(-1, 3)
            spec = self.operator_specifications[name]
            saved = spec['derivative_direction']
            saved_direction = np.asarray(saved['real'])+1j*np.asarray(saved['imaginary'])
            if (not np.array_equal(direction, saved_direction) or spec['force_sign'] != -1
                    or spec['derivative_scale'] != 1 or spec['quadrature_applied'] is not False
                    or spec['external_q'] != [0, 0, 0] or spec['time_inverse_ry'] != 0
                    or spec['basis_identity'] != self.context.pw_basis_identity
                    or spec['original_state_authority'] != self.context.state_authority
                    or spec['matched_model_identity'] != self.context.source['model_identity']
                    or spec['physical_sources'] != encode(list(self.context._coupling.factors.sources))
                    or not _same_binding(spec['source'], self.context.source['source'])
                    or spec['adjoint_api'] != 'BlochDirectionalAction.adjoint_apply'):
                raise ValueError('actual physical force specification differs: '+name)
            started = perf_counter()
            prepared = self.context._coupling.ionic.prepare(direction)
            if not np.array_equal(prepared.direction, direction):
                raise ValueError('physical preparation changed complex force direction: '+name)
            self.prepared[name] = prepared
            self.direction_preparations.append(dict(name=name,
                elapsed_seconds=perf_counter()-started, physical_actions=0,
                operator_identity=self.forces[name].operator_identity))

    def apply_extra(self, name, vectors, adjoint=False):
        """Apply exactly one prescribed independent leg to the SAME new692 bank."""
        index = len(self.extra_calls)
        if (index >= len(CALL_ORDER) or (name, adjoint) != CALL_ORDER[index]
                or vectors.shape != (self.context.psi.shape[0], self.extra_rank)
                or np.shares_memory(vectors, self.context.psi)):
            raise ValueError('only six ordered actions on new692 selected columns are admitted')
        if self._extra_vectors is None:
            self._extra_vectors = vectors
        elif vectors is not self._extra_vectors:
            raise ValueError('all six callbacks must use the identical selected extra bank')
        started = perf_counter()
        record = dict(name=name, adjoint=bool(adjoint), columns=vectors.shape[1],
            callback_index=index+1, status='started', operator_identity=(
                self.hamiltonian if name == 'hamiltonian' else self.forces[name]).operator_identity)
        self.extra_calls.append(record)
        self._emit(stage='new_physical_action', **record)
        try:
            operator = self.context._hamiltonian if name == 'hamiltonian' else self.prepared[name]
            action = operator.adjoint_apply(vectors) if adjoint else operator.apply(vectors)
            if name != 'hamiltonian':
                action *= -1
            record['status'] = 'complete'
            return action
        except BaseException as error:
            record.update(status='failed', error=repr(error))
            raise
        finally:
            record['elapsed_seconds'] = perf_counter()-started
            self._emit(stage='new_physical_action', **record)
