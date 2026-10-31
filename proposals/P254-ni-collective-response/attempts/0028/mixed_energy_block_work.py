"""Readout-only column grouping and authenticated completed scalar reuse.

The registered mean/dynamics block size is unchanged. Only each forked bare-H0
item temporarily uses the declared energy column grouping; the unchanged item
function owns every physical density dyad. Original QWork owns all publication,
worker supervision, cumulative action admissions and scalar barriers.
"""
import json

import numpy as np

from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import authenticate, file_binding


class EnergyBlockWork(QWork):
    def __init__(self, directory, budget, cpus, identity, *, dynamics,
                 column_block_size, prior_request, prior_receipts):
        super().__init__(directory, budget, cpus, identity)
        if type(column_block_size) is not int or column_block_size < 1:
            raise ValueError('positive explicit energy-only column grouping required')
        if dynamics.budget is not budget:
            raise ValueError('same original cumulative energy budget required')
        self.dynamics, self.column_block_size = dynamics, column_block_size
        self.prior_request = dict(prior_request)
        request = json.loads(authenticate(self.prior_request).read_text())
        if (request.get('schema') != 'partitioned-q-work-v2'
                or request.get('source_identity') != identity
                or request.get('operation') != 'native-mixed-energy-bare-h0'
                or request.get('count') != dynamics.count):
            raise ValueError('original full-source bare-energy request required')
        self.prior_input = request['input']
        self.donors = {}
        for binding in prior_receipts:
            record = json.loads(authenticate(binding).read_text())
            index = record.get('index')
            if (record.get('complete') is not True or record.get('request') != digest(request)
                    or type(index) is not int or not 0 <= index < dynamics.count
                    or index in self.donors or set(record.get('arrays', {})) != {'mixed_bare_terms'}):
                raise ValueError('original completed ordered bare-energy source required')
            array = record['arrays']['mixed_bare_terms']
            value = np.load(authenticate(array), allow_pickle=False)
            if value.shape != (4,) or value.dtype != np.complex128 or not np.isfinite(value).all():
                raise ValueError('four genuine original complex dyad scalars required')
            self.donors[index] = (dict(binding), array)

    def execute(self, operation, identity, count, function):
        if operation != 'native-mixed-energy-bare-h0':
            return super().execute(operation, identity, count, function)
        for name in ('input', 'prepared', 'mean_source'):
            if identity[name] != self.prior_input[name]:
                raise ValueError('retained scalar state/prepared/mean authority differs')
        old_energy = json.loads(authenticate(self.prior_input['request']).read_text())
        new_energy = json.loads(authenticate(identity['request']).read_text())
        if any(old_energy[name] != new_energy[name] for name in ('input', 'native_context', 'mean_ledger', 'work_identity', 'scalar')):
            raise ValueError('original scalar physical energy request differs')
        for role, binding in old_energy['sources'].items():
            authenticate(binding)
            if role != 'substrate_framework.partitioned_q_work' and new_energy['sources'].get(role) != binding:
                raise ValueError('original scalar physical implementation differs: '+role)
        import substrate_framework.partitioned_q_work as original_work
        if old_energy['sources']['substrate_framework.partitioned_q_work'] != file_binding(original_work.__file__):
            raise ValueError('original QWork implementation changed')
        execution = dict(identity, energy_only_column_block_size=self.column_block_size,
            retained_original_scalar_receipts=[self.donors[i][0] for i in sorted(self.donors)],
            original_bare_request=self.prior_request)

        def evaluate(index):
            if index in self.donors:
                return dict(mixed_bare_terms=np.load(authenticate(self.donors[index][1]), allow_pickle=False))
            original = self.dynamics.block_size
            self.dynamics.block_size = self.column_block_size
            try:
                return function(index)
            finally:
                self.dynamics.block_size = original

        return super().execute(operation, execution, count, evaluate)
