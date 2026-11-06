"""Resumable physical joint mean using the settled serial midpoint algorithm.

StreamedJointMean.step itself remains the solver: a state adapter replaces only
its final frame stream with independent saved per-slot outputs. No equations,
coefficients, convergence test, or frozen-stage selection are rederived here.
Completed boundaries/history retain the ordinary DiskJointMean interface.
"""
import json
from pathlib import Path
from uuid import uuid4

import numpy as np

from .partitioned_joint_slots import SlotExecutor, digest, publish, stage_identity
from .streamed_joint_mean import StreamedJointMean
from .streamed_joint_reduction import JointStage
from .streamed_joint_state import DiskJointMean
from .streamed_material_state import authenticate, file_binding


def save_stage(stage, directory, budget, *, metadata):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    identity = digest(dict(stage=stage_identity(stage), metadata=metadata))
    pointer = directory/'complete.json'
    if pointer.exists():
        saved = json.loads(pointer.read_text())
        if saved['identity'] != identity:
            raise ValueError('saved stage request changed')
        JointStage.restore(saved['record'])
        return saved['record']
    record = stage.save(directory/('arrays-'+uuid4().hex), budget, metadata=metadata)
    publish(pointer, dict(identity=identity, record=record), budget)
    return record


class _PartitionedSuccessor:
    def __init__(self, state, executor):
        self.state, self.executor = state, executor

    def __getattr__(self, name):
        return getattr(self.state, name)

    def successor(self, directory, frames, *, time, u, p, coordinate, momentum,
                  stage_record, evidence, budget):
        # The serial generator has not been entered. The same propagate() API
        # runs independently in workers with this exact saved frozen stage.
        stage = JointStage.restore(stage_record)
        key, work, coords, records = self.executor.execute(self.state, frozen_stage=stage,
            interval=stage_record['dt_ry_inverse'], operation='frames')
        slots = []
        for slot, record in zip(self.state.slots, records):
            frame = record['arrays']['frame']
            psi = np.load(frame['path'], mmap_mode='r', allow_pickle=False)
            if psi.shape != slot.shape or psi.dtype != np.complex128 or not np.isfinite(psi).all():
                raise ValueError('saved successor changed original raw frame')
            slots.append(dict(slot.metadata(), coefficient=slot.coefficient, frame=frame))
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=False)
        evidence.update(partitioned_frame_request=key, cpus=list(self.executor.cpus),
                        blas_threads=1, budget_at_publication=budget.snapshot(),
                        frames='authenticated completed per-slot files; original D/basis unchanged')
        DiskJointMean._publish(directory, slots, self.identity, time, u, p, coordinate, momentum,
            self.origin, self.history+[stage_record], evidence, budget)
        return DiskJointMean(directory/'boundary.json', expected_identity=self.identity, verify=False)


class PartitionedJointMean(StreamedJointMean):
    @classmethod
    def from_serial(cls, serial, *, directory, cpus, source_identity):
        result = cls(serial.reference, serial.factory, serial.coefficients, serial.field, serial.charges,
            model_identity=serial.identity, expected_slots=serial.expected_slots, grid_measure=serial.measure,
            budget=serial.budget, exponential=serial.exponential, reality_tolerance=serial.reality_tolerance)
        result.executor = SlotExecutor(result, directory, cpus=cpus, source_identity=source_identity)
        return result

    def reduce(self, state, *, frozen_stage=None, interval=0., coordinates=None):
        key, directory, coords, records = self.executor.execute(state, frozen_stage=frozen_stage,
            interval=interval, coordinates=coordinates)
        return self.executor.stage(directory, key, state, coords, records)

    def step(self, state, directory, *, dt, max_iterations, tolerance):
        # Reusing the settled solver also preserves its nonlinear failure receipt
        # and its precise tested-stage (rather than newly updated stage) choice.
        self.validate(state)
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        key, _, _ = self.executor.request(state)
        request = dict(input_request=key, dt=dt, max_iterations=max_iterations, tolerance=tolerance)
        manifest = directory/'request.json'
        if manifest.exists():
            if json.loads(manifest.read_text()) != request:
                raise ValueError('step directory belongs to a different state/solver request')
        else:
            publish(manifest, request, self.budget)
        complete = directory/'complete.json'
        if complete.exists():
            record = json.loads(complete.read_text())
            if record['request'] != digest(request):
                raise ValueError('completed step request changed')
            return DiskJointMean(authenticate(record['boundary']), expected_identity=self.identity)
        # A killed parent may have published a valid boundary but not the pointer.
        candidates = sorted(directory.glob('attempt-*/completed/boundary.json'))
        if candidates:
            result = DiskJointMean(candidates[0], expected_identity=self.identity)
        else:
            result = super().step(_PartitionedSuccessor(state, self.executor),
                directory/('attempt-'+uuid4().hex), dt=dt, max_iterations=max_iterations, tolerance=tolerance)
        publish(complete, dict(request=digest(request), boundary=file_binding(result.path)), self.budget)
        return result

    def evolve(self, state, directory, *, physical_interval, steps, max_steps,
               max_iterations, tolerance, refinement_label):
        if (type(steps) is not int or type(max_steps) is not int or not 1 <= steps <= max_steps
                or not np.isfinite(physical_interval) or physical_interval <= 0 or not refinement_label):
            raise ValueError('bounded steps, positive physical interval and refinement label required')
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        records = []
        for index in range(steps):
            state = self.step(state, directory/('step-%06d' % (index+1)), dt=physical_interval/steps,
                max_iterations=max_iterations, tolerance=tolerance)
            pointer = directory/('completed-%06d.json' % (index+1))
            value = dict(boundary=str(state.path), time_ry_inverse=state.time, requested_steps=steps,
                         completed_steps=index+1, refinement_label=refinement_label)
            if pointer.exists():
                if json.loads(pointer.read_text()) != value:
                    raise ValueError('completed boundary pointer changed')
            else:
                publish(pointer, value, self.budget)
            endpoint = self.reduce(state)
            record = save_stage(endpoint, directory/('endpoint-%06d' % (index+1)), self.budget,
                metadata=dict(kind='complete-endpoint-reduction', boundary=str(state.path),
                              time_ry_inverse=state.time, refinement_label=refinement_label))
            records.append(record)
        return state, records
