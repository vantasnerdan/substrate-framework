"""Bounded original-column propagation in the unchanged partitioned midpoint.

The previously exercised Arnoldi spaces/envelopes are reused. Physical factors
stay full-PW; no saved-band closure, phase reset, clipping or scalar change.
"""
from dataclasses import dataclass
from time import perf_counter
from uuid import uuid4

import numpy as np
from run_column_krylov_material import ambient_growth, operator_norm
from substrate_framework.column_frame_krylov import ColumnFrameKrylov
from substrate_framework.material_physical_operators import PhysicalFrozenAction
from substrate_framework.partitioned_joint_mean import PartitionedJointMean
from substrate_framework.partitioned_joint_slots import digest, stage_identity


@dataclass(frozen=True)
class ColumnLimits:
    max_dimension: int
    workspace_bytes: int
    absolute_error_target: float

    def __post_init__(self):
        if (type(self.max_dimension) is not int or self.max_dimension < 1
                or type(self.workspace_bytes) is not int or self.workspace_bytes < 1
                or not np.isfinite(self.absolute_error_target) or self.absolute_error_target <= 0):
            raise ValueError('positive explicit column dimension/workspace/error allocation required')


def field_norm_bound(field):
    """Triangle/Frobenius bound for both actual directional nonlocal legs.

    It also bounds the anti-Hermitian part conservatively. No field cutoff or
    assumed exact floating Hermiticity is needed. Ordinary floating evaluation.
    """
    if field is None:
        return 0.
    value = float(np.max(np.abs(2*(field.momenta @ field.direction))))
    for species, derivative in field.terms:
        value += float(2*np.linalg.norm(derivative)*np.linalg.norm(species.beta)
                       *np.linalg.norm(species.active_d, 2))
    return value


def bounded_columns(action, frame, interval, *, limits, norm_bound, growth_bound, identity):
    started = perf_counter()
    column = ColumnFrameKrylov(action, frame, max_dimension=limits.max_dimension,
        action_cap=limits.max_dimension*frame.shape[1], max_memory_bytes=limits.workspace_bytes,
        ambient_growth_bound=growth_bound, operator_norm_bound=norm_bound, operator_identity=identity)
    progress = []
    for generation in range(limits.max_dimension):
        evolved, envelope = column.propagate(interval)
        progress.append(dict(generation=generation,
            estimated_error=envelope['absolute_error_envelope_estimate'], seconds=perf_counter()-started))
        if envelope['absolute_error_envelope_estimate'] <= limits.absolute_error_target:
            return evolved, dict(envelope=envelope, stages=progress, seed=column.seed_receipt,
                dimensions=column.sizes.tolist(), calls=column.calls, memory_plan=column.memory_plan,
                seconds=perf_counter()-started, qualification='floating analytic envelope; not interval certified')
        if generation+1 == limits.max_dimension or column.expand()['accepted_columns'] == 0:
            break
    raise RuntimeError('original-column bound allocation exhausted without reaching requested error estimate: '+repr(progress))


class ColumnPartitionedMean(PartitionedJointMean):
    @classmethod
    def from_serial(cls, serial, *, directory, cpus, source_identity, limits):
        if not isinstance(limits, ColumnLimits):
            raise TypeError('explicit immutable ColumnLimits required')
        result = super().from_serial(serial, directory=directory, cpus=cpus, source_identity=source_identity)
        result.column_limits = limits
        result.column_receipts = result.executor.directory.parent/'column-propagations'
        result.column_receipts.mkdir(parents=True, exist_ok=True)
        result.exponential = bounded_columns
        result._column_stage = None
        return result

    def propagate(self, coupling, h0, vertex, slot, stage, psi, interval):
        frozen = PhysicalFrozenAction(h0, coupling, potential=stage.potential[:, slot.component],
                                      direction=stage.u, shift=stage.shift)
        norm, _ = operator_norm(frozen.frozen)
        growth, _ = ambient_growth(frozen.frozen)
        field_bound = field_norm_bound(frozen.field)
        if self._column_stage is not stage:
            self._column_stage = stage
            self._column_stage_identity = stage_identity(stage)
        identity = digest(dict(model=self.identity, slot=slot.metadata(), frame=slot.frame,
            stage=self._column_stage_identity, interval=interval, limits=self.column_limits.__dict__))
        evolved, receipt = self.exponential(self.counted(frozen.apply), psi, interval,
            limits=self.column_limits, norm_bound=norm+field_bound, growth_bound=growth+field_bound,
            identity=identity)
        self.exponentials += 1
        self.budget.json(self.column_receipts/(identity+'-'+uuid4().hex+'.json'),
            dict(operator_identity=identity, slot=slot.identity, interval=interval,
                 limits=self.column_limits.__dict__, **receipt))
        return evolved
