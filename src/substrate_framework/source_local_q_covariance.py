"""Two context-only cotangents, contracted by the original product formula."""
from pathlib import Path
import json

from .partitioned_q_storage import durable_publish
from .source_local_q_response import (SourceLocalQResponse, pullback_grid,
                                      prepare_pullback_frames)
from .source_local_q_contraction import source_contraction
from .streamed_material_state import file_binding


def _encoded(value):
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, dict):
        return {key: _encoded(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encoded(item) for item in value]
    return value


class SourceLocalTwoSidedCovariance:
    def __init__(self, preparation, executor):
        if preparation.response is not executor.r:
            raise ValueError('preparation and propagated q block differ')
        self.preparation, self.executor = preparation, executor

    def pair(self, left, right, times, second_times, directory,
             *, left_context=None, right_context=None):
        times, second_times = pullback_grid(times), pullback_grid(second_times)
        directory = Path(directory).resolve()
        frame_schedule = prepare_pullback_frames(self.executor, (times, second_times))
        states = []
        try:
            for side, probe, grid, binding in (
                    ('left', left, times, left_context),
                    ('right', right, second_times, right_context)):
                response = SourceLocalQResponse(self.executor, directory/side)
                state = (response.pullback(probe, grid) if binding is None
                         else response.restore(probe, grid, binding))
                states.append(state)
            result = source_contraction(self.preparation, *states, self.executor,
                                        directory/'contraction-work')
            result.update(time=float(times[-1]), second_time=float(second_times[-1]),
                dynamics=self.executor.dynamics, propagation_identity=self.executor.identity,
                execution_scheme='source-local-discrete-rk4-v1',
                left_context=states[0].binding, right_context=states[1].binding,
                mean_frame_times=frame_schedule, whole_spatial_reconstruction=False)
            record = dict(schema='source-local-q-pair-v1', complete=True,
                covariance=_encoded(result),
                restart='left_context/right_context are source-local pullbacks, not mapped Q checkpoints')
            path = directory/'complete.json'
            if path.exists():
                if json.loads(path.read_text()) != record:
                    raise ValueError('source-local pair publication differs')
            else:
                durable_publish(path, record, self.executor.budget)
            result['context_pair_binding'] = file_binding(path)
            self.executor.retire_frames_except((), result['context_pair_binding'])
            return result
        finally:
            for state in states:
                state.release()
