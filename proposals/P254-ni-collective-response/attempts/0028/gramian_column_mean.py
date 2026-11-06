"""Original independent-column Arnoldi with source-aware Duhamel evaluation.

The basis expansion, full-row operator, phase shifts, raw column ordering and
output are unchanged. Only the residual integral estimator changes. Existing
historical controllers/proofs remain immutable; current consumers use this one.
"""
from math import hypot, sqrt
from time import perf_counter

import numpy as np
from scipy.linalg import expm

from gramian_krylov_bounds import gramian_envelope
from substrate_framework.column_frame_krylov import ColumnFrameKrylov


def bounded_source_columns(action, frame, interval, *, limits, norm_bound, growth_bound, identity):
    started = perf_counter()
    column = ColumnFrameKrylov(action, frame, max_dimension=limits.max_dimension,
        action_cap=limits.max_dimension*frame.shape[1], max_memory_bytes=limits.workspace_bytes,
        ambient_growth_bound=growth_bound, operator_norm_bound=norm_bound, operator_identity=identity)
    progress = []
    for generation in range(limits.max_dimension):
        output = np.zeros(frame.shape, dtype=np.complex128)
        records = []
        for j, size in enumerate(column.sizes):
            if not size:
                records.append(dict(column=j, absolute_error_envelope_estimate=0.))
                continue
            q, projected, residual = column._decomposition(j)
            value = column.beta[j]*expm(-1j*interval*projected)[:,0]
            output[:,j] = np.exp(-1j*interval*column.shifts[j])*(q@value)
            record = gramian_envelope(projected, residual, column.beta[j],
                column.initial_errors[j], interval, ambient_growth=growth_bound,
                seed_action_defect=column.seed_defects[j],
                target=limits.absolute_error_target/sqrt(frame.shape[1]))
            records.append(dict(column=j, **record))
        estimate = hypot(*(record['absolute_error_envelope_estimate'] for record in records))
        progress.append(dict(generation=generation, estimated_error=estimate,
                             seconds=perf_counter()-started))
        if estimate <= limits.absolute_error_target:
            return output, dict(envelope=dict(absolute_error_envelope_estimate=estimate,
                columns=records, ambient_growth_bound=growth_bound,
                qualification='source-aware numerical Duhamel estimate; not interval certified'),
                stages=progress, seed=column.seed_receipt, dimensions=column.sizes.tolist(),
                calls=column.calls, memory_plan=column.memory_plan,
                seconds=perf_counter()-started,
                qualification='same original full-row Arnoldi output; changed residual integration only')
        if generation+1 == limits.max_dimension or column.expand()['accepted_columns'] == 0:
            break
    raise RuntimeError('source-aware original-column allocation exhausted: '+repr(progress))
