"""Memory-bounded grouping of independent original-column Arnoldi actions.

All full PW rows and every original column remain, in original order. The
concatenated Frobenius error is the Euclidean norm of the panel errors; a panel
of width w among n columns receives epsilon*sqrt(w/n). These are the existing
floating analytic envelopes, not outward-rounded enclosures.
"""
from dataclasses import replace
from math import hypot, sqrt

import numpy as np

from joint_column_mean import bounded_columns
from substrate_framework.column_frame_krylov import memory_bound_bytes


def panel_plan(shape, limits):
    if len(shape) != 2 or min(shape) < 1:
        raise ValueError('nonempty full-row frame required')
    rows, columns = shape
    whole = memory_bound_bytes(rows, columns, limits.max_dimension)
    if whole <= limits.workspace_bytes:
        return dict(columns_per_panel=columns, panel_count=1,
            maximum_numerical_bytes=whole, additional_output_bytes=0)
    output = np.dtype(np.complex128).itemsize*rows*columns
    available = limits.workspace_bytes-output
    lo, hi, width = 1, columns, 0
    while lo <= hi:
        middle = (lo+hi)//2
        if memory_bound_bytes(rows, middle, limits.max_dimension) <= available:
            width, lo = middle, middle+1
        else:
            hi = middle-1
    if width == 0:
        raise MemoryError('workspace cannot hold complete output and one original-column Arnoldi panel')
    return dict(columns_per_panel=width, panel_count=(columns+width-1)//width,
        maximum_numerical_bytes=output+memory_bound_bytes(rows, width, limits.max_dimension),
        additional_output_bytes=output)


def bounded_panel_columns(action, frame, interval, *, limits, norm_bound, growth_bound, identity):
    plan = panel_plan(frame.shape, limits)
    if plan['panel_count'] == 1:
        value, receipt = bounded_columns(action, frame, interval, limits=limits,
            norm_bound=norm_bound, growth_bound=growth_bound, identity=identity)
        return value, dict(panel_plan=plan, **receipt)
    result = np.empty(frame.shape, dtype=np.complex128)
    records, estimates = [], []
    columns = frame.shape[1]
    for start in range(0, columns, plan['columns_per_panel']):
        stop = min(columns, start+plan['columns_per_panel'])
        local = replace(limits,
            workspace_bytes=limits.workspace_bytes-plan['additional_output_bytes'],
            absolute_error_target=limits.absolute_error_target*sqrt((stop-start)/columns))
        value, receipt = bounded_columns(action, frame[:, start:stop], interval,
            limits=local, norm_bound=norm_bound, growth_bound=growth_bound,
            identity=identity+':columns:'+str(start)+':'+str(stop))
        result[:, start:stop] = value
        estimates.append(receipt['envelope']['absolute_error_envelope_estimate'])
        records.append(dict(start=start, stop=stop, absolute_error_target=local.absolute_error_target,
                            **receipt))
        del value
    combined = hypot(*estimates)
    if combined > limits.absolute_error_target:
        raise RuntimeError('combined original-column panel estimate exceeds allocation')
    return result, dict(panel_plan=plan, panels=records,
        absolute_error_envelope_estimate=combined,
        qualification='concatenated floating analytic envelope estimates; no band/rank cut or certification')
