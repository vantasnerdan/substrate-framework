"""Complete source-order reduction, kernel publication and safe moment retirement."""
import json

import numpy as np

from .partitioned_q_storage import (bounded_array_identity, durable_publish,
                                   release_mapped_pages)
from .partitioned_q_work import arrays
from .source_local_q_state import small_state, encode_small
from .streamed_material_state import file_binding


def complete_context(executor, work, path, request, function, small):
    budget, r = executor.budget, executor.r
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        record = json.loads(path.read_text())
        if (record.get('schema') != 'source-local-q-context-v1'
                or record.get('complete') is not True or record['request'] != request):
            raise ValueError('source-local context request differs')
    else:
        records = work.execute('source-local-dual-moments', request, len(r.edges), function)
        moment_directory = work.last[1]
        n = np.zeros(r.density_shape, complex)
        g, ell = np.zeros(r.ndof, complex), np.zeros(3, complex)
        for item in records:
            value = arrays(item)
            n += value['n']
            g += value['g']
            ell += value['ell']
            for array in value.values():
                release_mapped_pages(array)
        kernel = None
        if executor.dynamics == 'screened-joint':
            kernel = work.execute('source-local-complete-kernel',
                dict(context=request, density=bounded_array_identity(n)), 1,
                lambda index: dict(potential=executor.ka(n, q_bohr_inv=r.q_bohr_inv)))[0]
        derivative = executor.small_rhs(small, small_state(r), g, ell, dual=True)
        record = dict(schema='source-local-q-context-v1', complete=True, request=request,
            stage_time=request['stage_time'], input_small=encode_small(small),
            derivative_small=encode_small(derivative), kernel=kernel,
            moment_directory=str(moment_directory))
        durable_publish(path, record, budget)
    binding = file_binding(path)
    potential = None
    if record['kernel'] is not None:
        for array in record['kernel']['arrays'].values():
            work.verify(array)
        potential = arrays(record['kernel'])['potential']
    work.retire(record['moment_directory'], binding)
    return dict(record=record, binding=binding, potential=potential)
