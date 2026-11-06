"""Real Cartesian/fermion supplied contour-Hermitian matrix VIE2 transport.

No Gaussian production, thermal initialization, q/-q packing or interpolation.
Full original arrays are retained in input.npz. Unsupported independent lesser
legs are rejected rather than identified, averaged or silently discarded.
"""
from pathlib import Path
import json
import struct
import numpy as np
from command_receipts import binding, run_recorded

_HEADER = struct.Struct('<8siiiid')
_SCHEMA = 'p254-supplied-vie2-v1'


def _validate(times, arrays, statistics, order, beta):
    times = np.asarray(times, dtype='<f8')
    if (times.ndim != 1 or times.size < 2 or order not in range(1, 6)
            or times.size <= order or statistics not in (-1, 1)
            or not np.isfinite(beta) or beta <= 0):
        raise ValueError('invalid times/order/statistics/inert beta')
    dt = times[1] - times[0]
    if (not np.isfinite(times).all() or times[0] != 0 or dt <= 0
            or not np.allclose(times, np.arange(times.size)*dt, rtol=1e-12, atol=1e-12)):
        raise ValueError('unsupported nonuniform time history; no interpolation')
    shape = np.shape(arrays[0])
    if len(shape) != 4 or shape[:2] != (times.size, times.size) or shape[2] != shape[3] or shape[2] < 1:
        raise ValueError('expected full arrays [time,time,row,column]')
    for i, array in enumerate(arrays):
        if np.shape(array) != shape or not np.isfinite(array).all():
            raise ValueError('invalid component shape or nonfinite entries')
        if i % 2:
            defect = np.max(np.abs(array + array.swapaxes(0, 1).swapaxes(2, 3).conj()))
            if defect > 1e-10 * (1 + np.max(np.abs(array))):
                raise ValueError('independent non-contour-Hermitian lesser legs unsupported; no repair')
        elif np.any(np.asarray(array)[np.triu_indices(times.size, 1)] != 0):
            raise ValueError('R input must be causal, upper triangle exactly zero')
    return times, shape[2]


def _write_triangle(stream, array, lesser):
    for n in range(array.shape[0]):
        # Lesser stored as (m,n); retarded stored as (n,m).
        block = array[:n+1, n] if lesser else array[n, :n+1]
        stream.write(np.asarray(block, dtype='<c16', order='C').tobytes())


def _read_triangle(stream, count, dimension, lesser):
    array = np.zeros((count, count, dimension, dimension), dtype=np.complex128)
    for n in range(count):
        size = (n+1)*dimension*dimension
        data = stream.read(size*16)
        if len(data) != size*16:
            raise ValueError('truncated solver output')
        block = np.frombuffer(data, dtype='<c16').reshape(n+1, dimension, dimension)
        if lesser:
            array[:n+1, n] = block
            for m in range(n):
                array[n, m] = -block[m].conj().T
        else:
            array[n, :n+1] = block
    return array


def solve_supplied(times, source_ret, source_less, kernel_ret, kernel_less, *,
                   statistics, order, binary, workdir, beta=1.):
    """Return saved solution path; workdir is a fresh runtime directory outside Git.

    Coupled callers refresh supplied kernels themselves; each call solves their
    linear fixed-kernel VIE, not a selfconsistent iteration. The C++ Solver also
    exposes refresh/start/step for caller-owned incremental contour feedback.
    """
    arrays = [np.asarray(a) for a in (source_ret, source_less, kernel_ret, kernel_less)]
    times, dimension = _validate(times, arrays, statistics, order, beta)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=False)
    raw = workdir / 'input.npz'
    np.savez(raw, times=times, source_ret=arrays[0], source_less=arrays[1],
             kernel_ret=arrays[2], kernel_less=arrays[3], statistics=statistics, order=order, beta=beta)
    input_path, output_path = workdir / 'input.bin', workdir / 'output.bin'
    with input_path.open('wb') as stream:
        stream.write(_HEADER.pack(b'P254VIE1', len(times)-1, dimension, order, statistics, beta))
        stream.write(times.tobytes())
        for i, array in enumerate(arrays):
            _write_triangle(stream, array, bool(i % 2))
    # Input and implementation bindings are durable before the command.
    manifest = {'schema': _SCHEMA, 'input': binding(raw), 'transport': binding(input_path),
                'binary': binding(binary), 'statistics': statistics, 'order': order,
                'scope': 'single contour-Hermitian source/solution; real Cartesian POSITION or fermion',
                'diagonal_policy': 'copy supplied source R diagonal verbatim; no canonical replacement',
                'source_ret_diagonal_max': float(max(np.max(np.abs(arrays[0][n,n])) for n in range(len(times))))}
    (workdir / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    run_recorded([str(Path(binary).resolve()), str(input_path.resolve()), str(output_path.resolve())], workdir, 'native')
    with output_path.open('rb') as stream:
        header = stream.read(_HEADER.size)
        if header != _HEADER.pack(b'P254VIE1', len(times)-1, dimension, order, statistics, beta):
            raise ValueError('solver output header mismatch')
        restored_times = np.frombuffer(stream.read(len(times)*8), dtype='<f8')
        if not np.array_equal(times, restored_times):
            raise ValueError('solver output time mismatch')
        names = ('ret', 'less', 'f_ret', 'f_less', 'fcc_ret', 'fcc_less')
        output = {name: _read_triangle(stream, len(times), dimension, bool(i % 2))
                  for i, name in enumerate(names)}
        if stream.read(1):
            raise ValueError('trailing solver output bytes')
    solution = workdir / 'solution.npz'
    np.savez(solution, schema=_SCHEMA, times=times, statistics=statistics,
             order=order, beta=beta, **output)
    manifest.update(output=binding(output_path), solution=binding(solution))
    (workdir / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return solution


def load_solution(path):
    """Consume actual persisted numeric results; causal R, full reconstructed <."""
    with np.load(path, allow_pickle=False) as saved:
        if saved['schema'].item() != _SCHEMA:
            raise ValueError('unsupported saved solution schema')
        result = {key: saved[key] for key in saved.files if key != 'schema'}
    _validate(result['times'], [result['ret'], result['less']], int(result['statistics']),
              int(result['order']), float(result['beta']))
    return result
