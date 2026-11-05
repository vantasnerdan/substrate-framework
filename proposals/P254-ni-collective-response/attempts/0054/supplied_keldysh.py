"""Structure-preserving fixed supplied Dyson boundary; reuse A49 strict transport.

The saved schema stays A49-compatible; algorithm identity and raw diagnostics are
bound separately. No thermal initialization, inverse G0, tolerance change or repair.
"""
from pathlib import Path
import sys
import json
import numpy as np

A49 = Path(__file__).resolve().parent.parent / '0049'
if str(A49) not in sys.path:
    sys.path.append(str(A49))
import supplied_vie2 as transport
from command_receipts import binding, run_recorded


def solve_supplied(times, source_ret, source_less, kernel_ret, kernel_less, *,
                   statistics, order, binary, workdir, beta=1., omission=0):
    """Solve a complete fixed history; omission=1/2/3 is explicit oracle-only.

    Returned NPZ uses the unchanged strict A49 reader. F/Fcc are independently
    computed original -G0*Sigma/-Sigma*G0, not interchangeable companions.
    """
    if omission not in (0, 1, 2, 3):
        raise ValueError('unsupported explicit diagnostic omission')
    arrays = [np.asarray(a) for a in (source_ret, source_less, kernel_ret, kernel_less)]
    times, dimension = transport._validate(times, arrays, statistics, order, beta)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=False)
    raw, packed, output = (workdir / p for p in ('input.npz', 'input.bin', 'output.bin'))
    np.savez(raw, times=times, source_ret=arrays[0], source_less=arrays[1],
             kernel_ret=arrays[2], kernel_less=arrays[3], statistics=statistics, order=order, beta=beta)
    header = transport._HEADER.pack(b'P254VIE1', len(times)-1, dimension, order, statistics, beta)
    with packed.open('wb') as stream:
        stream.write(header)
        stream.write(times.tobytes())
        for index, array in enumerate(arrays):
            transport._write_triangle(stream, array, bool(index % 2))
    manifest = {'algorithm': 'p254-a54-two-sided-keldysh-v1', 'omission': omission,
                'input': binding(raw), 'transport': binding(packed), 'binary': binding(binary),
                'statistics': statistics, 'order': order, 'reader': binding(A49 / 'supplied_vie2.py')}
    receipt = workdir / 'manifest.json'
    receipt.write_text(json.dumps(manifest, indent=2)+'\n')
    argv = [str(Path(binary).resolve()), str(packed.resolve()), str(output.resolve())]
    if omission:
        argv.append(str(omission))
    run_recorded(argv, workdir, 'native', 120)
    # Native status/streams have been persisted before decode or strict consumer.
    with output.open('rb') as stream:
        if stream.read(transport._HEADER.size) != header:
            raise ValueError('solver output header mismatch')
        restored = np.frombuffer(stream.read(times.size*8), dtype='<f8')
        if not np.array_equal(restored, times):
            raise ValueError('solver output time mismatch')
        names = ('ret', 'less', 'f_ret', 'f_less', 'fcc_ret', 'fcc_less')
        components = {name: transport._read_triangle(stream, times.size, dimension, bool(n % 2))
                      for n, name in enumerate(names)}
        if stream.read(1):
            raise ValueError('trailing solver output bytes')
    # Products F and Fcc are a general companion pair, not self-Hermitian.
    # Only the result lesser uses its own exchange relation. No averaging.
    for n in range(times.size):
        for m in range(n):
            components['f_less'][n,m] = -components['fcc_less'][m,n].conj().T
            components['fcc_less'][n,m] = -components['f_less'][m,n].conj().T
    solution = workdir / 'solution.npz'
    np.savez(solution, schema=transport._SCHEMA, times=times, statistics=statistics,
             order=order, beta=beta, **components)
    manifest.update(output=binding(output), solution=binding(solution),
                    diagnostics=binding(str(output)+'.diagnostics.json'))
    receipt.write_text(json.dumps(manifest, indent=2)+'\n')
    return solution


def load_solution(path):
    """Unchanged strict saved-numeric consumer, including raw diagonal exchange."""
    return transport.load_solution(path)
