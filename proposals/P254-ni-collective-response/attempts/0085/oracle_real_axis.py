#!/usr/bin/env python3
"""Small same-discrete-operator oracle; NOT the old physical qualification suite.

Imports the admitted asymmetric noncommuting physical source/kernel fixtures.
Compares every packed R/< entry with the actual old A61 executable. The original
A54 executable source provides explicit wrong-order/drop-initial/drop-noise
mutants; they must fail the SAME finite-aware entrywise comparison predicate.
This tests the changed representation/operator, not continuum accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '0049'))
from finite_gaussian_reference import fermion_case, boson_case

OLD_SHA = 'c702cd49b7da7bc3a7d054988e3c116c148084d24092ba638154431339cbc064'
ATOL = 2e-11
RTOL = 2e-11


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def packed_input(path, times, case, statistics, order):
    nt, d = len(times)-1, case['source'][0].shape[-1]
    header = struct.pack('<8siiiid', b'P254VIE1', nt, d, order, statistics, 1.)
    with path.open('xb') as stream:
        stream.write(header)
        stream.write(np.asarray(times, dtype='<f8').tobytes())
        for component in (*case['source'], *case['kernel']):
            for n in range(nt+1):
                for m in range(n+1):
                    # R uses (n,m); < uses (m,n), exactly A49/A61's codec.
                    lesser = component is case['source'][1] or component is case['kernel'][1]
                    matrix = component[m,n] if lesser else component[n,m]
                    stream.write(np.asarray(matrix, dtype='<c16').tobytes(order='C'))
    return header + np.asarray(times, dtype='<f8').tobytes()


def result_entries(path, prefix, nt, dimension, extra_components=False):
    raw = path.read_bytes()
    entries = (nt+1)*(nt+2)//2*dimension*dimension
    expected = len(prefix) + entries*16*(6 if extra_components else 2)
    if len(raw) != expected or raw[:len(prefix)] != prefix:
        raise ValueError(f'{path}: output extent/header/times ABI mismatch')
    result = np.frombuffer(raw, dtype='<c16', offset=len(prefix), count=2*entries).copy()
    if not np.isfinite(result).all():
        raise ValueError(f'{path}: nonfinite R/< result (NaN cannot pass)')
    return result.reshape(2, entries)


def compare(actual, expected):
    if actual.shape != expected.shape or not np.isfinite(actual).all() or not np.isfinite(expected).all():
        raise ValueError('comparison requires same shape and every entry finite')
    delta = np.abs(actual-expected)
    allowed = ATOL + RTOL*np.abs(expected)
    return dict(accepted=bool(np.all(delta <= allowed)),
                retarded_max_abs=float(delta[0].max()),
                lesser_max_abs=float(delta[1].max()),
                max_tolerance_ratio=float(np.max(delta/allowed)),
                compared_complex_entries=int(actual.size))


def run(argv, work):
    completed = subprocess.run([str(x) for x in argv], check=True, capture_output=True,
                               text=True, timeout=60)
    # Persist original program output, including failed scientific predicates
    # evaluated afterward; no actual-Gamma input is ever consumed here.
    (work / (Path(argv[2]).name + '.stdout')).write_text(completed.stdout)
    (work / (Path(argv[2]).name + '.stderr')).write_text(completed.stderr)


def rejection_checks(binary, source, native, prefix_bytes, work):
    """Small unchanged-ABI rejection checks, never a material input replay."""
    outputs = {}
    old_result = binding(native)
    truncated = work/'truncated-input.bin'
    truncated.write_bytes(source.read_bytes()[:-1])
    nonfinite = work/'nonfinite-input.bin'
    raw = bytearray(source.read_bytes())
    raw[len(prefix_bytes):len(prefix_bytes)+8] = struct.pack('<d', float('nan'))
    nonfinite.write_bytes(raw)
    for name, argv in [
            ('existing-result', [binary,source,native]),
            ('caller-byte-cap', [binary,source,work/'rejected-cap.bin',source.stat().st_size-1]),
            ('truncated-extent', [binary,truncated,work/'rejected-truncated.bin']),
            ('nonfinite-input', [binary,nonfinite,work/'rejected-nonfinite.bin'])]:
        completed = subprocess.run([str(x) for x in argv], capture_output=True,
                                   text=True, timeout=60)
        outputs[name] = dict(returncode=completed.returncode, stderr=completed.stderr)
        if completed.returncode != 2:
            raise AssertionError(f'{name}: expected old-ABI rejection, got {completed.returncode}')
        if name != 'existing-result' and Path(argv[2]).exists():
            raise AssertionError(f'{name}: rejection must precede result creation')
    if binding(native) != old_result:
        raise AssertionError('existing-result rejection mutated the persisted valid result')
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--old-binary', type=Path,
                        default=Path('/var/tmp/p254-a61-lean-002/build/solve_result'))
    parser.add_argument('--mutant-binary', required=True, type=Path)
    parser.add_argument('--workdir', required=True, type=Path)
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=False)
    args.binary = args.binary.resolve()
    args.old_binary = args.old_binary.resolve()
    args.mutant_binary = args.mutant_binary.resolve()
    old_binding = binding(args.old_binary)
    if old_binding['sha256'] != OLD_SHA:
        raise ValueError('actual original A61 binary binding mismatch')
    sources = [HERE/'lean_result.cpp', HERE/'real_axis_keldysh.hpp', HERE/'CMakeLists.txt',
               Path(__file__).resolve(), HERE.parent/'0061'/'solve_result.cpp',
               HERE.parent/'0054'/'supplied_keldysh.hpp',
               HERE.parent/'0054'/'solve_components.cpp',
               HERE.parent/'0049'/'finite_gaussian_reference.py',
               HERE.parent/'0049'/'component_io.hpp',
               HERE.parent/'0049'/'supplied_vie2.hpp',
               HERE.parent/'0040'/'nessi_gaussian_ivp.hpp']
    prefix = Path('/home/administrator/.local/opt/nessi-6bb4c832')
    sources += [prefix/'include/cntr/cntr_convolution_impl.hpp',
                prefix/'include/cntr/cntr_vie2_impl.hpp', prefix/'lib/libcntr.so']
    report = dict(schema='p254-A85-finite-real-axis-equivalence-v1',
                  scope='Every stored R/< entry; same finite A54/A61 operator, not physical requalification or a speed claim.',
                  atol=ATOL, rtol=RTOL,
                  old_binary=old_binding, binary=binding(args.binary),
                  mutant_binary=binding(args.mutant_binary),
                  sources=[binding(p) for p in sources], cases={})
    report_path = args.workdir/'oracle-report.json'
    def save():
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    save()
    for label, factory, statistics in [('fermion', fermion_case, -1), ('boson_position', boson_case, 1)]:
        for order in range(1,6):
            # Both the literal minimal startup-only history and a history with
            # subsequent propagation are exercised for every admitted order.
            for nt in (order, 12):
                key = f'{label}-order{order}-nt{nt}'
                work = args.workdir/key
                work.mkdir()
                times = np.linspace(0., .6, nt+1)
                case = factory(times)
                qr, _ = case['source']
                kr, _ = case['kernel']
                asymmetric = float(np.max(np.abs(qr[-1,0]-qr[-1,0].T)))
                # Noncommutativity is measured on the actual admitted input,
                # rather than asserted from the fixture's descriptive name.
                q, k = qr[-1,0], kr[-1,0]
                commutator = float(np.max(np.abs(q@k-k@q)))
                if not np.isfinite(commutator) or commutator <= 1e-12:
                    raise AssertionError(f'{key}: fixture fails noncommuting-input gate')
                source = work/'input.bin'
                prefix_bytes = packed_input(source,times,case,statistics,order)
                native, old = work/'native.bin', work/'old.bin'
                run([args.binary,source,native,source.stat().st_size],work)
                run([args.old_binary,source,old,source.stat().st_size],work)
                actual = result_entries(native,prefix_bytes,nt,2)
                expected = result_entries(old,prefix_bytes,nt,2)
                verdict = compare(actual,expected)
                report['cases'][key] = dict(input=binding(source), native=binding(native),
                                          old=binding(old), noncommutator_max=commutator,
                                          source_ret_transpose_asymmetry_max=asymmetric,
                                          comparison=verdict, mutations={})
                save()
                if not verdict['accepted']:
                    raise AssertionError(f'{key}: changed-operator mismatch: {verdict}')
                if label == 'fermion' and order == 1 and nt == 1:
                    report['abi_rejections'] = rejection_checks(
                        args.binary,source,native,prefix_bytes,work)
                    save()
                # One exposing nt12 history per statistic/order is sufficient;
                # do not replay the broad physical proof/refinement suite.
                if nt == 12:
                    for omission,name in [(1,'drop-initial'),(2,'drop-noise'),(3,'wrong-order')]:
                        mutant = work/f'{name}.bin'
                        run([args.mutant_binary,source,mutant,omission],work)
                        mutated = result_entries(mutant,prefix_bytes,nt,2,extra_components=True)
                        sensitivity = compare(mutated,expected)
                        report['cases'][key]['mutations'][name] = dict(
                            output=binding(mutant), comparison=sensitivity)
                        save()
                        if sensitivity['accepted'] or sensitivity['max_tolerance_ratio'] <= 100.:
                            raise AssertionError(f'{key}: {name} not decisively exposed: {sensitivity}')
    report['verdict'] = 'PASS'
    save()
    print('A85_FINITE_REAL_AXIS_ORACLE_COMPLETE cases='+str(len(report['cases'])), flush=True)


if __name__ == '__main__':
    main()
