#!/usr/bin/env python3
"""Finite read-only material postprocessing; no Native computation or mode claim."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback
import numpy as np
from substrate_framework.qe_material_frozen_source import load_material_frozen
from substrate_framework.qe_frozen_curvature import verify_buffer

HERE = Path(__file__).resolve().parent


def bind(path, raw=None):
    raw = path.read_bytes() if raw is None else raw
    return {'path': str(path), 'bytes': str(len(raw)), 'sha256': hashlib.sha256(raw).hexdigest()}


def implementation():
    files = {str(Path(__file__).resolve())}
    for name, module in list(sys.modules.items()):
        filename = getattr(module, '__file__', None)
        if name.startswith(('substrate_framework', 'catdft')) and filename and Path(filename).suffix == '.py':
            files.add(str(Path(filename).resolve()))
    return [bind(Path(path)) for path in sorted(files)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    target = args.output_dir.resolve()
    target.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {'started_at': datetime.now(timezone.utc).isoformat(), 'status': 'running',
              'preregistration_commit': '25286acb', 'registration': bind(HERE/'registration.json'),
              'implementation': implementation(), 'command': [sys.executable, *sys.argv],
              'environment': {'cpu_affinity': [str(x) for x in sorted(os.sched_getaffinity(0))],
                              'BLAS_threads': os.environ.get('OPENBLAS_NUM_THREADS'), 'python': sys.version},
              'scope': 'Actual material frozen ionic+fixed-wavefunction explicit contribution, NOT screened phonons, harmonic stability, calibrated material dynamics, quantum state or coherence.'}
    def deadline(*_):
        raise TimeoutError('27s work deadline inside30s external allocation')
    signal.signal(signal.SIGALRM, deadline)
    signal.setitimer(signal.ITIMER_REAL, 27)
    try:
        registration = json.loads((HERE/'registration.json').read_bytes())
        inventory_raw = (HERE/'source-inventory.json').read_bytes()
        verify_buffer(inventory_raw, registration['source_inventory'])
        inventory = json.loads(inventory_raw)
        verify_buffer((HERE/'material-input-bundle.json').read_bytes(), inventory['actual_inputs'])
        for record in inventory['files']:
            verify_buffer(Path(record['path']).read_bytes(), record)
        value, inputs = load_material_frozen(HERE/'material-input-bundle.json')
        archive = target/'loaded-ni31h6-frozen.npz'
        arrays = {'pattern_matrix': value.pattern_matrix, 'cartesian_matrix': value.cartesian_matrix,
                  'patterns_U': value.patterns, 'PH_atom_mass_amu': np.asarray(value.identity['PH_atom_mass_amu'])}
        with archive.open('xb') as stream:
            np.savez_compressed(stream, **arrays)
        with np.load(archive, allow_pickle=False) as stored:
            for name, expected in arrays.items():
                actual = stored[name]
                if actual.dtype != expected.dtype or actual.shape != expected.shape or not np.array_equal(actual, expected):
                    raise ValueError('saved material component changed: '+name)
        originals = []
        for name, raw in zip(('native-dynmat.1.0.xml', 'native-patterns.1.xml'), value.buffers):
            path = target/name
            path.write_bytes(raw)
            originals.append(bind(path, raw))
        report.update(status='complete', actual_inputs=inputs, material_component=value.record(),
                      archive=bind(archive), preserved_native_buffers=originals,
                      saved_output_readback='all four arrays retain exact dtype/shape/values; no projection/mass weighting')
    except Exception as error:
        report.update(status='stopped_requires_root_decision', failure={'type': type(error).__name__,
                      'message': str(error), 'traceback': traceback.format_exc()})
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        report['elapsed_seconds'] = time.monotonic()-started
        report['finished_at'] = datetime.now(timezone.utc).isoformat()
        with (target/'result.json').open('x') as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write('\n')
        if report['status'] == 'complete':
            (target/'.gate_done_material_frozen').write_text(report['finished_at']+'\n')
        print(json.dumps({'status': report['status'], 'elapsed_seconds': report['elapsed_seconds'],
                          'output': str(target/'result.json'), 'shape': report.get('material_component', {}).get('shape'),
                          'failure': report.get('failure')}), flush=True)
    return 0 if report['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
