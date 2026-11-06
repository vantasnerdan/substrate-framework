#!/usr/bin/env python3
"""Original-control H0/K background intake, never a material propagation/ranking.

Run one case on one CPU in 6..10 with all numerical thread variables set to 1.
The immutable registration binds existing preparation, A25 backend/oracle and
loaded archive. --array-root must be outside this Git worktree.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def binding(path):
    path = Path(path).resolve()
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(path=str(path), bytes=str(path.stat().st_size), sha256=digest)


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, choices=(
        'pristine-ni32', 'empty-vacancy-ni31', 'loaded-ni31h6', 'dispersed-ni32h6'))
    parser.add_argument('--array-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    affinity = sorted(os.sched_getaffinity(0))
    if len(affinity) != 1 or affinity[0] not in range(6, 11):
        raise ValueError('one explicitly allocated CPU in 6..10 required')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise ValueError(name + '=1 required before Python startup')
    array_root, output = args.array_root.resolve(), args.output.resolve()
    if array_root == ROOT or ROOT in array_root.parents:
        raise ValueError('--array-root must be outside Git worktree')
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    def deadline(signum, frame):
        raise TimeoutError('180-second control intake deadline; completed archives remain reusable')
    signal.signal(signal.SIGALRM, deadline)
    signal.alarm(180)
    result = dict(schema='material-control-reference-intake-v1', status='failed', case=args.case,
        cpu_affinity=affinity, scope='one selected full original frame; bounded H0/K intake only',
        exclusions=['no all-column H0 replay', 'no full mesh propagation', 'no comparative ranking',
                    'no new SCF/PH/export', 'no new equilibrium or occupation repair'])
    try:
        # Snapshot raw source bytes BEFORE imports or scientific consumption.
        paths = set((ROOT / 'src/substrate_framework').glob('*.py'))
        for attempt in ('0024', '0025', '0028'):
            paths.update((HERE.parent / attempt).glob('*.py'))
        snapshots = {str(p.resolve()): binding(p) for p in paths}
        result['registration'] = binding(HERE / 'control-reference-registration.json')
        for directory in (ROOT / 'src', HERE.parent / '0024', HERE.parent / '0025', HERE):
            sys.path.insert(0, str(directory))
        import numpy as np
        from substrate_framework.material_control_reference import (
            bound_json, retain_control_background, restore_control_background,
            selected_control_hamiltonian)
        from substrate_framework.qe_saved_native import read_bound_bytes
        from material_pilot_io import commit_block, canonical
        from run_material_reference_pilot import probe_parameters
        from run_fused_action_material import probes
        registration = bound_json(result['registration'])
        for source in registration['unchanged_scientific_sources']:
            read_bound_bytes(source)
        for proof in registration['inherited_proof']:
            read_bound_bytes(proof, max_bytes=32 << 20)
        prep_binding = registration['cases'][args.case]
        result['preparation'] = prep_binding
        receipt_binding, disposition = retain_control_background(prep_binding,
            registration['library_binding'], array_root,
            loaded_receipt_binding=(registration['loaded_background_receipt']
                                    if args.case == 'loaded-ni31h6' else None))
        background, prep, receipt = restore_control_background(prep_binding, receipt_binding)
        if prep['provenance']['case'] != args.case:
            raise ValueError('registered control case mismatch')
        result.update(background_receipt=receipt_binding, background_archive=receipt['archive'],
            disposition=disposition, original_inputs=prep['provenance'],
            background_sources=list(background.sources), selected=prep['selected'],
            selected_arrays=prep['selected_arrays'],
            loaded_content_equivalence=('raw XML/density/UPF byte digests, not path names'
                if args.case == 'loaded-ni31h6' else 'not applicable'))
        result['background_ready_seconds'] = time.monotonic() - start
        h0, state = selected_control_hamiltonian(background, prep)
        vectors, probe_meta = probes(state['coefficients'])
        h_action = h0.apply(vectors)
        (modes, _), (coefficients, _) = probe_parameters(registration['kernel_probe_inputs'])
        found = np.flatnonzero(np.all(background.fft.miller == modes[:, 0, None], axis=0))
        if found.size != 1:
            raise ValueError('original deterministic probe mode absent from this physical grid')
        fourier = np.zeros((background.fft.miller.shape[1], 2), dtype=np.complex128)
        fourier[found[0]] = 1e-4 * coefficients[:, 0]
        field = background.fft.from_g(fourier)
        k_action = background.apply_kernel(field, q_bohr_inv=(0., 0., 0.))
        if not np.isfinite(h_action).all() or not np.isfinite(k_action).all():
            raise ValueError('nonfinite actual H0/K action; no clipping or repair')
        diagnostics = dict(h0_probes=vectors, h0_action_ry=h_action,
            delta_n_updown_bohr_minus3=field, kernel_action_ry=k_action)
        diagnostic_root = array_root / args.case / output.name
        diagnostic_root.mkdir(parents=True, exist_ok=False)
        meta = dict(schema='control-reference-intake-probes-v1', preparation=prep_binding,
            selected_arrays=prep['selected_arrays'], selected=prep['selected'],
            background_receipt=receipt_binding, h0_probes=probe_meta,
            kernel_probe=dict(inputs=registration['kernel_probe_inputs'], column_zero_based=0,
                amplitude_bohr_minus3=1e-4, q_bohr_inv=[0., 0., 0.],
                qualification='same admitted deterministic complex mathematical probe, not a physical preparation'))
        diagnostic_receipt = commit_block(diagnostic_root, 'intake-diagnostics', diagnostics,
            meta, hashlib.sha256(canonical(meta)).hexdigest(), 1 << 30)
        result.update(diagnostics=diagnostic_receipt,
            h0_action_column_norms_ry=np.linalg.norm(h_action, axis=0).tolist(),
            h0_probe_metadata=probe_meta,
            kernel_action_max_abs_ry=float(np.max(np.abs(k_action))),
            original_background_imaginary_maxima=dict(background.imaginary_maxima))
        used = []
        for module in tuple(sys.modules.values()):
            filename = getattr(module, '__file__', None)
            if filename and str(Path(filename).resolve()) in snapshots:
                source = snapshots[str(Path(filename).resolve())]
                if binding(filename) != source:
                    raise ValueError('scientific consumer source changed during execution: ' + filename)
                used.append(source)
        result['implementation_sources_before_consumption'] = sorted(
            {s['path']: s for s in used}.values(), key=lambda s: s['path'])
        result['status'] = 'completed-background-restoration-and-bounded-H0-K-intake'
    except BaseException as error:
        result['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        result['elapsed_seconds'] = time.monotonic() - start
        write_json(output / 'result.json', result)
    print(json.dumps(dict(result=binding(output / 'result.json'), status=result['status']), sort_keys=True))
    signal.alarm(0)


if __name__ == '__main__':
    main()
