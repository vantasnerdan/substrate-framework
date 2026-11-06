"""Admitted actual37atom Gamma ionic blocks; Root owns all-exit notification.

Usage: python run_material_ionic_pilot.py --output DIR --admission JSON
Add --context-only for actual manager context/source preflight, no material.
See material_pilot_io.validate_admission for the complete grant schema.
Saved-band overlaps are explicit bilinears, not a complete-space propagator;
raw Marzari-Vanderbilt occupations are metadata, never quantum probabilities.
"""
import time
BOOT = time.monotonic()
import argparse
import json
from pathlib import Path
import signal
import sys
import traceback
import uuid


def main():
    target = None
    summary = {'status': 'failed', 'blocks': [], 'first_useful_block_seconds': None,
               'scope': 'bare ionic Gamma only; no screening/H0/physical ensemble',
               'launch_id': uuid.uuid4().hex}
    code = 1
    try:
        import numpy as np
        from material_pilot_inputs import load_material, coordinate_arrays, DIRECTIONS
        from material_pilot_io import (runtime_context, implementation, validate_admission,
                                      commit_block, recover_block, save_json)
        from substrate_framework.qe_material_ionic_action import saved_band_overlap
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--output', required=True, type=Path)
        parser.add_argument('--admission', required=True, type=Path)
        parser.add_argument('--context-only', action='store_true')
        args = parser.parse_args()
        target = args.output.resolve()
        target.mkdir(parents=True, exist_ok=True)
        summary['context'] = runtime_context()
        summary['implementation_sources'] = implementation()
        grant, inventory, recipe, identity, admission = validate_admission(
            args.admission, summary['implementation_sources'], summary['context'], args.context_only)
        summary.update(recipe=recipe, recipe_sha256=identity, admission=admission,
                       remaining_seconds_granted=grant['remaining_seconds'])
        if args.context_only:
            summary['status'] = 'context-only'
            code = 0
        else:
            def interrupted(signum, frame):
                raise InterruptedError('signal '+str(signum)+'; only completed blocks recoverable')
            for signum in (signal.SIGTERM, signal.SIGINT):
                signal.signal(signum, interrupted)
            def budget():
                elapsed = time.monotonic()-BOOT
                if elapsed >= grant['remaining_seconds']:
                    raise TimeoutError('explicit remaining cumulative allocation exhausted')
                if summary['first_useful_block_seconds'] is None and elapsed >= grant['first_block_seconds']:
                    raise TimeoutError('no useful full-PW block within explicit first-block allowance')
            def block(key, arrays, metadata):
                budget()
                previous = recover_block(target, key, identity, metadata) if grant['resume'] else None
                result = previous if previous else commit_block(
                    target, key, arrays, metadata, identity, grant['output_max_bytes'])
                budget()
                return result
            if not grant['resume'] and any(target.glob('*.npz')):
                raise ValueError('fresh admission cannot consume existing blocks')
            budget()
            ref, waves, density, action, directions = load_material(inventory)
            if implementation() != summary['implementation_sources']:
                raise ValueError('source consumers changed after frozen preflight')
            arrays = coordinate_arrays(ref, waves, density)
            arrays['directions_dimensionless'] = np.stack(directions)
            coordinate_metadata = {'meaning': 'original full37atom and physical PW identity',
                                   'occupations': 'raw MV metadata, NOT probabilities',
                                   'energies': 'Ry', 'positions': 'bohr', 'reciprocal': 'bohr^-1',
                                   'direction_order': DIRECTIONS}
            summary['coordinates'] = block('coordinates', arrays, coordinate_metadata)
            prepared = [action.prepare(u) for u in directions]
            for spin, wave in enumerate(waves, 1):
                # One full saved-wavefunction conjugation per spin, not per block.
                adjoint = wave.coefficients.conj().T.copy()
                for spec, operation in zip(DIRECTIONS, prepared):
                    matrix = np.empty((346, 346), dtype=np.complex128)
                    prefix = 'spin'+str(spin)+'-atom'+str(spec['atom_one_based'])+'-x'
                    for start in range(0, 346, grant['block_bands']):
                        budget()
                        stop = min(start+grant['block_bands'], 346)
                        key = prefix+'-bands'+str(start+1).zfill(3)+'-'+str(stop).zfill(3)
                        metadata = {'spin_one_based': spin, 'k_spin_local_one_based': 1,
                                    'direction': spec, 'source_bands_one_based': list(range(start+1, stop+1)),
                                    'target_saved_bands_one_based': list(range(1, 347)),
                                    'full_pw_rows': 'coordinates.npz:pw_miller original column order',
                                    'atom_order': 'coordinates.npz:atom_one_based/atom_symbols',
                                    'units': 'Ry/bohr', 'projection': 'Psi_saved_dagger @ gPsi; NOT closure'}
                        receipt = recover_block(target, key, identity, metadata) if grant['resume'] else None
                        if receipt:
                            with np.load(target/(key+'.npz'), allow_pickle=False) as loaded:
                                matrix[:, start:stop] = loaded['saved_band_overlap_ry_per_bohr']
                        else:
                            full = operation.apply(wave.coefficients[:, start:stop])
                            overlap = saved_band_overlap(adjoint, full)
                            receipt = commit_block(target, key,
                                {'gpsi_ry_per_bohr': full, 'saved_band_overlap_ry_per_bohr': overlap,
                                 'source_bands_one_based': np.arange(start+1, stop+1, dtype=np.int64)},
                                metadata, identity, grant['output_max_bytes'])
                            matrix[:, start:stop] = overlap
                        summary['blocks'].append(receipt)
                        # A late first commit is still preserved but fails admission.
                        budget()
                        if summary['first_useful_block_seconds'] is None:
                            summary['first_useful_block_seconds'] = time.monotonic()-BOOT
                        summary['elapsed_seconds'] = time.monotonic()-BOOT
                        save_json(target/'summary.json', summary, rotate=True)
                    final_metadata = {'spin_one_based': spin, 'direction': spec,
                                      'rows': 'saved target bands1..346', 'columns': 'saved source bands1..346',
                                      'units': 'Ry/bohr', 'meaning': 'explicit saved-band bilinear; not full H0'}
                    summary.setdefault('saved_band_matrices', []).append(block(prefix+'-saved-band-matrix',
                        {'saved_band_overlap_ry_per_bohr': matrix}, final_metadata))
                del adjoint
            summary['status'] = 'complete'
            code = 0
    except BaseException as error:
        summary['failure'] = {'type': type(error).__name__, 'message': str(error),
                              'traceback': traceback.format_exc()}
        if isinstance(error, SystemExit):
            code = int(error.code or 0)
    finally:
        summary['elapsed_seconds'] = time.monotonic()-BOOT
        summary['exit_code'] = code
        if target is not None:
            try:
                # Unique terminal records preserve ALL failures across recovery;
                # summary/.prev is the current/previous view, never failure erasure.
                save_json(target/('terminal-'+summary['launch_id']+'.json'), summary)
                save_json(target/'summary.json', summary, rotate=True)
            except BaseException as error:
                summary['terminal_write_failure'] = repr(error)
                code = 1
                summary['exit_code'] = code
        print(json.dumps(summary, allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
