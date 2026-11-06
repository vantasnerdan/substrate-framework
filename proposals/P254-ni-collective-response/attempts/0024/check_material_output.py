"""Consume completed material archives; no ionic action is recomputed.

The Gamma adjoint allowance is carried from run_material_ionic_oracle.py.
Recovery covers committed archive/receipt pairs, not interrupted operations.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from material_pilot_io import binding, implementation, recover_block, save_json
from substrate_framework.qe_saved_native import read_bound_bytes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--terminal', required=True, type=Path)
    parser.add_argument('--bytes', required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    report = {'passed': False, 'implementation_sources': implementation(), 'matrices': []}
    try:
        raw, source = read_bound_bytes({'path': str(args.terminal), 'bytes': args.bytes,
                                       'sha256': args.sha256})
        original = json.loads(raw)
        report['terminal'] = dict(source)
        target = args.terminal.parent
        receipts = [original['coordinates'], *original['blocks'], *original['saved_band_matrices']]
        for receipt in receipts:
            key = Path(receipt['archive']['path']).stem
            recovered = recover_block(target, key, original['recipe_sha256'], receipt['metadata'])
            if recovered != receipt:
                raise ValueError('completed archive receipt changed: '+key)
        report['exact_completed_receipts'] = len(receipts)
        for receipt in original['saved_band_matrices']:
            with np.load(receipt['archive']['path'], allow_pickle=False) as archive:
                matrix = archive['saved_band_overlap_ry_per_bohr']
            error = float(np.max(np.abs(matrix-matrix.conj().T)))
            allowance = 2e-11+2e-10*float(np.max(np.abs(matrix)))
            if not np.isfinite(matrix).all() or error > allowance:
                raise ValueError('raw complete saved-band matrix violates Gamma adjoint')
            report['matrices'].append({'archive': receipt['archive'], 'metadata': receipt['metadata'],
                                      'adjoint_max_error': error, 'allowance': allowance})
        first = original['blocks'][0]
        key = Path(first['archive']['path']).stem
        wrong_metadata = dict(first['metadata'], spin_one_based=2)
        refusals = []
        for name, identity, metadata in [('changed-recipe', '0'*64, first['metadata']),
                                        ('wrong-spin', original['recipe_sha256'], wrong_metadata)]:
            try:
                recover_block(target, key, identity, metadata)
            except ValueError as error:
                refusals.append({'name': name, 'observed_error': str(error)})
            else:
                raise ValueError('unsafe completed archive accepted: '+name)
        report['identity_refusals'] = refusals
        report['passed'] = True
    except Exception as error:
        report['failure'] = {'type': type(error).__name__, 'message': str(error)}
    report['elapsed_seconds'] = time.monotonic()-started
    save_json(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k != 'implementation_sources'}), flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
