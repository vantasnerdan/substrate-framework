#!/usr/bin/env python3
"""Exercise only the warned local-G0 block; retain all other settled0024 checks."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import warnings

from scipy.integrate import IntegrationWarning
from substrate_framework.qe_saved_native import read_bound_bytes
from substrate_framework.qe_nc_upf import read_nc_upf
from substrate_framework.qe_nc_radial import build_radial_tables
from run_offline_nc_oracle import local_g0_quadrature_check

HERE = Path(__file__).resolve().parent


def binding(path):
    raw = path.read_bytes()
    return {'path':str(path),'bytes':str(len(raw)),
            'sha256':hashlib.sha256(raw).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    report = {'checks':[], 'mutations':[], 'quadrature_error_estimates':[],
              'scope':'Only local-G0 spline-knot quadrature repair; no Native/comparator replay.',
              'implementation':[binding(Path(__file__)), binding(HERE/'run_offline_nc_oracle.py')]}
    try:
        admission = json.loads((HERE/'source-wave-admission-001.json').read_bytes())
        registration = json.loads(read_bound_bytes(admission['registration'])[0])
        inventory = json.loads(read_bound_bytes(registration['access'])[0])
        background = next(item for item in inventory['existing_comparator_all_files']
                          if item['path'].endswith('/background.json'))
        omega = json.loads(read_bound_bytes(background)[0])['omega']
        report['inputs'] = [background]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', IntegrationWarning)
            for item in inventory['material_inputs']:
                if not item['path'].endswith('.upf'):
                    continue
                upf = read_nc_upf(item['path'], expected_sha256=item['sha256'],
                                  expected_bytes=item['bytes'])
                report['inputs'].append(item)
                table = build_radial_tables(upf, omega, .37)
                local_g0_quadrature_check(report, table)
            report['integration_warnings'] = [str(warning.message) for warning in caught]
        report['passed'] = (not report['integration_warnings']
                            and all(check['passed'] for check in report['checks'])
                            and all(item['exposed'] for item in report['mutations']))
    except Exception as error:
        report['passed'] = False
        report['exception'] = {'type':type(error).__name__, 'message':str(error)}
    report['elapsed_seconds'] = time.monotonic()-started
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2, allow_nan=False), flush=True)
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
