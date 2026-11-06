#!/usr/bin/env python3
"""A27 NEW finite ionic/composition oracle; execution is Root-only.

No old oracle, material input decoder or Native executable is called. The
scalar comparator uses literal ordered pair energies and five-point energy
differentiation; the composition comparator differentiates its own polynomial.
Maximum verdict: small finite numerical/symbolic evidence, never material
PBE equilibrium or convergence of arbitrary user-selected Ewald cutoffs.
Run only under an explicit source-bound Root budget: --result NEW_PATH.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
CHECKS = []


def compare(name, actual, expected, *, atol=2e-11, rtol=2e-10, mutation=False):
    import numpy as np
    a,b = np.asarray(actual),np.asarray(expected)
    compatible = a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    residual = float(np.max(np.abs(a-b),initial=0)) if compatible else None
    agrees = compatible and bool(np.allclose(a,b,atol=atol,rtol=rtol))
    passed = compatible and (not agrees if mutation else agrees)
    CHECKS.append(dict(name=('mutation_' if mutation else '')+name,passed=bool(passed),
                       mutation=mutation,max_abs_residual=residual,atol=atol,rtol=rtol))
    if not passed:
        raise AssertionError(name+': '+str(residual))


def mutation(name,wrong,expected,**kwargs):
    compare(name,wrong,expected,mutation=True,**kwargs)


def exact(name,difference,*,mutation=False):
    import sympy as sp
    entries = list(difference) if isinstance(difference,sp.MatrixBase) else [difference]
    residuals = [sp.simplify(entry) for entry in entries]
    zero = all(entry == 0 for entry in residuals)
    passed = not zero if mutation else zero
    CHECKS.append(dict(name=('mutation_' if mutation else '')+name,passed=bool(passed),
                       mutation=mutation,symbolic_residuals=[str(entry) for entry in residuals]))
    if not passed:
        raise AssertionError(name+': '+str(residuals))


def exact_mutation(name,difference):
    exact(name,difference,mutation=True)


def source_bindings():
    paths = [Path(__file__), Path(__file__).with_name('ionic_scalar_oracle.py'),
             Path(__file__).with_name('ionic_composition_oracle.py')]
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'qe_ionic_ewald.py','qe_ewald_sums.py','prepared_reference.py',
        'prepared_factor_dynamics.py','prepared_covariance_source.py',
        'prepared_time_observables.py','physical_factor_preparation.py')]
    result = []
    for path in paths:
        data = path.read_bytes()
        result.append(dict(path=str(path),bytes=str(len(data)),sha256=hashlib.sha256(data).hexdigest()))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result',required=True,type=Path)
    args = parser.parse_args()
    # Refuse to overwrite an earlier execution; original failures are evidence.
    with args.result.open('x') as output:
        report = dict(schema='p254-a27-ionic-reference-oracle-1',sources=source_bindings(),
                      scope='NEW synthetic finite scalar/composition only; no material or time trajectory',
                      blocks=['scalar constants','energy derivative and independent force',
                              'symmetry units truncation sensitivities','symbolic common scalar',
                              'nonstationary mean and affine energy/probe consumers'],checks=CHECKS)
        try:
            from ionic_scalar_oracle import check_scalar
            from ionic_composition_oracle import check_polynomial,check_composition
            check_scalar(compare,mutation)
            check_polynomial(exact,exact_mutation)
            check_composition(compare,mutation)
            report['status'] = 'passed_finite_oracles_only'
            status = 0
        except Exception:
            report['status'] = 'failed'
            report['traceback'] = traceback.format_exc()
            status = 1
        report['executed_checks'] = len(CHECKS)
        output.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return status


if __name__ == '__main__':
    raise SystemExit(main())
