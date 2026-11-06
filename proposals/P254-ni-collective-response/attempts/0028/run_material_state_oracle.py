#!/usr/bin/env python3
"""New raw-state/H0 moment oracle; no material data or old oracle replay."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))
from substrate_framework.material_reference_state import prepare_source, hamiltonian_state_moments


def check(name, actual, expected, report, *, mutation=False):
    a, b = np.asarray(actual), np.asarray(expected)
    error = float(np.max(np.abs(a-b), initial=0))
    allowance = 1e-11*(1+float(np.max(np.abs(b), initial=0)))
    passed = a.shape == b.shape and bool(np.isfinite(a).all() and np.isfinite(b).all())
    passed = passed and (error > 1e-6 if mutation else error <= allowance)
    report["checks"].append(dict(name=name, passed=passed, mutation=mutation, error=error, allowance=allowance))
    if not passed:
        raise AssertionError(name)


def calculate(report):
    # Independent exact block: all Hermitian H entries and D entries vary.
    h0,h1,h2,x,y,z,d0,d1,a,b = sp.symbols("h0 h1 h2 x y z d0 d1 a b", real=True)
    h = sp.Matrix([[h0,x+sp.I*y,z],[x-sp.I*y,h1,x],[z,x,h2]])
    psi = sp.Matrix([[1,sp.I],[sp.Rational(1,2),1],[sp.I,sp.Rational(2,3)]])
    d = sp.Matrix([[d0,a+sp.I*b],[a-sp.I*b,d1]])
    rho = psi*d*psi.H
    s,m,t = psi.H*psi, psi.H*h*psi, psi.H*h*h*psi
    drift = -sp.I*(h*rho-rho*h)
    compact = 2*sp.trace(d*s*d*t)-2*sp.trace(d*m*d*m)
    if sp.expand(compact-sp.trace(drift.H*drift)) != 0:
        raise AssertionError("exact full commutator norm identity")
    if (psi.H*drift*psi+sp.I*(m*d*s-s*d*m)).applyfunc(sp.expand) != sp.zeros(2):
        raise AssertionError("exact projected drift identity")
    report["symbolic"] = dict(full_commutator_norm=True, projected_density_drift=True)
    raw = np.asarray(psi).astype(complex)
    state = prepare_source(raw, [.2,.8], source_identity="explicit-finite-raw-frame",
                           basis_identity="explicit-three-row-space", authority="finite oracle, not a material preparation")
    hn = np.array([[.4,.3+.2j,-.1],[.3-.2j,-.6,.25],[-.1,.25,1.3]],complex)
    rn = raw @ state.density_coefficients @ raw.conj().T
    actual = hamiltonian_state_moments(raw, state, hn @ raw)
    dr = -1j*(hn@rn-rn@hn)
    check("physical_mean_number", np.trace(rn), 1., report)
    check("actual_Pauli_nonzero_spectrum", np.linalg.eigvalsh(rn), [0.,.2,.8], report)
    check("full_row_energy", actual["one_body_energy_ry"], np.trace(hn@rn), report)
    check("full_row_commutator_norm", actual["density_commutator_norm_squared_ry2"], np.sum(np.abs(dr)**2), report)
    check("full_row_projected_drift", actual["projected_density_drift_ry"], raw.conj().T@dr@raw, report)
    check("lost_cancellation_mutation", 2*np.trace(state.density_coefficients@state.gram@state.density_coefficients@actual["h0_action_gram_ry2"]), np.sum(np.abs(dr)**2), report, mutation=True)
    check("wrong_drift_sign_mutation", -actual["projected_density_drift_ry"], raw.conj().T@dr@raw, report, mutation=True)
    check("raw_probabilities_mutation", raw@np.diag([.2,.8])@raw.conj().T, rn, report, mutation=True)
    # Full-H identity is stationary, including the empty orthogonal complement.
    stationary = hamiltonian_state_moments(raw, state, raw)
    check("stationary_full_identity", stationary["density_commutator_norm_squared_ry2"], 0., report)
    check("stationary_projected_identity", stationary["projected_density_drift_ry"], np.zeros((2,2)), report)


def main(output):
    start = time.monotonic()
    report = dict(schema="P254-material-state-moment-oracle-v1", status="running", checks=[])
    code = 1
    try:
        calculate(report)
        report["status"], code = "passed", 0
    except BaseException:
        report.update(status="failed", exception=traceback.format_exc())
    report["elapsed_seconds"] = time.monotonic()-start
    report["sources"] = []
    for path in (Path(__file__), ROOT/"src/substrate_framework/material_reference_state.py"):
        raw = path.read_bytes()
        report["sources"].append(dict(path=str(path),bytes=str(len(raw)),sha256=hashlib.sha256(raw).hexdigest()))
    Path(output).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(status=report["status"], elapsed_seconds=report["elapsed_seconds"], output=output)))
    return code


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    raise SystemExit(main(parser.parse_args().output))
