#!/usr/bin/env python3
"""Verify finite quadratic response identities; supply no Ni/LENR parameters."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

import sympy as sp

from substrate_framework.material_matter_field import (
    eliminate_bath,
    minimal_coupling_hessian,
    static_matched_bath_hessian,
)
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent


def zero_matrix(matrix):
    return all(sp.cancel(value) == 0 for value in matrix)


def minimal_coupling_block(ledger):
    m, k, w2 = sp.symbols("m k omega_field_squared", positive=True)
    c = sp.symbols("c", real=True)
    hessian = minimal_coupling_hessian(m, k, w2, c)
    minors = [sp.factor(hessian[:n, :n].det()) for n in range(1, 5)]
    ledger.check("complete-square Hessian is symmetric", hessian == hessian.T)
    ledger.check("all leading principal minors positive", all(
        minor.is_positive is True for minor in minors))
    without_self_term = hessian.copy().as_mutable()
    without_self_term[2, 2] -= c**2 / m
    mutant = without_self_term.subs({m: 1, k: 1, w2: 1, c: 2})
    ledger.check("dropping self-term exposes a negative energy direction",
                 mutant[:3, :3].det() < 0)
    canonical = sp.diag(sp.Matrix([[0, 1], [-1, 0]]),
                        sp.Matrix([[0, 1], [-1, 0]]))
    drift = canonical * hessian
    ledger.check("Hamiltonian energy derivative vanishes",
                 zero_matrix(drift.T * hessian + hessian * drift))
    wrong_drift = drift.as_mutable()
    wrong_drift[0, 0] += 1
    ledger.check("added gain violates energy identity", not zero_matrix(
        wrong_drift.T * hessian + hessian * wrong_drift))
    return {"hessian": str(hessian), "leading_principal_minors": list(map(str, minors)),
            "scope": "Finite positive quadratic energy only; not a Ni charge or a phase theorem."}


def static_screening_block(ledger):
    k11, k12, k22, b1, b2 = sp.symbols("k11 k12 k22 b1 b2", positive=True)
    c11, c12, c21, c22, z = sp.symbols("c11 c12 c21 c22 z")
    retained = sp.Matrix([[k11, k12], [k12, k22]])
    bath = sp.diag(b1, b2)
    coupling = sp.Matrix([[c11, c12], [c21, c22]])
    full = static_matched_bath_hessian(retained, bath, coupling)
    zero = eliminate_bath(full, 2, 0)
    ledger.check("static elimination returns calibrated stiffness",
                 zero_matrix(zero - retained))
    expected = retained - z**2 * sp.eye(2)
    expected += coupling * bath.inv() * coupling.T
    expected -= coupling * (bath - z**2 * sp.eye(2)).inv() * coupling.T
    ledger.check("dynamic Schur complement retains static subtraction",
                 zero_matrix(eliminate_bath(full, 2, z) - expected))
    missing_counterterm = full.copy().as_mutable()
    missing_counterterm[:2, :2] = retained
    ledger.check("missing counterterm fails static calibration", not zero_matrix(
        eliminate_bath(missing_counterterm, 2, 0) - retained))
    transposed_coupling = static_matched_bath_hessian(retained, bath, coupling.T)
    ledger.check("asymmetric coupling exposes orientation error", not zero_matrix(
        eliminate_bath(transposed_coupling, 2, z) - expected))
    return {"effective_inverse_response": str(expected),
            "scope": "Finite nonsingular mass-weighted bath; physical bath spectrum must be calculated."}


def direct_resolvent_block(ledger):
    retained = sp.Matrix([[3, 1], [1, 5]])
    bath = sp.Matrix([[7, 2], [2, 11]])
    coupling = sp.Matrix([[1, 2], [3, -1]])
    full = static_matched_bath_hessian(retained, bath, coupling)
    z = (2 + sp.I) / 3
    direct = (full - z**2 * sp.eye(4)).inv()[:2, :2]
    reduced = eliminate_bath(full, 2, z).inv()
    ledger.check("direct full inverse agrees with reduced response",
                 zero_matrix(direct - reduced))
    imaginary = direct.applyfunc(sp.im)
    ledger.check("retarded positive-frequency imaginary response is positive",
                 imaginary[0, 0] > 0 and imaginary.det() > 0)
    advanced = (full - sp.conjugate(z)**2 * sp.eye(4)).inv()[:2, :2]
    ledger.check("advanced-sign mutation reverses absorptive sign",
                 advanced.applyfunc(sp.im)[0, 0] < 0)
    return {"fixture_scope": "Exact asymmetric rational example, not material parameters.",
            "regularizer_scope": "Upper-half-plane frequency is a resolvent prescription, not a physical linewidth.",
            "direct_response": str(direct)}


def save_receipt(payload):
    target = HERE / "receipt.json"
    temporary = target.with_suffix(".json.tmp")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, target)
    descriptor = os.open(HERE, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def main():
    started = time.monotonic()
    ledger = CheckLedger("P254-quadratic-blocks")
    blocks = {}
    for construct in (minimal_coupling_block, static_screening_block,
                      direct_resolvent_block):
        blocks[construct.__name__] = construct(ledger)
        print(f"BLOCK_COMPLETE {construct.__name__}", flush=True)
    ledger.finish()
    module = Path(__import__("substrate_framework.material_matter_field",
                             fromlist=["__file__"]).__file__)
    save_receipt({"schema_version": 1, "status": "complete",
                  "evidence_scope": "symbolic finite quadratic identities and exposing mutations only",
                  "does_not_license": ["Ni material parameters", "spontaneous coherence", "nuclear enhancement"],
                  "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in (Path(__file__), module)},
                  "sympy_version": sp.__version__, "elapsed_seconds": time.monotonic() - started,
                  "checks": ledger.passed, "blocks": blocks})
    print("P254_BLOCKS_COMPLETE receipt.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
