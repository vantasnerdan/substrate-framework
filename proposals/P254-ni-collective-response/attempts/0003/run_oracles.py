#!/usr/bin/env python3
"""Run exact convention oracles and native raw-matrix comparisons; no simulations."""

import argparse
import json
from pathlib import Path
import platform
import time
import traceback

import numpy as np
import sympy as sp

from substrate_framework import material_force_lattice, material_modes
from substrate_framework.verification import CheckLedger
from force_lattice_oracle import force_lattice_block
from native_pilot_oracle import native_pilot_block
from units_modes_oracle import unit_block, mixed_mass_block, native_rounding_block
from native_qe_replay import replay_serialized_operators
from oracle_artifacts import sha256, save_receipt
from cache_reuse_oracle import cache_reuse_block

HERE = Path(__file__).resolve().parent
PROPOSAL = HERE.parents[1]
SOURCES = PROPOSAL / "sources/qe75-phonon-operator"




def source_hashes(ledger):
    manifest = json.loads((SOURCES/"manifest.json").read_text())
    hashes = {}
    for entry in manifest["files"] + manifest["pilot"]:
        path = PROPOSAL/entry["path"]
        digest = sha256(path)
        ledger.check(f"snapshot hash {entry['path']}", digest == entry["sha256"])
        hashes[str(path)] = digest
    for path in [*HERE.glob("*.py"),Path(material_modes.__file__),Path(material_force_lattice.__file__),SOURCES/"manifest.json"]:
        hashes[str(path)] = sha256(path)
    return hashes




def main():
    parser = argparse.ArgumentParser(description=__doc__)
    pilot = SOURCES/"pilot-attempt-02"
    parser.add_argument("--result",type=Path,default=pilot/"result.json")
    parser.add_argument("--native-directory",type=Path,default=pilot)
    parser.add_argument("--receipt",type=Path,default=HERE/"receipt.json")
    args = parser.parse_args()
    ledger = CheckLedger("P254-raw-phonon-operator")
    started = time.monotonic()
    receipt = {"schema_version":1,"status":"running","blocks":{},
               "scope":"Source-backed finite harmonic units/modes and discrete periodic Fourier conventions only.",
               "does_not_license":["converged Ni physics","infinite pair law","phonon lifetime from eta","phase/BEC","nuclear enhancement"],
               "python":platform.python_version(),"numpy":np.__version__,"sympy":sp.__version__}
    try:
        receipt["source_sha256"] = source_hashes(ledger)
        for block in (unit_block,mixed_mass_block,native_rounding_block,force_lattice_block):
            receipt["blocks"][block.__name__] = block(ledger)
            print(f"BLOCK_COMPLETE {block.__name__}",flush=True)
        receipt["native_input_sha256"] = {str(args.result):sha256(args.result)}
        data = json.loads(args.result.read_text())["data"]
        for record in data["dynamical_matrices"]:
            path = args.native_directory/Path(record["artifact"]).name
            receipt["native_input_sha256"][str(path)] = sha256(path)
        receipt["native_replay"] = replay_serialized_operators(args.result,args.native_directory)
        receipt["blocks"]["native_pilot_block"] = native_pilot_block(ledger,args.result,args.native_directory,receipt["native_replay"])
        receipt["blocks"]["cache_reuse_block"] = cache_reuse_block(ledger,args.result,args.native_directory,receipt["native_replay"])
        ledger.finish()
        receipt["status"] = "complete"
    except Exception:
        receipt["status"] = "failed"
        receipt["raw_exception"] = traceback.format_exc()
        raise
    finally:
        receipt["passed_checks"] = ledger.passed
        receipt["elapsed_seconds"] = time.monotonic()-started
        save_receipt(args.receipt,receipt)
        if receipt["status"] == "failed":
            print(f"P254_RAW_PHONON_OPERATOR_FAILED {args.receipt}",flush=True)
    print(f"P254_RAW_PHONON_OPERATOR_COMPLETE {args.receipt}",flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
