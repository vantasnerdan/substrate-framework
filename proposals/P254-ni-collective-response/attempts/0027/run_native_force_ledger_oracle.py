#!/usr/bin/env python3
"""Consume original Native force components under Root's explicit input binding.

No SCF/PH execution and no reconstructed physical force. Original coordinates
are independently decoded; semantic mutants carry their own valid digests.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time
import traceback
import xml.etree.ElementTree as ET

import numpy as np

import native_force_ledger_cases as cases
from substrate_framework.qe_native_force_ledger import load_native_force_components
from substrate_framework.qe_saved_native import read_bound_bytes
from substrate_framework.qe_saved_reference import load_saved_reference


def loaded_consumer_paths():
    paths = {Path(__file__).resolve(), Path(cases.__file__).resolve()}
    for name, module in tuple(sys.modules.items()):
        if name.startswith("substrate_framework") and getattr(module, "__file__", None):
            paths.add(Path(module.__file__).resolve())
    return paths


def admit_sources(bindings):
    admitted = {}
    for binding in bindings:
        _, actual = read_bound_bytes(binding, max_bytes=16 << 20)
        path = Path(actual["path"]).resolve()
        if path in admitted:
            raise ValueError("duplicate declared consumer source")
        admitted[path] = actual
    missing = loaded_consumer_paths() - admitted.keys()
    if missing:
        raise ValueError("unbound preloaded consumer source: " + ",".join(map(str, sorted(missing))))
    return admitted


def consume(inputs, require):
    reference = load_saved_reference(inputs["xml_binding"])
    species = {record.name: i for i, record in enumerate(reference.species, start=1)}
    atom_types = tuple(species[name] for name in reference.atom_names)
    raw, _ = read_bound_bytes(inputs["stdout_binding"], max_bytes=64 << 20)
    expected = cases.independent_original_rows(raw.decode("utf-8").splitlines(), atom_types)
    ledger = load_native_force_components(inputs["stdout_binding"],
                                         atom_type_indices=atom_types,
                                         evaluation_index=inputs["evaluation_index"])
    cases.compare_original(ledger, expected, require)
    cases.exercise_boundaries(require)
    return ledger.as_dict(), dict(reference.source)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--inputs-bytes", required=True)
    parser.add_argument("--inputs-sha256", required=True)
    parser.add_argument("--result", required=True, type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks, failures = [], []
    result = {"schema": "p254-0027-original-force-ledger-1", "checks": checks,
              "failures": failures, "exception": None, "ledger": None}

    def require(name, condition):
        passed = bool(condition)
        checks.append(dict(name=name, passed=passed))
        if not passed:
            failures.append(name)
            raise AssertionError(name + ": failed input-bridge prerequisite")

    try:
        # Exercise the selected lazy XML path BEFORE the source inventory.
        ET.Element("preload").find("child")
        raw, control = read_bound_bytes(dict(path=str(args.inputs), bytes=args.inputs_bytes,
                                            sha256=args.inputs_sha256), max_bytes=1 << 20)
        inputs = json.loads(raw)
        admitted = admit_sources(inputs["source_bindings"])
        result["input_binding"] = dict(control)
        result["original_source_bindings"] = [dict(admitted[path]) for path in sorted(admitted)]
        result["ledger"], result["xml_binding"] = consume(inputs, require)
        require("cold_consumer_source_closure", loaded_consumer_paths() <= admitted.keys())
    except Exception:
        result["exception"] = traceback.format_exc()
        failures.append("input_bridge_exception")
    result.update(passed=not failures, elapsed_seconds=time.monotonic() - started,
                  runtime=dict(python=sys.version, numpy=np.__version__),
                  scope="Original rounded stage-labelled component input; no SCF, XML-force reconstruction, physical preparation or exact-zero inference.",
                  independence="Original whitespace-column decoding versus fixed-width production parsing; explicit full-width, signed-zero, evaluation and malformed-record cases with separately rebound semantic mutants.")
    args.result.parent.mkdir(parents=True, exist_ok=True)
    with args.result.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(result=str(args.result), passed=result["passed"], checks=len(checks), failures=failures)))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
