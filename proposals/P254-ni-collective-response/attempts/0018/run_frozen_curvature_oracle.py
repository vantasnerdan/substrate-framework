#!/usr/bin/env python3
"""Import the two pinned generation-003 frozen curvatures and prove convention sensitivity.

Only --output is configurable: actual inputs are fixed by attempts/0018's source
inventory. Open the output exclusively; never overwrite an earlier receipt.
Engineering import/proof only: no phonon, coherence, loss or nuclear conclusion.
"""
from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import traceback
from typing import Any

import numpy as np
import sympy as sp
import catdft.qe_phonon_io as xml

import substrate_framework.qe_frozen_curvature as frozen
from substrate_framework.qe_response_bundle import load_bundle

ATTEMPT = Path(__file__).resolve().parent
BUNDLE_PINS = {
    "gamma": ("221099", "d2725814eabbb06ad9e662b4dbff28fe339f7f3612b7589f124700f485d99278"),
    "half-b1": ("221099", "c0ef5d5c4ec78e2a7f7a4fb8e82a129b205a95365591a339c13ce5a0eff5818b"),
}
AUTHORITY_NAMES = {
    "dynmat0.f90", "dynmat_us.f90", "d2ionq.f90", "rotate_pattern_add.f90",
    "ph_restart.f90", "xmltools.f90", "dyndia.f90",
}
OPTIONAL_FLAGS = ("okvan", "nlcc_any", "noncolin", "do_comp_mt", "do_cutoff_2D",
                  "llondon", "lxdm", "ldftd3", "lda_plus_u")


def receipt(path: Path, data: bytes) -> dict[str, str]:
    return {"path": str(path), "bytes": str(len(data)),
            "sha256": hashlib.sha256(data).hexdigest()}


def bound_read(record: dict[str, Any]) -> tuple[bytes, dict[str, str]]:
    data = Path(record["path"]).read_bytes()
    return data, frozen.verify_buffer(data, record)


def source_receipts() -> tuple[dict[str, Any], dict[str, Any]]:
    inventory_path = ATTEMPT / "source-inventory.json"
    raw = inventory_path.read_bytes()
    inventory = json.loads(raw)
    sources: dict[str, Any] = {
        "preregistration_commit": "b3f358db", "source_inventory": receipt(inventory_path, raw),
        "registration": receipt(ATTEMPT / "registration.json", (ATTEMPT / "registration.json").read_bytes()),
        "authority_inventories": [], "actual_native_authorities": [],
    }
    # Verify settled report and inventories before interpreting their contents.
    for record in inventory["settled_prior_inputs"]:
        name = Path(record["path"]).name
        if name not in ("native-action-assessment.capture.json", "native-action-source-inventory.json",
                        "native-action-additional-source-inventory.json"):
            continue
        data, binding = bound_read(record)
        sources["authority_inventories"].append(binding)
        if name == "native-action-assessment.capture.json":
            continue
        for authority in json.loads(data)["files"]:
            if Path(authority["path"]).name in AUTHORITY_NAMES:
                _, authority_binding = bound_read(authority)
                sources["actual_native_authorities"].append(authority_binding)
    # The complete, hash-pinned NativeAction report table additionally binds
    # ph_restart; its table receipt is not repeated in the two shorter inventories.
    restart_pin = {
        "path": "/opt/qebuild/q-e-qe-7.5/PHonon/PH/ph_restart.f90",
        "bytes": "50465",
        "sha256": "d39efbc92b7151686e7b3dfd8656f7805adbeb52bafaad93cdfae8916ba58e25",
    }
    _, restart_binding = bound_read(restart_pin)
    sources["actual_native_authorities"].append(restart_binding)
    names = {Path(record["path"]).name for record in sources["actual_native_authorities"]}
    required = AUTHORITY_NAMES
    if not required <= names:
        raise ValueError(f"missing preregistered frozen-curvature authority pins: {sorted(required - names)}")
    sources["implementation"] = [receipt(path, path.read_bytes()) for path in
        (Path(frozen.__file__).resolve(), Path(xml.__file__).resolve(), Path(__file__).resolve())]
    return inventory, sources


def flags(metadata: Any) -> dict[str, Any]:
    result = {}
    for name in OPTIONAL_FLAGS:
        if name in metadata:
            value = metadata[name]
            if type(value) is not bool:
                raise ValueError(f"actual Native applicability flag is not Boolean: {name}")
            result[name] = {"status": "available", "value": value,
                            "source": "closed bundle background.json"}
        else:
            result[name] = {"status": "unavailable", "expected_type": "native Boolean",
                            "reason": "not serialized in the existing generation-003 background metadata; no value inferred"}
    return result


def source_identity(bundle: Any, manifest: dict[str, str], label: str) -> dict[str, Any]:
    bg, meta = bundle.background, bundle.background.metadata
    if meta["nat"] != 2 or meta["iq"] != 1:
        raise ValueError("registered two-atom / iq=1 identity changed")
    q = bg.array("q_cart_2pi_alat", (3,), "float64")
    expected = np.zeros(3) if label == "gamma" else np.array(
        [-0.35355339059327379, 0.35355339059327379, 0.35355339059327379])
    if not np.allclose(q, expected, rtol=0, atol=1e-14):
        raise ValueError("actual q differs from the settled producer representative")
    at, reciprocal = bg.array("at", (3, 3), "float64"), bg.array("bg", (3, 3), "float64")
    return {
        "label": label, "nat": meta["nat"], "iq": meta["iq"], "bundle_manifest": manifest,
        "q_cart_2pi_alat": q.tolist(), "q_units": meta["q_units"],
        "q_crystal_representative": (at.T @ q).tolist(), "q_direction": meta["q_direction"],
        "q_identity_binding": "same producer directory; patterns QPOINT_NUMBER agrees with closed bundle iq; dynmat.1.0.xml is irrep zero",
        "cell": {"alat_bohr": meta["alat"], "omega_bohr3": meta["omega"],
                 "at_columns_alat": at.tolist(), "bg_columns_2pi_alat": reciprocal.tolist(),
                 "atom_positions": {"status": "unavailable", "expected_type": "Cartesian atomic coordinates",
                                    "reason": "not serialized in this bundle background; manifest retains actual producer identity"}},
        "functional": meta["functional"], "applicability_flags": flags(meta),
    }


def import_actual(label: str, records: list[dict[str, Any]]) -> dict[str, Any]:
    dyn = [record for record in records if f"/{label}/" in record["path"] and
           Path(record["path"]).name == "dynmat.1.0.xml"]
    patterns = [record for record in records if f"/{label}/" in record["path"] and
                Path(record["path"]).name == "patterns.1.xml"]
    if len(dyn) != 1 or len(patterns) != 1:
        raise ValueError("inventory must bind exactly one original dynmat and patterns pair per q")
    producer = Path(dyn[0]["path"]).parents[3]
    if Path(patterns[0]["path"]).parent != Path(dyn[0]["path"]).parent:
        raise ValueError("curvature and patterns are not the same producer's phsave")
    bundle = load_bundle(producer / "bundle" / "q-1")
    size, sha256 = BUNDLE_PINS[label]
    manifest = frozen.verify_buffer(bundle.buffers["manifest.json"],
        {"path": str(bundle.root / "manifest.json"), "bytes": size, "sha256": sha256})
    dyn_bytes, _ = bound_read(dyn[0])
    pattern_bytes, _ = bound_read(patterns[0])
    imported = frozen.parse_frozen_curvature(dyn_bytes, pattern_bytes,
        dynmat_source=dyn[0], patterns_source=patterns[0],
        identity=source_identity(bundle, manifest, label))
    # Compare to existing typed complex native pattern columns, not an identity
    # replacement. Every bare record must agree with the XML mode/irrep ordering.
    max_pattern_error = 0.0
    modes = set()
    for bare in bundle.bare:
        mode = bare.metadata["mode"]
        if (bare.metadata["irrep"], bare.metadata["perturbation_index"]) != imported.pattern_order[mode - 1]:
            raise ValueError("XML and closed-bundle irrep/perturbation ordering differ")
        difference = np.max(np.abs(bare.array("Cartesian_complex_pattern", (6,), "complex128") -
                                   imported.patterns[:, mode - 1]))
        max_pattern_error = max(max_pattern_error, float(difference))
        modes.add(mode)
    if modes != set(range(1, 7)) or max_pattern_error > 1e-12:
        raise ValueError("existing bundle patterns do not corroborate complete XML pattern basis")
    result = imported.record()
    result["original_xml_base64"] = {
        "dynmat.1.0.xml": base64.b64encode(imported.buffers[0]).decode("ascii"),
        "patterns.1.xml": base64.b64encode(imported.buffers[1]).decode("ascii"),
    }
    result["bundle_pattern_comparison"] = {"records": len(bundle.bare), "modes": sorted(modes),
        "maxabs_difference": max_pattern_error, "absolute_allowance": 1e-12,
        "scope": "decimal-XML basis corroboration, not a Native in-memory bit-exact claim"}
    result["import_pass"] = True
    return result


def fixture_buffers(c: sp.Matrix, u: sp.Matrix) -> tuple[bytes, bytes]:
    def pair(value: Any) -> str:
        return f"{sp.N(sp.re(value), 17)} {sp.N(sp.im(value), 17)}"
    # Explicit column-then-row enumeration is independent of the production
    # reader's reshape implementation, and the fixture is not transpose invariant.
    serialized = "\n".join(pair(c[row, col]) for col in range(6) for row in range(6))
    dyn = f"<Root><PM_HEADER><DONE_IRR>true</DONE_IRR></PM_HEADER><PARTIAL_MATRIX><PARTIAL_DYN>{serialized}</PARTIAL_DYN></PARTIAL_MATRIX></Root>"
    columns = []
    for col in range(6):
        values = "\n".join(pair(u[row, col]) for row in range(6))
        columns.append(f"<PERTURBATION.{col + 1}><DISPLACEMENT_PATTERN>{values}</DISPLACEMENT_PATTERN></PERTURBATION.{col + 1}>")
    patterns = ("<Root><IRREPS_INFO><QPOINT_NUMBER>1</QPOINT_NUMBER><QPOINT_GROUP_RANK>1</QPOINT_GROUP_RANK>"
                "<MINUS_Q_SYM>false</MINUS_Q_SYM><NUMBER_IRR_REP>1</NUMBER_IRR_REP>"
                "<REPRESENTION.1><NUMBER_OF_PERTURBATIONS>6</NUMBER_OF_PERTURBATIONS>" +
                "".join(columns) + "</REPRESENTION.1></IRREPS_INFO></Root>")
    return dyn.encode(), patterns.encode()


def exact_matrix(array: np.ndarray) -> sp.Matrix:
    return sp.Matrix(array.shape[0], array.shape[1], lambda i, j:
        sp.Rational(str(float(array[i, j].real))) + sp.I * sp.Rational(str(float(array[i, j].imag))))


def convention_oracle() -> dict[str, Any]:
    """Exact independent quadratic-energy oracle, plus exposed deliberate defects."""
    c = sp.diag(2, 5, 7, 11, 13, 17)
    for row, col, value in ((0, 1, 1 + 2 * sp.I), (0, 2, 2 - sp.I), (3, 5, -1 + 3 * sp.I)):
        c[row, col], c[col, row] = value, sp.conjugate(value)
    u = sp.eye(6)
    u[0, 0], u[0, 1] = sp.Rational(3, 5), 4 * sp.I / 5
    u[1, 0], u[1, 1] = sp.Rational(4, 5), -3 * sp.I / 5
    x = sp.Matrix([1 + sp.I, 2 - sp.I, 1 - 2 * sp.I, -1 + sp.I, 2, 3 + sp.I])
    if c != c.H or u.H * u != sp.eye(6):
        raise AssertionError("invalid exact Hermitian/unitary sensitivity fixture")
    # Independent explicit amplitudes and scalar quadratic form: neither calls
    # production transformation nor obtains an expected value by roundtripping.
    p = [(3 * x[0] + 4 * x[1]) / 5, (-4 * sp.I * x[0] + 3 * sp.I * x[1]) / 5, *x[2:]]
    diagonal_energy = sum(weight * sp.conjugate(amplitude) * amplitude
                          for weight, amplitude in zip((2, 5, 7, 11, 13, 17), p))
    cross = ((1 + 2 * sp.I) * sp.conjugate(p[0]) * p[1] +
             (2 - sp.I) * sp.conjugate(p[0]) * p[2] +
             (-1 + 3 * sp.I) * sp.conjugate(p[3]) * p[5])
    expected_energy = sp.simplify((diagonal_energy + 2 * sp.re(cross)) / 2)
    dyn, patterns = fixture_buffers(c, u)
    identity = {"nat": 2, "iq": 1, "q_cart_2pi_alat": [0, 0, 0],
                "cell": {"scope": "exact algebra fixture, not a material input"},
                "bundle_manifest": {"scope": "synthetic proof fixture only"}, "applicability_flags": {}}
    imported = frozen.parse_frozen_curvature(dyn, patterns,
        dynmat_source=receipt(Path("exact-fixture/dynmat.1.0.xml"), dyn),
        patterns_source=receipt(Path("exact-fixture/patterns.1.xml"), patterns), identity=identity)
    parsed_c, parsed_u = exact_matrix(imported.pattern_matrix), exact_matrix(imported.patterns)
    if parsed_c != c or parsed_u != u:
        raise AssertionError("production byte parser changed exact fixture entries")
    exact_cart = parsed_u * parsed_c * parsed_u.H
    exact_energy = sp.simplify((x.H * exact_cart * x)[0] / 2)
    if exact_energy != expected_energy:
        raise AssertionError("complex pattern quadratic energy does not equal independent scalar energy")
    numeric_error = float(np.max(np.abs(imported.cartesian_matrix - np.array(exact_cart, dtype=complex))))
    if numeric_error > 1e-12:
        raise AssertionError("production Cartesian transform disagrees with exact energy matrix")
    wrong_matrices = {
        "C_order_instead_of_F": parsed_u * parsed_c.T * parsed_u.H,
        "transpose_instead_of_adjoint": parsed_u * parsed_c * parsed_u.T,
        "reversed_U_orientation": parsed_u.H * parsed_c * parsed_u,
        "identity_pattern_assumption": parsed_c,
        "swapped_real_imag_pairs": parsed_u * (parsed_c.applyfunc(sp.im) + sp.I * parsed_c.applyfunc(sp.re)) * parsed_u.H,
    }
    mutations = {}
    for name, matrix in wrong_matrices.items():
        energy = sp.simplify((x.H * matrix * x)[0] / 2)
        delta = sp.simplify(energy - expected_energy)
        mutations[name] = {"energy": str(energy), "difference_from_independent_energy": str(delta),
                           "exposed": delta != 0}
    if not all(result["exposed"] for result in mutations.values()):
        raise AssertionError("a required convention defect was not exposed by the energy oracle")
    # Byte binding and DONE_IRR are exercised against the same pure entry point.
    rejections = {}
    for name, bad_dyn, binding in (
        ("changed_bound_byte", dyn + b" ", receipt(Path("fixture"), dyn)),
        ("DONE_IRR_false", dyn.replace(b">true<", b">false<"),
         receipt(Path("fixture"), dyn.replace(b">true<", b">false<"))),
    ):
        try:
            frozen.parse_frozen_curvature(bad_dyn, patterns, dynmat_source=binding,
                patterns_source=receipt(Path("fixture-patterns"), patterns), identity=identity)
        except ValueError as error:
            rejections[name] = {"rejected": True, "error": str(error)}
        else:
            raise AssertionError(f"required invalid-input rejection missing: {name}")
    return {"proof_pass": True, "fixture_C_pattern": str(c), "fixture_U": str(u),
            "fixture_displacement": str(x), "energy_definition": "one-half xdagger C_cart x",
            "independent_scalar_energy": str(expected_energy), "parsed_exact_energy": str(exact_energy),
            "production_cartesian_maxabs_error": numeric_error, "numerical_allowance": 1e-12,
            "mutations": mutations, "invalid_input_rejections": rejections,
            "scope": "exact complex-Hermitian/unitary convention sensitivity, not a material energy calculation"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new JSON receipt, no overwrites")
    args = parser.parse_args()
    if args.output.suffix.lower() != ".json":
        parser.error("--output must name a new .json file")
    # Reserve before consuming inputs: an existing receipt cannot be overwritten,
    # including a competing process racing to use the same output name.
    with args.output.open("x", encoding="utf-8") as output:
        report: dict[str, Any] = {"schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
            "command": [sys.executable, *sys.argv], "actual_inputs": {},
            "scope": "engineering-only absolute frozen-curvature import and convention proof; no phonon/coherence interpretation"}
        try:
            inventory, report["sources"] = source_receipts()
            for label in ("gamma", "half-b1"):
                try:
                    report["actual_inputs"][label] = import_actual(label, inventory["frozen_curvature_inputs"])
                except Exception as error:
                    report["actual_inputs"][label] = {"import_pass": False, "error_type": type(error).__name__,
                        "error": str(error), "traceback": traceback.format_exc()}
            report["complex_energy_oracle"] = convention_oracle()
        except Exception as error:
            report["first_failure"] = {"error_type": type(error).__name__, "error": str(error),
                                       "traceback": traceback.format_exc()}
        report["acceptance_pass"] = (len(report["actual_inputs"]) == 2 and
            all(result["import_pass"] for result in report["actual_inputs"].values()) and
            report.get("complex_energy_oracle", {}).get("proof_pass", False) and "first_failure" not in report)
        text = json.dumps(report, indent=2, allow_nan=False) + "\n"
        output.write(text)
    print(text, end="")
    return 0 if report["acceptance_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
