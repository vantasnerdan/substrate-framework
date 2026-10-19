"""Explicit local archive construction and isolated scientific reproduction.

Nothing executes on extraction. Run --help for the opt-in commands.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib import metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from typing import Any
import zipfile

HERE = Path(__file__).resolve().parent
BUNDLE = "journal_flux_axial_mobility_20261019"
PACKAGE = "journal_flux_axial_mobility_package_20261019.py"
MANIFEST = "journal_flux_axial_mobility_manifest_20261019.json"
CARRIER = "rich_journal_flux_carrier_20261019"
BAND = "rich_journal_baryon_band_20261019"
COMPACT = "rich_journal_compact_hierarchy_20261019"
FIGURES = "journal_flux_axial_mobility_figures_20261019"
DATA = "journal_flux_axial_mobility_data_20261019.json"
TABLES = "journal_flux_axial_mobility_tables_20261019.tex"
BIB = "journal_campaign_20261019.bib"
REQUIREMENTS = "journal_flux_axial_mobility_requirements_20261019.txt"
DISTRIBUTIONS = ("numpy", "scipy", "threadpoolctl", "matplotlib")
REFERENCE_ENVIRONMENT = {
    "python": "3.12.3", "numpy": "2.5.3", "scipy": "1.18.1",
    "threadpoolctl": "3.7.0", "matplotlib": "3.11.2",
}
SCIENTIFIC_SOURCES = (
    CARRIER + ".py", BAND + ".py", COMPACT + ".py", "rich_nucleon_hierarchy_20261019.py",
    "rich_composite_local_flux_binding_20261017_model.py",
    "rich_colour_spin_root_20261019.py",
)
CALCULATIONS = (CARRIER + "_result.json", BAND + "_result.json", COMPACT + "_result.json")
FIGURE_OUTPUTS = (FIGURES + ".pdf", DATA, TABLES)
PAPER_INPUTS = (BUNDLE + ".tex", BIB, "journal_campaign_20261019_references.json")
# No directory traversal/globbing: changing publication scope requires editing this list.
PAYLOAD = (
    PACKAGE, "journal_flux_axial_mobility_REPRODUCE_20261019.md", REQUIREMENTS,
    *SCIENTIFIC_SOURCES, *CALCULATIONS, FIGURES + ".py", *FIGURE_OUTPUTS,
    *PAPER_INPUTS, BUNDLE + ".pdf",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def environment() -> dict[str, str]:
    result = {"python": platform.python_version()}
    for name in DISTRIBUTIONS:
        try:
            result[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            result[name] = "not installed"
    return result


def read_payload() -> dict[str, bytes]:
    missing = [name for name in PAYLOAD if not (HERE / name).is_file()]
    if missing:
        raise ValueError("Missing final package prerequisites:\n  " + "\n  ".join(missing))
    if any((HERE / name).is_symlink() for name in PAYLOAD):
        raise ValueError("Allowlisted files must be regular files, not symlinks")
    return {name: (HERE / name).read_bytes() for name in PAYLOAD}


def check_bindings(payload: dict[str, bytes]) -> None:
    carrier = json.loads(payload[CALCULATIONS[0]])
    compact = json.loads(payload[CALCULATIONS[2]])
    figures = json.loads(payload[DATA])
    expected_sources = set(SCIENTIFIC_SOURCES) - {BAND + ".py", COMPACT + ".py"}
    if set(carrier["source_sha256"]) != expected_sources:
        raise ValueError("Carrier provenance file set differs from the explicit scientific allowlist")
    expected_compact = (expected_sources - {"rich_colour_spin_root_20261019.py"}) | {
        COMPACT + ".py", CALCULATIONS[0],
    }
    if set(compact["source_sha256"]) != expected_compact:
        raise ValueError("Compact hierarchy provenance differs from the explicit scientific allowlist")
    expected_figures = expected_sources | expected_compact | {
        *CALCULATIONS, BAND + ".py", FIGURES + ".py",
    }
    if set(figures["source_sha256"]) != expected_figures:
        raise ValueError("Figure provenance file set differs from the explicit allowlist")
    for label, bindings in (("carrier", carrier["source_sha256"]),
                            ("compact hierarchy", compact["source_sha256"]),
                            ("figure/table data", figures["source_sha256"])):
        for name, expected in bindings.items():
            if sha256(payload[name]) != expected:
                raise ValueError(f"Stale {label} source binding: {name}; regenerate its dependent output")
    if figures["carrier_recorded_source_sha256"] != carrier["source_sha256"]:
        raise ValueError("Figure data and carrier recorded provenance disagree")
    # These private path roots are assembled so this check does not match itself.
    private_roots = [('/'.join(('', part, ''))).encode() for part in ("home", "root", "Users", "tmp")]
    for name, data in payload.items():
        if any(prefix in data for prefix in private_roots):
            raise ValueError(f"Private absolute path in public payload {name}; inspect rather than silently redact")


def build(output: Path) -> None:
    payload = read_payload()
    check_bindings(payload)
    output = output.resolve()
    if output in {(HERE / name).resolve() for name in (*PAYLOAD, MANIFEST)}:
        raise ValueError("Archive destination must not overwrite an input or its manifest")
    if not output.parent.is_dir():
        raise ValueError("Archive destination parent directory must already exist")
    manifest = {
        "schema": 1, "bundle_directory": BUNDLE, "hash_algorithm": "SHA256",
        "files": {name: {"sha256": sha256(data), "bytes": len(data)}
                  for name, data in payload.items()},
        "manifest_scope": "Every payload file; manifest excludes itself and the enclosing ZIP to avoid circular hashes.",
        "scientific_source_policy": "Scientific sources and reference results are byte-preserved; no public-copy transformations.",
        "provenance_labels": "Historical source-provenance labels are opaque identifiers, not required external services. All computation is local.",
        "reference_environment": REFERENCE_ENVIRONMENT,
        "construction_environment": environment(),
        "reference_toolchain": {
            "pdflatex": "pdfTeX 3.141592653-2.6-1.40.25 (TeX Live 2023/Debian)",
            "bibtex": "BibTeX 0.99d (TeX Live 2023/Debian)",
        },
        "reproduction": {
            "entrypoint": PACKAGE, "full_mode": "all", "frozen_input_mode": "paper",
            "work_directory": "Required new directory separate from the frozen bundle",
            "numeric_comparison": {"relative_tolerance": 2e-7, "absolute_tolerance": 2e-8},
        },
    }
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix=output.name + ".", suffix=".tmp",
                                         dir=output.parent, delete=False) as stream:
            temporary = Path(stream.name)
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=9) as archive:
            for name, data in sorted({**payload, MANIFEST: manifest_bytes}.items()):
                entry = zipfile.ZipInfo(BUNDLE + "/" + name, date_time=(2026, 10, 19, 0, 0, 0))
                entry.compress_type = zipfile.ZIP_DEFLATED
                entry.external_attr = 0o100644 << 16
                archive.writestr(entry, data)
        temporary.replace(output)
        (HERE / MANIFEST).write_bytes(manifest_bytes)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    print(json.dumps({"archive": str(output), "sha256": sha256(output.read_bytes()),
                      "payload_files": len(payload), "archive_members": len(payload) + 1}, indent=2))


def verify() -> dict[str, bytes]:
    if not (HERE / MANIFEST).is_file():
        raise ValueError("Frozen manifest missing; run build only after all final paper prerequisites exist")
    manifest = json.loads((HERE / MANIFEST).read_text())
    if manifest.get("schema") != 1 or manifest.get("bundle_directory") != BUNDLE:
        raise ValueError("Unsupported manifest schema or bundle identity")
    if set(manifest["files"]) != set(PAYLOAD):
        raise ValueError("Manifest file set differs from the entrypoint's explicit allowlist")
    payload = read_payload()
    for name, data in payload.items():
        recorded = manifest["files"][name]
        if len(data) != recorded["bytes"] or sha256(data) != recorded["sha256"]:
            raise ValueError(f"Frozen payload mismatch: {name}")
    check_bindings(payload)
    return payload


def observable_summary(carrier: dict[str, Any], band: dict[str, Any],
                       compact: dict[str, Any]) -> dict[str, Any]:
    """Compare scientific observables, not timings, solver calls or PDF timestamps."""
    keys = ("loop_representation", "total_representation_degree_cutoff", "mass_m", "I",
            "energy", "source_only_energy", "axial_neutral_tau3_sigma3",
            "vacuum_seed_overlap_Z0", "probability_law_defect", "inelastic_axial",
            "canonical_degree3_exterior")
    mergers = [{
        "mass_m": row["mass_m"],
        "compact_energy": row["compact_N6"]["energy"],
        "one_N3_energy": row["one_N3_cell"]["energy"],
        "separated_pn_I0_T1_energy": row["separated_pn_I0_T1_energy"],
        "compact_minus_separated": row["compact_minus_separated"],
        "one_source_only_compact_minus_separated": row["one_source_only_compact_minus_separated"],
        "leading_compact_minus_separated": row["leading_compact_minus_separated"],
        "physical_sign": row["physical_sign"],
    } for row in compact["comparisons"]]
    return {"carrier_rows": [{key: row[key] for key in keys if key in row}
                             for row in carrier["rows"]], "band": band,
            "compact_merger_comparisons": mergers}


def compare(reference: Any, actual: Any, path: str = "") -> list[dict[str, Any]]:
    if isinstance(reference, dict) and isinstance(actual, dict):
        if reference.keys() != actual.keys():
            return [{"field": path, "error": "dictionary keys differ"}]
        return [item for key in reference for item in compare(reference[key], actual[key], path + "/" + key)]
    if isinstance(reference, list) and isinstance(actual, list):
        if len(reference) != len(actual):
            return [{"field": path, "error": "list lengths differ"}]
        return [item for i, (left, right) in enumerate(zip(reference, actual))
                for item in compare(left, right, path + "/" + str(i))]
    if (isinstance(reference, (int, float)) and isinstance(actual, (int, float))
            and not isinstance(reference, bool) and not isinstance(actual, bool)):
        equal = math.isclose(reference, actual, rel_tol=2e-7, abs_tol=2e-8)
    else:
        equal = type(reference) is type(actual) and reference == actual
    return [] if equal else [{"field": path, "reference": reference, "actual": actual}]


def run_child(argv: list[str], work: Path, log_name: str, child_env: dict[str, str]) -> None:
    print("Running " + log_name, flush=True)
    with (work / (log_name + ".stdout.txt")).open("w") as stdout, \
            (work / (log_name + ".stderr.txt")).open("w") as stderr:
        result = subprocess.run(argv, cwd=work, env=child_env, stdout=stdout, stderr=stderr, check=False)
    if result.returncode:
        raise ValueError(f"{log_name} exited {result.returncode}; inspect its stdout/stderr files in the work directory")


def reproduce(work: Path, mode: str) -> None:
    payload = verify()
    actual_environment = environment()
    missing = [name for name in DISTRIBUTIONS if actual_environment[name] == "not installed"]
    missing += [name for name in ("pdflatex", "bibtex") if shutil.which(name) is None]
    if missing:
        raise ValueError("Missing reproduction prerequisites: " + ", ".join(missing))
    work = work.resolve()
    if work == HERE or work in HERE.parents or HERE in work.parents:
        raise ValueError("Work directory must be separate from, not inside or above, the frozen bundle")
    if work.exists():
        raise ValueError("Work directory already exists; choose a new path to preserve earlier runs")
    work.mkdir(parents=True)
    inputs = (*SCIENTIFIC_SOURCES, FIGURES + ".py", *PAPER_INPUTS)
    if mode == "paper":
        inputs += CALCULATIONS
    for name in inputs:
        (work / name).write_bytes(payload[name])
    child_env = dict(os.environ)
    child_env.update({name: "1" for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")})
    child_env["MPLCONFIGDIR"] = str(work / ".matplotlib")
    child_env["PYTHONNOUSERSITE"] = "1"
    child_env.pop("PYTHONPATH", None)
    report: dict[str, Any] = {"mode": mode, "success": False, "environment": actual_environment,
                              "reference_environment": REFERENCE_ENVIRONMENT,
                              "reference_files_modified": False}
    try:
        if mode == "all":
            for supplier in (CARRIER, BAND, COMPACT):
                run_child([sys.executable, supplier + ".py"], work, supplier, child_env)
            reference = observable_summary(*(json.loads(payload[name]) for name in CALCULATIONS))
            actual = observable_summary(*(json.loads((work / name).read_text()) for name in CALCULATIONS))
            differences = compare(reference, actual)
            report["scientific_comparison"] = {"relative_tolerance": 2e-7, "absolute_tolerance": 2e-8,
                                               "differences": differences, "passed": not differences}
            if differences:
                raise ValueError("Scientific observables differ; inspect reproduction_report.json")
        else:
            report["scientific_comparison"] = {"executed": False, "reason": "Frozen calculation inputs selected"}
        run_child([sys.executable, FIGURES + ".py"], work, "figures", child_env)
        latex = ["pdflatex", "-no-shell-escape", "-interaction=nonstopmode", "-halt-on-error", BUNDLE + ".tex"]
        run_child(latex, work, "pdflatex-1", child_env)
        run_child(["bibtex", BUNDLE], work, "bibtex", child_env)
        run_child(latex, work, "pdflatex-2", child_env)
        run_child(latex, work, "pdflatex-3", child_env)
        generated = (*CALCULATIONS, *FIGURE_OUTPUTS, BUNDLE + ".pdf")
        report["output_sha256"] = {name: sha256((work / name).read_bytes()) for name in generated}
        # Detect accidental writes back into the frozen reference, not just successful subprocesses.
        verify()
        report["success"] = True
    except (ValueError, OSError) as error:
        report["error"] = str(error)
        raise
    finally:
        (work / "reproduction_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"success": True, "mode": mode, "work_directory": str(work),
                      "report": "reproduction_report.json"}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("status", help="List exact required files and missing final prerequisites; write nothing")
    builder = commands.add_parser("build", help="Require complete final inputs, then write manifest and explicit-file ZIP")
    builder.add_argument("--output", type=Path, default=HERE / (BUNDLE + ".zip"))
    commands.add_parser("verify", help="Check frozen file hashes and source/result/figure bindings; no calculations")
    runner = commands.add_parser("reproduce", help="Opt in to isolated calculations/figures/LaTeX; never overwrite frozen files")
    runner.add_argument("--work-dir", type=Path, required=True)
    runner.add_argument("--mode", choices=("all", "paper"), default="all")
    args = parser.parse_args()
    try:
        if args.command == "status":
            missing = [name for name in PAYLOAD if not (HERE / name).is_file()]
            print(json.dumps({"required_files": list(PAYLOAD), "missing_prerequisites": missing,
                              "ready_for_build": not missing}, indent=2))
            if missing:
                raise SystemExit(1)
        elif args.command == "build":
            build(args.output)
        elif args.command == "verify":
            verify()
            print(json.dumps({"verified_payload_files": len(PAYLOAD), "source_bindings": "matched"}))
        else:
            reproduce(args.work_dir, args.mode)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()
