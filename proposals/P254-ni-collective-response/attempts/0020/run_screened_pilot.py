#!/usr/bin/env python3
"""Actual P254 same-frequency screened ionic curvature; bounded, column-checkpointed."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback

import numpy as np
import scipy

from substrate_framework.pw_independent_response import PWIndependentResponse
from substrate_framework.pw_screened_response import solve_screened_response
from substrate_framework.qe_pw_response_frames import build_pw_response_models
from substrate_framework.spin_pbe_kernel import apply_spin_pbe_kernel
from pilot_inputs import HERE, binding, load_inputs, loaded_sources
from pilot_checkpoints import checkpoint_report, column_keys, recover_columns, save_json


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def interrupted(signum, _frame) -> None:
    raise InterruptedError(f"received signal {signum}; settled column checkpoints retained")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    target = args.output_dir.resolve()
    if target.exists() and not args.resume:
        parser.error("output directory exists; use explicit --resume with unchanged source identity")
    target.mkdir(parents=True, exist_ok=args.resume)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    prior_elapsed = 0.
    prior_completed = False
    if args.resume:
        prior_result = json.loads((target / "result.json").read_bytes())
        prior_elapsed = float(prior_result["cumulative_elapsed_seconds"])
        prior_completed = bool(prior_result["completed_columns"])
    started = time.monotonic()
    registration = json.loads((HERE / "registration.json").read_bytes())
    settings = registration["actual_pilot"]
    recovered = {}
    admitted = False
    allocation_deadline = time.time() + settings["total_wall_cap_seconds"]
    report = {"started_at": now(), "status": "running", "completed_columns": [],
              "registration": binding(HERE / "registration.json"),
              "preregistration_commit": "eacf18e0", "scope": registration["scientific_scope"],
              "settings": settings, "command": [sys.executable, *sys.argv],
              "environment": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                  "cpu_affinity": [str(cpu) for cpu in sorted(os.sched_getaffinity(0))],
                  "thread_variables": {name: os.environ.get(name) for name in
                       ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}}}

    def event(kind: str, **values) -> None:
        record = {"at": now(), "elapsed_seconds": time.monotonic() - started,
                  "event": kind, **values}
        with (target / "telemetry.jsonl").open("a") as stream:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
            stream.flush()
        print(json.dumps(record, allow_nan=False), flush=True)

    def check_budget() -> None:
        elapsed = prior_elapsed + time.monotonic() - started
        if elapsed > settings["total_wall_cap_seconds"] or time.time() >= allocation_deadline:
            raise TimeoutError("registered cumulative 600s pilot allocation exhausted")
        if not prior_completed and not report["completed_columns"] and elapsed > settings["first_completed_column_due_seconds_from_launch"]:
            raise TimeoutError("no useful completed column within registered 150s; Root method decision required")

    def checkpoint() -> None:
        checkpoint_report(target, report, recovered, time.monotonic() - started, prior_elapsed)

    try:
        event("pilot_started", output_dir=str(target))
        report["implementation"] = loaded_sources()
        identity = {"implementation": report["implementation"], "settings": settings,
                    "input_inventory": binding(HERE / "actual-input-inventory.json")}
        fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        identity_path = target / "identity.json"
        if args.resume:
            previous = json.loads(identity_path.read_bytes())
            if previous["fingerprint"] != fingerprint:
                raise ValueError("resume source/input/settings identity differs; no mismatched column reuse")
            allocation_deadline = float(previous["allocation_deadline_unix"])
            if not np.isfinite(allocation_deadline):
                raise ValueError("invalid original allocation deadline")
            recovered = recover_columns(target, fingerprint, settings)
            prior_completed = bool(recovered)
        else:
            save_json(identity_path, {"fingerprint": fingerprint,
                                     "allocation_deadline_unix": allocation_deadline, **identity})
        admitted = True
        bundles, frozen, receipts = load_inputs()
        report["actual_inputs"] = receipts
        # The XML consumer was imported before source inventory, not lazily after it.
        event("actual_inputs_bound")
        models = None
        if len(recovered) != len(column_keys(settings)):
            check_budget()
            build_started = time.monotonic()
            models = build_pw_response_models(bundles["gamma"], bundles["half"], bundles["minus"])
            report["frame_build_seconds"] = time.monotonic() - build_started
            report["reference"] = models["gamma"].metadata
            event("frames_constructed", seconds=report["frame_build_seconds"])
        else:
            report["reference"] = prior_result["reference"]
            report["frame_build_seconds"] = 0.
        checkpoint()
        matrices = {}
        for point in settings["q"]:
            kernel = lambda density, point=point: apply_spin_pbe_kernel(bundles[point], density)
            for frequency_index, (real, imaginary) in enumerate(settings["frequency_ry"]):
                response = (PWIndependentResponse(models[point], complex(real, imaginary))
                            if models is not None else None)
                columns = []
                key = f"{point}-z{frequency_index}"
                for axis in range(6):
                    name = f"{key}-column{axis}"
                    receipt_path, archive_path = target / (name + ".json"), target / (name + ".npz")
                    if name in recovered:
                        _, column = recovered[name]
                        columns.append(column)
                        event("matched_column_reused", name=name)
                        continue
                    check_budget()
                    event("column_started", name=name)
                    displacement = np.eye(6, dtype=complex)[axis]
                    column_started = time.monotonic()
                    remaining = min(allocation_deadline - time.time(),
                                    settings["total_wall_cap_seconds"] - (prior_elapsed + column_started - started))
                    value = solve_screened_response(response, kernel, displacement,
                        rtol=settings["rtol"], max_matvec=settings["max_matvec_per_column"],
                        wall_seconds=min(settings["column_wall_cap_seconds"], remaining))
                    curvature = frozen[point].cartesian_matrix[:, axis] + value.ionic_gradient
                    if not np.all(np.isfinite(curvature)):
                        raise ValueError("nonfinite coupled curvature cannot be checkpointed")
                    temporary = archive_path.with_suffix(".npz.tmp")
                    with temporary.open("wb") as stream:
                        np.savez_compressed(stream, density=value.density, ionic_gradient=value.ionic_gradient,
                            induced_potential=value.induced_potential, coupled_curvature_column=curvature,
                            frozen_curvature_column=frozen[point].cartesian_matrix[:, axis])
                    os.replace(temporary, archive_path)
                    receipt = {"name": name, "point": point, "frequency_ry": [real, imaginary],
                               "axis_zero_based": str(axis), "finished_at": now(),
                               "wall_seconds": time.monotonic() - column_started,
                               "fingerprint": fingerprint, "diagnostics": value.diagnostics,
                               "archive": binding(archive_path)}
                    save_json(receipt_path, receipt)
                    report["completed_columns"].append(receipt)
                    columns.append(curvature)
                    checkpoint()
                    event("column_completed", name=name, seconds=receipt["wall_seconds"],
                          completed=str(len(report["completed_columns"])), diagnostics=value.diagnostics)
                    if len(report["completed_columns"]) == 1 and not prior_completed and prior_elapsed + time.monotonic() - started > settings["first_completed_column_due_seconds_from_launch"]:
                        raise TimeoutError("first completed column arrived after150s; saved result retained, further work stopped")
                    check_budget()
                matrices[key] = np.column_stack(columns)
                event("matrix_completed", name=key)
        archive = target / "coupled-matrices.npz"
        temporary = archive.with_suffix(".npz.tmp")
        with temporary.open("wb") as stream:
            np.savez_compressed(stream, **matrices)
        if archive.exists():
            os.replace(archive, archive.with_suffix(".npz.prev"))
        os.replace(temporary, archive)
        report["matrix_archive"] = binding(archive)
        report["status"] = "complete"
        report["matrix_definition"] = "Cfrozen + physical delta<g>/delta_u; Ry/bohr2, unmassweighted, no ionic kinetic term"
    except Exception as error:
        report["status"] = "stopped_requires_root_decision"
        report["failure"] = {"type": type(error).__name__, "message": str(error),
                             "traceback": traceback.format_exc(),
                             "diagnostics": getattr(error, "diagnostics", None)}
    finally:
        report["finished_at"] = now()
        if admitted:
            checkpoint()
        else:
            save_json(target / "resume-refusal.json", report)
        if report["status"] == "complete":
            (target / ".gate_done_screened_response").write_text(report["finished_at"] + "\n")
        event("pilot_finished", status=report["status"], completed=str(len(report["completed_columns"])))
        print("[dft-done]" if report["status"] == "complete" else "[dft-fail]",
              "p254-screened-pilot-0020", str(target / ("result.json" if admitted else "resume-refusal.json")), flush=True)
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
