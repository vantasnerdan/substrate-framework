"""Atomic state and exact-byte completed-column recovery; no scientific recomputation."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from pilot_inputs import binding


def save_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    if path.exists():
        os.replace(path, path.with_suffix(path.suffix + ".prev"))
    os.replace(temporary, path)


def column_keys(settings: dict) -> tuple[str, ...]:
    return tuple(f"{point}-z{frequency}-column{axis}"
                 for point in settings["q"]
                 for frequency in range(len(settings["frequency_ry"]))
                 for axis in range(6))


def recover_columns(target: Path, fingerprint: str, settings: dict) -> dict:
    """Recover only fully committed, same-generation original-byte columns.

    A receipt is the commit boundary: orphan .tmp/NPZ files are not results.
    The result summary is not the column index, so an interrupted summary update
    cannot erase a committed column. Original absolute archive ownership remains
    bound; this is not permission to move/re-pin a scientific generation.
    """
    recovered = {}
    for name in column_keys(settings):
        receipt_path, archive_path = target / (name + ".json"), target / (name + ".npz")
        if not receipt_path.exists():
            continue
        receipt = json.loads(receipt_path.read_bytes())
        if (receipt["name"] != name or receipt["fingerprint"] != fingerprint
                or receipt["archive"] != binding(archive_path)):
            raise ValueError(f"column checkpoint identity/bytes differ: {name}")
        with np.load(archive_path, allow_pickle=False) as stored:
            column = stored["coupled_curvature_column"].copy()
        if column.shape != (6,) or not np.iscomplexobj(column) or not np.all(np.isfinite(column)):
            raise ValueError(f"invalid physical coupled-curvature column: {name}")
        recovered[name] = (receipt, column)
    return recovered


def checkpoint_report(target: Path, report: dict, recovered: dict, elapsed: float,
                      prior_elapsed: float) -> None:
    """Publish a current summary without forgetting previously committed columns."""
    by_name = {name: value[0] for name, value in recovered.items()}
    by_name.update({receipt["name"]: receipt for receipt in report["completed_columns"]})
    report["completed_columns"] = list(by_name.values())
    report["elapsed_seconds"] = elapsed
    report["cumulative_elapsed_seconds"] = prior_elapsed + elapsed
    save_json(target / "result.json", report)
