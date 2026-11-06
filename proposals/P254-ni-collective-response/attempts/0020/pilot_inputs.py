"""Content-bound actual inputs for the registered P254 screened-response pilot."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from substrate_framework.qe_frozen_curvature import parse_frozen_curvature, verify_buffer
from substrate_framework.qe_response_bundle import load_bundle


HERE = Path(__file__).resolve().parent


def binding(path: Path, data: bytes | None = None) -> dict:
    raw = path.read_bytes() if data is None else data
    return {"path": str(path), "bytes": str(len(raw)), "sha256": hashlib.sha256(raw).hexdigest()}


def read_bound(record: dict) -> bytes:
    raw = Path(record["path"]).read_bytes()
    verify_buffer(raw, record)
    return raw


def load_inputs() -> tuple[dict, dict, dict]:
    inventory_raw = (HERE / "actual-input-inventory.json").read_bytes()
    inventory = json.loads(inventory_raw)
    bundles, frozen, receipts = {}, {}, {}
    for label, record in inventory["points"].items():
        bundle = load_bundle(record["bundle_root"])
        verify_buffer(bundle.buffers["manifest.json"], record["bundle_manifest"])
        # Bind producer identity as well as the actual loaded bundle buffers.
        read_bound(record["case"])
        read_bound(record["owner_result"])
        bundles[label] = bundle
        receipts[label] = {"bundle_buffers": [binding(bundle.root / name, data)
                                               for name, data in bundle.buffers.items()],
                           "case": record["case"], "owner_result": record["owner_result"]}
        if label == "minus":
            continue
        bg, meta = bundle.background, bundle.background.metadata
        identity = {"nat": meta["nat"], "iq": meta["iq"], "label": label,
                    "q_cart_2pi_alat": bg.array("q_cart_2pi_alat").tolist(),
                    "bundle_manifest": record["bundle_manifest"],
                    "cell": {"alat_bohr": meta["alat"], "omega_bohr3": meta["omega"],
                             "at_columns_alat": bg.array("at").tolist(),
                             "bg_columns_2pi_alat": bg.array("bg").tolist()},
                    "applicability_flags": {key: meta[key] for key in
                        ("okvan", "nlcc_any", "noncolin", "do_comp_mt", "do_cutoff_2D", "lda_plus_u")}}
        value = parse_frozen_curvature(read_bound(record["frozen_dynmat"]),
            read_bound(record["frozen_patterns"]), dynmat_source=record["frozen_dynmat"],
            patterns_source=record["frozen_patterns"], identity=identity)
        for bare in bundle.bare:
            mode = bare.metadata["mode"]
            if np.max(np.abs(value.patterns[:, mode - 1] - bare.array("Cartesian_complex_pattern"))) > 1e-12:
                raise ValueError("actual own-producer XML/bundle displacement patterns disagree")
        frozen[label] = value
        receipts[label]["frozen_curvature"] = value.record()
    return bundles, frozen, {"inventory": binding(HERE / "actual-input-inventory.json", inventory_raw),
                             "points": receipts}


def loaded_sources() -> list[dict]:
    """Bind actual imported implementation files, including the lazy XML consumer."""
    records = {}
    for name, module in list(sys.modules.items()):
        if name.startswith(("substrate_framework", "catdft")):
            filename = getattr(module, "__file__", None)
            if filename and Path(filename).suffix == ".py":
                path = Path(filename).resolve()
                records[str(path)] = binding(path)
    for path in (Path(__file__).resolve(), HERE / "run_screened_pilot.py", HERE / "pilot_checkpoints.py"):
        records[str(path)] = binding(path)
    return list(records.values())
