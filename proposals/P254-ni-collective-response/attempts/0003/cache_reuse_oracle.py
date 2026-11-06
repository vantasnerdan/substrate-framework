"""Exercise recovery through the actual native replay cache, without rerunning QE."""

from pathlib import Path

from native_qe_replay import replay_serialized_operators
from oracle_artifacts import sha256


def cache_reuse_block(ledger, result_path, native_directory, first):
    receipt_path = Path(first["frequency_file"]).parent.parent / "receipt.json"
    paths = [Path(path) for path in first["artifact_sha256"]] + [receipt_path]
    before = {str(path): {"sha256":sha256(path), "mtime_ns":path.stat().st_mtime_ns} for path in paths}
    recovered = replay_serialized_operators(result_path,native_directory)
    after = {str(path): {"sha256":sha256(path), "mtime_ns":path.stat().st_mtime_ns} for path in paths}
    ledger.check("completed native replay recovered without changing outputs or mtimes", before == after)
    ledger.check("recovered native replay retains bound original inputs and decoded spectra", recovered == first)
    return {"receipt":str(receipt_path), "unchanged_artifacts":before,
            "scope":"Actual native completed-cache recovery; no native process or replacement attempt on reuse."}
