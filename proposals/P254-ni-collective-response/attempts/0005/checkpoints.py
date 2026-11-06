"""Atomically persist checked blocks and require exact source/runtime binding on resume."""
import json
import os
from provenance import fingerprint, loaded_runtime, sha


def atomic(path, value):
    temporary = path.with_suffix(path.suffix+".tmp")
    with temporary.open("w") as stream:
        json.dump(value, stream, indent=2, default=str)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def save(directory, block, result, provenance):
    atomic(directory / (block+".checkpoint.json"),
           {"block": block, "fingerprint": fingerprint(provenance), "result": result,
            "loaded_runtime_sha256": loaded_runtime()})


def recover(directory, block, provenance):
    path = directory / (block+".checkpoint.json")
    if not path.exists():
        return None
    checkpoint = json.loads(path.read_text())
    assert checkpoint["fingerprint"] == fingerprint(provenance), "checkpoint input/runtime binding changed"
    assert checkpoint["block"] == block
    for runtime_path, expected in checkpoint["loaded_runtime_sha256"].items():
        assert sha(runtime_path) == expected, ("checkpoint runtime file changed", runtime_path)
    return {"result": checkpoint["result"], "source": str(path), "sha256": sha(path)}
