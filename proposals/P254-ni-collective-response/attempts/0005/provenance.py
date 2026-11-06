"""Hash declared inputs, Python sources, imported modules, and mapped runtime libraries."""
import hashlib
import os
from pathlib import Path
import platform
import sys

BASE = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest():
    import numpy
    import scipy
    import sympy
    files = list(BASE.rglob("*.py")) + [BASE / "preregistration.json", BASE / "derivation.txt"]
    files += [BASE / "execution-protocol.json"]
    files += [BASE.parents[1] / "W2-readiness.json"]
    files += [BASE.parents[1] / "sources" / name for name in (
        "dal-corso2001-provenance.json", "electronic-model-source-correction.json")]
    return {
        "inputs_and_code_sha256": {str(p): sha(p) for p in sorted(files)},
        "python": {"executable": sys.executable, "sha256": sha(sys.executable),
                   "version": sys.version, "platform": platform.platform()},
        "packages": {m.__name__: {"version": m.__version__, "file": m.__file__,
                                  "sha256": sha(m.__file__)} for m in (numpy, scipy, sympy)},
        "affinity": sorted(os.sched_getaffinity(0)),
        "environment": {k: os.environ.get(k) for k in (
            "PYTHONPATH", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "INVOCATION_ID", "XDG_RUNTIME_DIR", "DBUS_SESSION_BUS_ADDRESS")},
        "command": [sys.executable] + sys.argv,
    }


def loaded_runtime():
    paths = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
             if getattr(m, "__file__", None) and Path(m.__file__).is_file()}
    for row in Path("/proc/self/maps").read_text().splitlines():
        fields = row.split(maxsplit=5)
        if len(fields) == 6 and fields[5].startswith("/") and Path(fields[5]).is_file():
            paths.add(Path(fields[5]).resolve())
    return {str(p): sha(p) for p in sorted(paths)}


def fingerprint(provenance):
    return {key: provenance[key] for key in ("inputs_and_code_sha256", "python", "packages")}
