"""Describe, then explicitly build the PRIVATE original-XClib PBE wrapper.

No build by default. Root saves --describe JSON, closes/adopts its input set,
then passes that document as --admission with --build. Only A25/build/internal-pbe
is created; an existing target is refused. No QE/old-build changes or fallback.
The unchanged archive supplies every PBE routine. Its non-PIC internal data
references are localized; only the unused BEEF adjunct is recompiled PIC,
because the original beefun.o contains an external stdout PC32 relocation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent
QE = Path("/opt/qebuild/q-e-qe-7.5/XClib")
BUILD = HERE / "build" / "internal-pbe"
FC = "/usr/bin/gfortran"
CC = "/usr/bin/cc"
MPI = "/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so"
BLAS = "/usr/lib/x86_64-linux-gnu/libopenblas.so"
ARCHIVE_SHA256 = "a8f7bdfdfd3c94fa460643511aca9d7e6493ee41a4c7c8466d9b0f19eabab245"


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return {"path": str(path), "bytes": str(len(data)),
            "sha256": hashlib.sha256(data).hexdigest()}


def fact(argv):
    return subprocess.run(argv, check=True, text=True, capture_output=True, timeout=10).stdout


def dependencies(path):
    text = fact(["/usr/bin/ldd", str(path)])
    if "not found" in text:
        raise RuntimeError(f"unresolved dynamic dependency: {text}")
    return sorted({Path(match).resolve() for match in re.findall(r"(/[^\s()]+)", text)})


def commands():
    # All generated modules, object files, copied inputs and ELF stay private.
    return [
        [CC, "-O3", "-fPIC", "-c", "inputs/beefun.c", "-o", "beefun.pic.o"],
        [FC, "-O3", "-fPIC", "-Iinputs", "-J.", "-c", "inputs/internal_pbe_bridge.f90",
         "-o", "internal_pbe_bridge.o"],
        [FC, "-shared", "-Wl,-z,defs", "-Wl,-z,text", "-Wl,-Bsymbolic",
         "-Wl,--exclude-libs,ALL", "internal_pbe_bridge.o", "beefun.pic.o",
         "inputs/xc_lib.a", MPI, BLAS, "-o", "libinternal_pbe.so"],
    ]


def describe():
    inventory = json.loads((HERE / "source-access-inventory.json").read_text())
    registered = inventory["new_sources_and_compiled_interface_inputs"]
    source_paths = {QE / name for name in ("xc_lib.a", "beefun.c", "beefleg.h", "pbecor.h")}
    source_paths.update(Path(item["path"]) for item in registered
                        if Path(item["path"]).parent == QE and item["path"].endswith(".mod"))
    source_paths.add(HERE / "internal_pbe_bridge.f90")
    inputs = set(source_paths) | {Path(__file__), HERE / "source-access-inventory.json"}
    # Dependency-only preprocessing: records real system headers, not a guessed list.
    header_text = fact([CC, "-M", "-O3", "-fPIC", str(QE / "beefun.c")])
    inputs.update(Path(p) for p in header_text.replace("\\\n", " ").split(":", 1)[1].split())
    tools = {Path(FC), Path(CC), Path("/usr/bin/as"), Path("/usr/bin/ld"), Path("/usr/bin/ldd")}
    for compiler, name in ((FC, "f951"), (CC, "cc1"), (FC, "collect2")):
        tools.add(Path(fact([compiler, f"-print-prog-name={name}"]).strip()))
    inputs.update(tools)
    # GCC startup/archive inputs and the real MPI/BLAS/runtime library closure.
    libraries = {Path(MPI), Path(BLAS)}
    for name in ("crti.o", "crtn.o", "crtbeginS.o", "crtendS.o", "libgfortran.so",
                 "libm.so", "libmvec.so.1", "libgcc.a", "libgcc_s.so", "libgcc_s.so.1",
                 "libc.so", "libc_nonshared.a"):
        value = fact([FC, f"-print-file-name={name}"]).strip()
        if value == name:
            raise RuntimeError(f"compiler did not resolve required link input {name}")
        libraries.add(Path(value))
    inputs.update(libraries)
    # ldd is a script; linker scripts/startup objects are not dynamic ELFs.
    for path in (Path(FC), Path(CC), *[t for t in tools if t.name != "ldd"],
                 Path(MPI), Path(BLAS), Path(fact([FC, "-print-file-name=libgfortran.so"]).strip())):
        inputs.update(dependencies(path))
    bindings = [binding(path) for path in sorted({p.resolve() for p in inputs})]
    by_path = {item["path"]: item for item in bindings}
    for item in registered:
        resolved = str(Path(item["path"]).resolve())
        if resolved in by_path and any(by_path[resolved][k] != item[k] for k in ("bytes", "sha256")):
            raise RuntimeError(f"registered source changed: {resolved}")
    if by_path[str((QE / "xc_lib.a").resolve())]["sha256"] != ARCHIVE_SHA256:
        raise RuntimeError("not the original Native XClib archive")
    return {"schema": "internal-pbe-build-v1", "cwd": str(BUILD), "commands": commands(),
            "private_copies": [str(p.resolve()) for p in sorted(source_paths)],
            "inputs": bindings, "compiler_versions": {FC: fact([FC, "--version"]), CC: fact([CC, "--version"])},
            "scope": "original internal PBE archive plus original BEEF PIC adjunct; no numerical execution"}


def build(plan, admission, seconds):
    if plan != admission:
        raise ValueError("current build inputs/argv differ from Root's admitted description")
    if not 0 < seconds <= 120:
        raise ValueError("build deadline must be in (0,120] seconds")
    BUILD.parent.mkdir(exist_ok=True)
    BUILD.mkdir()  # never overwrite/rebuild while any consumer could be live
    (BUILD / "inputs").mkdir()
    admitted = {item["path"]: item for item in plan["inputs"]}
    for name in plan["private_copies"]:
        data = Path(name).read_bytes()
        expected = admitted[name]
        if str(len(data)) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
            raise ValueError(f"input changed before private copy: {name}")
        (BUILD / "inputs" / Path(name).name).write_bytes(data)
    deadline = time.monotonic() + seconds
    executed = []
    for index, argv in enumerate(plan["commands"]):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise TimeoutError("cumulative private-build deadline exhausted")
        start = time.monotonic()
        result = subprocess.run(argv, cwd=BUILD, capture_output=True, timeout=remaining)
        (BUILD / f"step-{index}.stdout").write_bytes(result.stdout)
        (BUILD / f"step-{index}.stderr").write_bytes(result.stderr)
        executed.append({"argv": argv, "seconds": time.monotonic()-start, "returncode": result.returncode})
        result.check_returncode()
    library = binding(BUILD / "libinternal_pbe.so")
    library["dependencies"] = [binding(path) for path in dependencies(library["path"])]
    # Every actual linked runtime must already be in the Root-closed input set.
    for item in library["dependencies"]:
        if item != admitted.get(item["path"]):
            raise RuntimeError(f"unadmitted linked runtime dependency: {item['path']}")
    os.chmod(library["path"], 0o444)
    receipt = {"plan": plan, "executed": executed, "library": library,
               "numerical_execution": False}
    (BUILD / "build-receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    (BUILD / "library-binding.json").write_text(json.dumps(library, indent=2)+"\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--describe", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--admission", type=Path)
    parser.add_argument("--seconds", type=float, default=120)
    args = parser.parse_args()
    if args.build and (args.describe or args.admission is None):
        parser.error("--build requires --admission and excludes --describe")
    plan = describe()
    result = build(plan, json.loads(args.admission.read_text()), args.seconds) if args.build else plan
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
