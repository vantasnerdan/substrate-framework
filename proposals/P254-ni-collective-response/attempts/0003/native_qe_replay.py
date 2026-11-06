"""Re-diagonalize serialized PH matrices through the native q2r/matdyn path.

This is a linear-operator oracle, not another DFT simulation. PH's original
OMEGA tags refer to its in-memory matrix, whereas q2qstar_ph:118-140 writes
an additional symmetry average. Never identify those two finite-precision
objects. Keep raw inputs immutable and use no ASR or polar correction.
"""

import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import uuid

from oracle_artifacts import save_receipt, sha256

HERE = Path(__file__).resolve().parent
BINARIES = {
    "q2r": ("/opt/qebuild/q-e-qe-7.5/bin/q2r.x",
            "0c5f75d35cca50cae60a25dbd8f33a4437f5c95bf45bffda635229f94089b4e9"),
    "matdyn": ("/opt/qebuild/q-e-qe-7.5/bin/matdyn.x",
               "32d735261925e78b5cd9c548f063d242d841b04d80dcf8145ab8c8fbbee9f1eb"),
}


def fixed_fields(line, width=10):
    return [float(line[start:start + width])
            for start in range(0, len(line), width) if line[start:start + width].strip()]


def frequency_table(path, nat, grid):
    """Decode matdyn:776-779, including fixed-width adjacent negative fields."""
    lines = Path(path).read_text().splitlines()
    header = re.fullmatch(r"\s*&plot\s+nbnd=\s*(\d+),\s*nks=\s*(\d+)\s*/\s*", lines[0])
    indices = list(itertools.product(*(range(n) for n in grid)))
    if header is None or tuple(map(int, header.groups())) != (3 * nat, len(indices)):
        raise ValueError("matdyn frequency-table dimensions disagree with native input")
    rows = []
    cursor = 1
    for index in indices:
        coordinates = fixed_fields(lines[cursor][10:])
        if len(coordinates) != 4:
            raise ValueError("matdyn q-coordinate row must contain three coordinates and a weight")
        cursor += 1
        frequencies = []
        while len(frequencies) < 3 * nat:
            frequencies.extend(fixed_fields(lines[cursor]))
            cursor += 1
        if len(frequencies) != 3 * nat or not all(map(math.isfinite, coordinates + frequencies)):
            raise ValueError("invalid native frequency record")
        rows.append({"grid_index": list(index),
                     "q_fractional": [i / n for i, n in zip(index, grid)],
                     "q_cartesian_2pi_alat": coordinates[:3],
                     "frequency_cm_1": frequencies})
    if any(line.strip() for line in lines[cursor:]):
        raise ValueError("unexpected trailing matdyn frequency records")
    return rows


def run_native(name, directory, input_text):
    executable, _ = BINARIES[name]
    input_path = directory / f"{name}.in"
    input_path.write_text(input_text)
    command = {"argv": [executable], "stdin": str(input_path), "cwd": str(directory)}
    save_receipt(directory / f"{name}-command.json", command)
    with input_path.open("r") as source, (directory / f"{name}.stdout").open("w") as out, \
            (directory / f"{name}.stderr").open("w") as err:
        subprocess.run([executable], cwd=directory, stdin=source, stdout=out, stderr=err, check=True)
    if "JOB DONE." not in (directory / f"{name}.stdout").read_text():
        raise ValueError(f"native {name} did not publish normal completion")
    print(f"NATIVE_SERIALIZED_OPERATOR_STAGE_COMPLETE {name} {directory}", flush=True)


def replay_serialized_operators(result_path, native_directory):
    """Return a deeply bound native replay, preserving interrupted attempts."""
    result_path, native_directory = Path(result_path).resolve(), Path(native_directory).resolve()
    envelope = json.loads(result_path.read_text())
    if envelope["status"] != "complete":
        raise ValueError("native PH stage is not complete")
    data = envelope["data"]
    grid = data["q_manifest"]["grid"]
    nat = data["dynamical_matrices"][0]["geometry"]["number_of_atoms"]
    inputs = {str(path): sha256(path) for path in (
        result_path, native_directory / "dyn0", Path(__file__), HERE / "oracle_artifacts.py",
        *(native_directory / Path(record["artifact"]).name for record in data["dynamical_matrices"]),
    )}
    for name, (executable, expected) in BINARIES.items():
        digest = sha256(executable)
        if digest != expected:
            raise ValueError(f"{name} binary differs from the registered QE7.5 build")
        inputs[executable] = digest
    import hashlib
    fingerprint = hashlib.sha256(json.dumps(inputs, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    root = HERE / "native-replay" / fingerprint
    root.mkdir(parents=True, exist_ok=True)
    receipt_path = root / "receipt.json"
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt["input_sha256"] != inputs or receipt["status"] != "complete":
            raise ValueError("native replay receipt is not bound to these complete inputs")
        directory = Path(receipt["frequency_file"]).parent
        required = {"q2r.in", "q2r.stdout", "q2r.stderr", "q2r-command.json",
                    "matdyn.in", "matdyn.stdout", "matdyn.stderr", "matdyn-command.json",
                    "ifc.xml", "sampled.freq"}
        actual = {str(path) for path in directory.rglob("*") if path.is_file()}
        if (directory.parent != root or Path(receipt["frequency_file"]).name != "sampled.freq"
                or not required <= {Path(path).name for path in actual}
                or actual != set(receipt["artifact_sha256"])
                or receipt["grid"] != grid or receipt["number_of_atoms"] != nat):
            raise ValueError("native replay cache has incomplete or displaced artifacts/metadata")
        for path, digest in receipt["artifact_sha256"].items():
            if sha256(path) != digest:
                raise ValueError(f"native replay artifact changed: {path}")
        rows = frequency_table(receipt["frequency_file"], nat, grid)
        if rows != receipt["rows"]:
            raise ValueError("native replay metadata differs from the retained raw table")
        return receipt
    directory = root / uuid.uuid4().hex
    directory.mkdir()
    # has_xml strips the explicit .xml marker; do_q2r then reads dyn0 and dynN.xml.
    fildyn = native_directory / "dyn.xml"
    run_native("q2r", directory,
               f"&input\n fildyn='{fildyn}', flfrc='ifc', zasr='no', write_lr=.false.\n/\n")
    indices = list(itertools.product(*(range(n) for n in grid)))
    points = "\n".join(" ".join(format(i / n, ".17g") for i, n in zip(index, grid))
                       for index in indices)
    run_native("matdyn", directory,
               "&input\n flfrc='ifc.xml', flfrq='sampled.freq', flvec='sampled.modes',\n"
               " asr='no', q_in_cryst_coord=.true., q_in_band_form=.false.,\n"
               " loto_disable=.true., na_ifc=.false., read_lr=.false., fd=.false.\n/\n"
               f"{len(indices)}\n{points}\n")
    frequency_path = directory / "sampled.freq"
    receipt = {"schema_version": 1, "status": "complete", "input_sha256": inputs,
               "frequency_file": str(frequency_path), "grid": grid, "number_of_atoms": nat,
               "rows": frequency_table(frequency_path, nat, grid),
               "native_frequency_format": "matdyn.f90:779, (6f10.4), inverse centimetres",
               "frequency_rounding_half_width_cm_1": 0.00005,
               "q_coordinate_rounding_half_width": 0.0000005,
               "processing": "q2r explicit Hermitian average; no ASR, no nonanalytic/polar terms; original sampled q only",
               "artifact_sha256": {str(path): sha256(path) for path in directory.rglob("*") if path.is_file()}}
    save_receipt(receipt_path, receipt)
    return receipt
