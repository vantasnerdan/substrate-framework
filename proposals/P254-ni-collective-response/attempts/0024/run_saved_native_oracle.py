"""New saved-reader oracle against original 0018 Native exports and material XML.

No H0 solve, SCF/PH execution, original oracle replay, or normalization repair.
Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 and --output PATH. Every
check (including failures) is emitted before nonzero exit. Tolerances below
are absolute serialization/normalization checks, not material calibration.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import struct
import tempfile
import time
import xml.etree.ElementTree as ET

import numpy as np

from substrate_framework.qe_saved_native import load_density, load_wavefunctions, read_bound_bytes
from substrate_framework.qe_saved_reference import load_saved_reference

HERE = Path(__file__).resolve().parent
INVENTORY = dict(path=str(HERE / "source-access-inventory.json"), bytes="210197",
                 sha256="daa1d1509c50e01bff328f8ae215ad0bf4a71815da0aa1ab06295776f36824f1")
CHECKS = []
SOURCES = {}


def check(name, passed, **metrics):
    CHECKS.append(dict(name=name, passed=bool(passed), **metrics))


def error(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


def close(name, a, b, tolerance):
    delta = error(a, b)
    check(name, delta <= tolerance, max_abs_error=delta, tolerance=tolerance)


def bound_read(binding):
    data, receipt = read_bound_bytes(binding)
    SOURCES[receipt["path"]] = dict(receipt)
    return data


def source(binding, loader):
    result = loader(binding)
    SOURCES[result.source["path"]] = dict(result.source)
    return result


def array(document, name, base, bindings):
    descriptor = document["datasets"][name]
    path = str(base / descriptor["path"])
    binding = bindings[path]
    if str(descriptor["bytes"]) != binding["bytes"] or descriptor["sha256"] != binding["sha256"]:
        raise ValueError("independent exported descriptor differs from original inventory")
    dtype = {"float64": "<f8", "complex128": "<c16", "int32": "<i4"}[descriptor["dtype"]]
    if descriptor["order"] != "F":
        raise ValueError("unsupported comparator array order")
    return np.frombuffer(bound_read(binding), dtype=dtype).reshape(tuple(int(n) for n in descriptor["shape"]), order="F")


def miller_order(saved, target):
    rows = {tuple(int(x) for x in row): i for i, row in enumerate(saved.T)}
    if len(rows) != saved.shape[1]:
        raise ValueError("duplicate saved Miller rows")
    return np.array([rows[tuple(int(x) for x in row)] for row in target.T], dtype=np.int64)


def consistent(wfc, ref, expected_ik, expected_spin):
    point, energy, _ = ref.select(expected_ik, expected_spin)
    return (wfc.ik == expected_ik and wfc.ispin == expected_spin
            and wfc.igwx == point.npw and wfc.nbnd == energy.size
            and wfc.gamma_only == ref.gamma_only and wfc.npol == 1 and wfc.scalef == 1
            and error(wfc.k_cart_bohr_inv, point.k_cart_bohr_inv) < 2e-13
            and error(wfc.reciprocal_bohr_inv, ref.reciprocal_bohr_inv) < 2e-13)


def fixture_checks(inv, bindings):
    base = Path(inv["existing_comparator_manifest"]["path"]).parent
    saved = base.parent.parent / "scratch/scf.save"
    bg = json.loads(bound_read(bindings[str(base / "background.json")]))
    ref = source(bindings[str(saved / "data-file-schema.xml")], load_saved_reference)
    mill = array(bg, "mill", base, bindings)
    close("fixture physical reciprocal columns", ref.reciprocal_bohr_inv,
          array(bg, "bg", base, bindings) * bg["tpiba"], 2e-13)
    close("fixture geometry columns", ref.cell_bohr,
          array(bg, "at", base, bindings) * bg["alat"], 2e-13)
    unequal = 0
    witness = None
    for spin in (1, 2):
        for ik in range(1, ref.nks + 1):
            global_ik = (spin - 1) * ref.nks + ik
            metadata = json.loads(bound_read(bindings[str(base / f"state-{spin}-{global_ik}.json")]))
            filename = f"wfc{'up' if spin == 1 else 'dw'}{ik}.dat"
            binding = bindings[str(saved / filename)]
            wfc = source(binding, load_wavefunctions)
            label = f"fixture spin={spin} ik={ik}"
            check(label + " local-index/header/XML", consistent(wfc, ref, ik, spin),
                  ik=wfc.ik, ispin=wfc.ispin, ngw=wfc.ngw, igwx=wfc.igwx, nbnd=wfc.nbnd)
            unequal += int(wfc.ngw != wfc.igwx)
            igk = array(metadata, "igk_source", base, bindings).astype(np.int64) - 1
            target_mill = mill[:, igk]
            order = miller_order(wfc.miller, target_mill)
            exported = array(metadata, "evc", base, bindings)[:int(metadata["ngw_source"])]
            close(label + " original Native evc Miller-reordered", wfc.coefficients[order], exported, 2e-12)
            close(label + " physical k units", wfc.k_cart_bohr_inv,
                  array(metadata, "k_source", base, bindings) * bg["tpiba"], 2e-13)
            point, energies, raw = ref.select(ik, spin)
            close(label + " Ha to Ry eigenvalues", energies,
                  array(metadata, "eigenvalues_source", base, bindings), 5e-11)
            weight = array(metadata, "source_kweights", base, bindings)[0]
            close(label + " unchanged raw occupation convention", raw * weight,
                  array(metadata, "weighted_native_wg", base, bindings), 5e-12)
            if spin == 2 and ik == 1:
                check("global spin-stacked index mutation rejected", not consistent(replace(wfc, ik=global_ik), ref, ik, spin))
                check("implicit spin mutation rejected", not consistent(replace(wfc, ispin=1), ref, ik, spin))
            if spin == 1 and ik == 2:
                witness = (binding, wfc, ref, exported, order)
                check("extra tpiba k scaling mutation rejected", not consistent(replace(wfc, k_cart_bohr_inv=wfc.k_cart_bohr_inv * bg["tpiba"]), ref, ik, spin))
                check("unconverted XML energy mutation rejected", error(point.eigenvalues_ha[spin - 1], energies) > 5e-11)
    check("general k distinguishes global ngw from compact igwx", unequal > 0, unequal_header_count=unequal)
    density = source(bindings[str(saved / "charge-density.dat")], load_density)
    order = miller_order(density.miller, mill)
    real_density = array(bg, "rho_total_magnetization", base, bindings)
    shape = tuple(int(bg[f"nr{i}"]) for i in (1, 2, 3))
    nodes = tuple(np.mod(mill[i], shape[i]) for i in range(3))
    for component in range(density.nspin):
        independent = np.fft.fftn(real_density[:, component].reshape(shape, order="F")) / np.prod(shape)
        close(f"density original total/magnetization component {component}", density.components[order, component], independent[nodes], 2e-11)
    close("density original physical reciprocal basis", density.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, 2e-13)
    return witness


def fft_checks(wfc, ref, label):
    shape = ref.fft_smooth
    indices = tuple(np.mod(wfc.miller[i], shape[i]) for i in range(3))
    linear = np.ravel_multi_index(indices, shape)
    check(label + " injective physical FFT embedding", len(np.unique(linear)) == wfc.igwx)
    gram = wfc.coefficients.conj().T @ wfc.coefficients
    close(label + " all-band overlap", gram, np.eye(wfc.nbnd), 2e-8)
    selected = sorted({0, wfc.nbnd // 2, wfc.nbnd - 1})
    volume = abs(float(np.linalg.det(ref.cell_bohr)))
    for band in selected:
        grid = np.zeros(shape, dtype=np.complex128)
        grid[indices] = wfc.coefficients[:, band]
        psi = np.fft.ifftn(grid) * np.prod(shape) / np.sqrt(volume)
        physical_norm = float(volume * np.mean(np.abs(psi) ** 2))
        check(label + f" physical FFT norm band {band + 1}", abs(physical_norm - 1) <= 2e-8,
              integral_abs_psi_squared=physical_norm, tolerance=2e-8)
        roundtrip = np.fft.fftn(psi) * np.sqrt(volume) / np.prod(shape)
        close(label + f" FFT roundtrip band {band + 1}", roundtrip[indices], wfc.coefficients[:, band], 2e-13)
        if band == 0:
            wrong_norm = physical_norm / float(np.prod(shape)) ** 2
            check(label + " omitted inverse FFT factor mutation rejected", abs(wrong_norm - 1) > 2e-8)


def material_checks(inv):
    inputs = {Path(b["path"]).name: b for b in inv["material_inputs"]}
    xml_binding = inputs["data-file-schema.xml"]
    ref = source(xml_binding, load_saved_reference)
    original = ET.fromstring(bound_read(xml_binding)).find("output")
    check("material original full atom identity", ref.nat == 37 and ref.atom_names == ("Ni",) * 31 + ("H",) * 6,
          nat=ref.nat, atom_indices=list(ref.atom_indices), species_order=[s.name for s in ref.species])
    close("material XML cell reciprocal duality", ref.cell_bohr.T @ ref.reciprocal_bohr_inv, 2 * np.pi * np.eye(3), 2e-12)
    check("material physical input sizes", ref.fft_grid == (72, 72, 72) and ref.fft_smooth == (72, 72, 72)
          and ref.nbnd_by_spin == (346, 346) and ref.ngm == 148597,
          fft_grid=list(ref.fft_grid), nbnd_by_spin=list(ref.nbnd_by_spin), ngm=ref.ngm)
    original_atoms = original.findall("atomic_structure/atomic_positions/atom")
    positions = np.array([[float(x) for x in a.text.split()] for a in original_atoms]).T
    close("material all original atom positions", ref.positions_bohr, positions, 0.0)
    for parsed, raw in zip(ref.kpoints, original.findall("band_structure/ks_energies")):
        original_occ = np.array([float(x) for x in raw.findtext("occupations").split()])
        close(f"material original raw MV ik={parsed.ik}", np.concatenate(parsed.raw_occupations), original_occ, 0.0)
    occ = np.concatenate(ref.kpoints[0].raw_occupations)
    check("material raw MV is not a probability preparation", ref.smearing == "mv" and np.max(occ) > 1,
          raw_MV_min=float(np.min(occ)), raw_MV_max=float(np.max(occ)))
    check("MV probability clipping mutation detected", error(occ, np.clip(occ, 0, 1)) > 0)
    density = source(inputs["charge-density.dat"], load_density)
    check("material density dimensions", density.ngm == ref.ngm and density.nspin == len(ref.nbnd_by_spin))
    close("material density physical reciprocal columns", density.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, 2e-13)
    zero = np.flatnonzero(np.all(density.miller == 0, axis=0))
    check("material unique density zero mode", zero.size == 1)
    volume = abs(float(np.linalg.det(ref.cell_bohr)))
    if zero.size == 1:
        close("material density total charge normalization", density.components[zero[0], 0] * volume, ref.nelec, 2e-8)
        close("material density magnetization normalization", density.components[zero[0], 1] * volume,
              float(original.findtext("magnetization/total")), 2e-8)
    for spin, filename in ((1, "wfcup1.dat"), (2, "wfcdw1.dat")):
        wfc = source(inputs[filename], load_wavefunctions)
        label = f"material spin={spin}"
        check(label + " selected exact Gamma/header/XML", consistent(wfc, ref, 1, spin) and np.all(wfc.k_cart_bohr_inv == 0)
              and wfc.igwx == 18613, ik=wfc.ik, ispin=wfc.ispin, ngw=wfc.ngw, igwx=wfc.igwx, nbnd=wfc.nbnd)
        order = miller_order(density.miller, wfc.miller)
        check(label + " every PW maps into density Miller basis", len(np.unique(order)) == wfc.igwx)
        check(label + " redundant tpiba basis mutation detected", not consistent(replace(wfc, reciprocal_bohr_inv=wfc.reciprocal_bohr_inv * (2 * np.pi / ref.alat_bohr)), ref, 1, spin))
        fft_checks(wfc, ref, label)


def rejection(name, function):
    try:
        function()
    except (ValueError, ET.ParseError) as exc:
        check(name, True, exception=type(exc).__name__, message=str(exc))
    else:
        check(name, False, message="mutant accepted")


def mutation_checks(witness):
    binding, wfc, ref, exported, order = witness
    rejection("independent source digest mutation rejected", lambda: load_wavefunctions(dict(binding, sha256="0" * 64)))
    data = bound_read(binding)
    check("Miller row offset mutation detected", error(wfc.coefficients[np.roll(order, 1)], exported) > 2e-12)
    check("band phase repair mutation detected", error(-wfc.coefficients[order], exported) > 2e-12)
    with tempfile.TemporaryDirectory(prefix="p254-saved-reader-oracle-") as directory:
        path = Path(directory) / "mutated-wfc.dat"
        def parse_mutated(payload):
            path.write_bytes(payload)
            return load_wavefunctions(dict(path=str(path), bytes=str(len(payload)), sha256=hashlib.sha256(payload).hexdigest()))
        rejection("independently rebound truncated record rejected", lambda: parse_mutated(data[:-7]))
        for name, offset, fmt, value in (
            ("malformed trailing marker", 48, "<i", 40),
            ("big-endian marker ABI", 0, ">i", 44),
            ("gamma compressed", 36, "<i", 1),
            ("spinor header", 64, "<i", 2),
            ("scale factor", 40, "<d", 2.0),
        ):
            payload = bytearray(data)
            struct.pack_into(fmt, payload, offset, value)
            rejection(name + " mutation rejected", lambda p=bytes(payload): parse_mutated(p))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    fatal = None
    try:
        inv = json.loads(bound_read(INVENTORY))
        bindings = {b["path"]: b for key in ("existing_comparator_all_files", "existing_comparator_saved_input_sources", "material_inputs") for b in inv[key]}
        # Verify consulted original primary records before any dependent numerical use.
        for binding in inv["primary_QE_sources"]:
            if Path(binding["path"]).name in ("io_base.f90", "pw_restart_new.f90"):
                bound_read(binding)
        for binding in (
            dict(path="/opt/qebuild/q-e-qe-7.5/PW/src/io_rho_xml.f90", bytes="12640", sha256="400488ecf2e48d4ba04d2d7692ec0cfd5ca1d96003054ae40761ea2beb6f7a9d"),
            dict(path="/opt/qebuild/q-e-qe-7.5/Modules/qexsd_init.f90", bytes="74976", sha256="aa4d0cab034ec79d39a55afa1ae5f6b78fb08c5864789a05e16eae69dd59231d"),
        ):
            bound_read(binding)
        witness = fixture_checks(inv, bindings)
        mutation_checks(witness)
        material_checks(inv)
    except Exception as exc:
        fatal = dict(type=type(exc).__name__, message=str(exc))
    passed = fatal is None and bool(CHECKS) and all(c["passed"] for c in CHECKS)
    report = dict(kind="P254-0024-new-saved-reader-oracle", passed=passed, checks=CHECKS,
                  fatal=fatal, sources=list(SOURCES.values()), elapsed_seconds=time.monotonic() - started,
                  scope="Independent actual0018 full-band/native-PW comparisons, actual selected material records/XML and full-band overlap; FFT norms at first/middle/last bands. No material H0 solve or quantum preparation.")
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text)
    print(text, end="", flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
