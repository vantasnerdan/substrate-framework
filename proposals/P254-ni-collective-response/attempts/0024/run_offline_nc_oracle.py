#!/usr/bin/env python3
"""New UPF->beta/D/local-field comparator; never starts QE or prior oracles.

Run only after Root source closure. All comparisons use existing actual0018
Gamma-q exports over every stored k/spin and both source/target PW maps.
The independent elementary checks are convention/quadrature checks, not
material accuracy claims. Fixed-PW displaced ionic-action checks belong to Root.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.special import lpmv, spherical_jn

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from substrate_framework.qe_nc_upf import read_nc_upf, expand_channels
from substrate_framework.qe_nc_radial import build_radial_tables, qe_simpson_weights, qe_spherical_bessel
from substrate_framework.qe_nc_projectors import centered_projectors, place_projectors, qe_real_harmonics
from substrate_framework.qe_response_bundle import load_bundle

BASE = Path("/home/administrator/projects/catalyst-nc-operator/runs/2026-10-19_collective-response")
DEFAULT_BUNDLE = BASE / "full-ionic-operators-0018/gamma/bundle/q-1"
DEFAULT_MATERIAL = BASE / "nc-material/gamma-continuation-002/loaded-ni31h6/native/scratch/scf.save"
RTOL, ATOL = 2e-10, 2e-11


def compare(report, name, actual, reference, *, rtol=RTOL, atol=ATOL):
    actual, reference = np.asarray(actual), np.asarray(reference)
    difference = actual-reference
    max_abs = float(np.max(np.abs(difference), initial=0))
    scale = float(np.max(np.abs(reference), initial=0))
    limit = atol + rtol*scale
    record = {"name": name, "max_abs": max_abs, "reference_max_abs": scale,
              "relative_l2": float(np.linalg.norm(difference)/max(np.linalg.norm(reference), 1e-300)),
              "limit": limit, "passed": bool(np.isfinite(max_abs) and max_abs <= limit)}
    report["checks"].append(record)
    return record


def mutation(report, name, wrong, reference):
    wrong, reference = np.asarray(wrong), np.asarray(reference)
    error = float(np.max(np.abs(wrong-reference), initial=0))
    limit = ATOL + RTOL*float(np.max(np.abs(reference), initial=0))
    report["mutations"].append({"name": name, "max_abs": error,
                                "comparison_limit": limit, "exposed": bool(error > 100*limit)})


def records(value):
    if isinstance(value, dict):
        if {"path", "bytes", "sha256"} <= value.keys():
            yield value
        for child in value.values():
            yield from records(child)
    elif isinstance(value, list):
        for child in value:
            yield from records(child)


def bound_bytes(path, bindings, report):
    path = Path(path).resolve()
    binding = bindings[str(path)]
    data = path.read_bytes()
    if len(data) != int(binding["bytes"]) or hashlib.sha256(data).hexdigest() != binding["sha256"]:
        raise ValueError(f"original source binding mismatch: {path}")
    report["inputs"].append({"path": str(path), "bytes": str(len(data)), "sha256": binding["sha256"]})
    return data


def elementary_checks(report):
    # Polynomial moments are independent analytic integrals. Even QE Simpson
    # is quadratic-exact, not cubic-exact; no false cubic claim is made.
    for n, degree in ((9, 3), (8, 2)):
        x = np.linspace(0, 1, n)
        w = qe_simpson_weights(np.full(n, 1/(n-1)))
        for power in range(degree+1):
            compare(report, f"simpson_n{n}_moment{power}", w @ x**power,
                    1/(power+1), rtol=0, atol=3e-15)
        mutation(report, f"simpson_n{n}_omitted_endpoint", w[:-1] @ x[:-1]**2, 1/3)
    x = np.array([0, 1e-8, .03, .0499, .0501, .13, 1.1, 4.2, 15.0])
    for l in range(3):
        compare(report, f"bessel_l{l}_independent_scipy", qe_spherical_bessel(l, x),
                spherical_jn(l, x), rtol=0, atol=2e-12)
    p = np.array([[.31, -.72, .43], [-.5, .3, -.7], [0, 0, 1], [1, 0, 0], [0, -1, 0]])
    radius = np.linalg.norm(p, axis=1)
    cost, phi = p[:, 2]/radius, np.arctan2(p[:, 1], p[:, 0])
    independent = []
    for l in range(3):
        for m in range(l+1):
            radial = np.sqrt((2*l+1)/(4*np.pi)*math.factorial(l-m)/math.factorial(l+m))*lpmv(m, l, cost)
            independent.append(radial if m == 0 else np.sqrt(2)*radial*np.cos(m*phi))
            if m:
                independent.append(np.sqrt(2)*radial*np.sin(m*phi))
    y = qe_real_harmonics(p, 2)
    compare(report, "harmonics_independent_associated_legendre", y, np.column_stack(independent), atol=2e-14, rtol=0)
    for l in range(3):
        compare(report, f"harmonics_addition_theorem_l{l}", np.sum(y[:, l*l:(l+1)**2]**2, axis=1),
                np.full(len(p), (2*l+1)/(4*np.pi)), atol=2e-14, rtol=0)
    # An off-diagonal radial D forces the lm-selection rule to be tested even
    # though the actual selected UPFs happen to have diagonal raw PP_DIJ.
    ell = np.array([1, 1, 0])
    raw = np.array([[2., .7, .3], [.7, 3., -.2], [.3, -.2, 5.]])
    radial, _, _, d = expand_channels(ell, raw)
    expected = np.zeros((7, 7))
    expected[:3, :3], expected[3:6, 3:6] = 2*np.eye(3), 3*np.eye(3)
    expected[:3, 3:6], expected[3:6, :3], expected[6, 6] = .7*np.eye(3), .7*np.eye(3), 5
    compare(report, "D_independent_channel_selection", d, expected, rtol=0, atol=0)
    mutation(report, "D_without_angular_selection", raw[np.ix_(radial, radial)], expected)


def quadrature_checks(report, table):
    upf, omega = table.upf, table.omega_bohr3
    n = upf.beta_mesh
    # Independent continuous integration of a cubic interpolant through the
    # original radial samples; this tests finite-grid Simpson accuracy rather
    # than demanding it equal continuum quadrature to machine precision.
    for nb, l in enumerate(upf.angular_l):
        spline = CubicSpline(upf.r[:n], upf.beta_r[:n, nb])
        q = .37
        integral, error = quad(lambda r: float(spline(r))*r*spherical_jn(int(l), q*r),
                               upf.r[0], upf.r[n-1], epsabs=1e-10, epsrel=1e-10, limit=150)
        direct = integral*4*np.pi/np.sqrt(omega)
        compare(report, f"{upf.element}_beta{nb}_adaptive_radial_quadrature", table.beta(np.array([q]))[0, nb],
                direct, rtol=2e-6, atol=2e-7)
        report["quadrature_error_estimates"].append({"element": upf.element, "radial_projector": nb,
                                                    "quad_absolute_estimate": error})
    local_g0_quadrature_check(report, table)


def local_g0_quadrature_check(report, table):
    """Integrate each known spline interval; do not exhaust adaptive discovery."""
    upf, omega = table.upf, table.omega_bohr3
    n = upf.local_mesh
    spline = CubicSpline(upf.r[:n], upf.r[:n]*(upf.r[:n]*upf.local_ry[:n]+2*upf.z_valence))
    direct, error = quad(lambda r: float(spline(r)), upf.r[0], upf.r[n-1],
                         epsabs=1e-9, points=upf.r[1:n-1], limit=n+1)
    compare(report, f"{upf.element}_local_G0_adaptive_quadrature", table.local_g0,
            4*np.pi/omega*direct, rtol=2e-6, atol=2e-7)
    report["quadrature_error_estimates"].append({
        "element": upf.element, "component": "local_G0",
        "quad_absolute_estimate": error, "known_interior_knots": n-2})
    mutation(report, f"{upf.element}_wrong_G0_short_range_limit",
             table.local_short_values[0], table.local_g0)


def load_inputs(args, report):
    inventory = json.loads(args.inventory.read_bytes())
    bindings = {str(Path(item["path"]).resolve()): item for item in records(inventory)
                if Path(item["path"]).is_absolute()}
    saved = args.bundle.parent.parent / "scratch/scf.save"
    xml = ET.fromstring(bound_bytes(saved / "data-file-schema.xml", bindings, report))
    if xml.attrib.get("Units") != "Hartree atomic units":
        raise ValueError("actual saved XML must retain Hartree atomic units")
    atoms = xml.findall("output/atomic_structure/atomic_positions/atom")
    tau = np.array([[float(x) for x in atom.text.split()] for atom in atoms])
    species = xml.findall("output/atomic_species/species")
    upfs = []
    for sp in species:
        name = sp.findtext("pseudo_file").strip()
        material_bytes = bound_bytes(args.material_saved / name, bindings, report)
        fixture_bytes = bound_bytes(saved / name, bindings, report)
        if material_bytes != fixture_bytes:
            raise ValueError("actual material/comparator UPF originals differ")
        record = bindings[str((saved/name).resolve())]
        upfs.append(read_nc_upf(saved/name, expected_sha256=record["sha256"], expected_bytes=record["bytes"]))
    bound_bytes(args.bundle / "manifest.json", bindings, report)
    bundle = load_bundle(args.bundle)
    bg = bundle.background
    if not bg.metadata["lgamma"] or np.any(bg.array("q_cart_2pi_alat") != 0):
        raise ValueError("this oracle is the existing Gamma-q comparator only")
    species_names = [sp.attrib["name"] for sp in species]
    expected_types = [species_names.index(atom.attrib["name"])+1 for atom in atoms]
    compare(report, "XML_atom_species_order", bg.array("ityp"), expected_types, rtol=0, atol=0)
    cell = np.column_stack([[float(x) for x in xml.findtext(f"output/atomic_structure/cell/a{i}").split()]
                            for i in (1, 2, 3)])
    compare(report, "XML_volume", abs(np.linalg.det(cell)), bg.metadata["omega"], rtol=1e-13, atol=0)
    compare(report, "XML_lattice", cell, bg.metadata["alat"]*bg.array("at"), rtol=1e-13, atol=0)
    return bundle, tau, upfs


def compare_projectors(report, bundle, tau, tables):
    bg = bundle.background
    types = bg.array("ityp")
    g = bg.metadata["tpiba"]*bg.array("g").T
    cache = {}
    for state in bundle.states:
        spin = state.metadata["spin"]-1
        for side in ("source", "target"):
            indices = state.array(f"igk_{side}").astype(np.int64)-1
            p = g[indices] + bg.metadata["tpiba"]*state.array(f"k_{side}")
            native = state.array(f"vkb_{side}")[:len(p)]
            offset = 0
            for species, table in enumerate(tables, start=1):
                key = (species, p.tobytes())
                if key not in cache:
                    cache[key] = centered_projectors(table, p)
                centered = cache[key]
                for atom in np.flatnonzero(types == species):
                    nh = table.upf.nh
                    b = place_projectors(centered, p, tau[atom])
                    ref = native[:, offset:offset+nh]
                    label = f"{state.path}:{side}:atom{atom+1}"
                    compare(report, label+":beta", b, ref)
                    compare(report, label+":active_D", table.upf.expanded_d,
                            state.array("deeq")[:nh, :nh, atom, spin], rtol=0, atol=0)
                    if state is bundle.states[0] and side == "source":
                        mutation(report, label+":wrong_beta_normalization", b/(4*np.pi), ref)
                        if np.any(table.upf.channel_l % 2):
                            mutation(report, label+":wrong_angular_sign", b*((-1.)**table.upf.channel_l), ref)
                        if np.linalg.norm(tau[atom]) > 0:
                            mutation(report, label+":wrong_atomic_phase", place_projectors(centered, p, -tau[atom]), ref)
                        mutation(report, label+":D_Ry_to_Hartree_error", table.upf.expanded_d/2,
                                 state.array("deeq")[:nh, :nh, atom, spin])
                    offset += nh
            compare(report, f"{state.path}:{side}:channel_count", offset, state.metadata["nkb"], rtol=0, atol=0)
    report["centered_cache_entries"] = len(cache)


def compare_local_fields(report, bundle, tau, tables):
    bg, meta = bundle.background, bundle.background.metadata
    shape = tuple(meta[f"nr{i}"] for i in (1, 2, 3))
    storage = tuple(meta[f"nr{i}x"] for i in (1, 2, 3))
    if shape != storage:
        raise ValueError("actual0018 comparator requires unpadded serial FFT storage")
    mill = bg.array("mill")
    slots = np.ravel_multi_index(tuple(mill[i] % shape[i] for i in range(3)), shape, order="F")
    compare(report, "actual_mill_to_nl", slots+1, bg.array("nl"), rtol=0, atol=0)
    p = meta["tpiba"]*(bg.array("bg") @ mill).T
    compare(report, "actual_mill_to_cartesian_G", p, meta["tpiba"]*bg.array("g").T)
    q = np.linalg.norm(p, axis=1)
    vloc = [table.vloc(q) for table in tables]
    def field(coefficients):
        spectrum = np.zeros(np.prod(shape), dtype=np.complex128)
        spectrum[slots] = coefficients
        return np.fft.ifftn(spectrum.reshape(shape, order="F"), norm="forward").reshape(-1, order="F")
    for atom, species in enumerate(bg.array("ityp")):
        v = vloc[species-1]
        phase = np.exp(-1j*(p @ tau[atom]))
        table = tables[species-1]
        tail = np.zeros_like(q)
        nz = q != 0
        tail[nz] = 8*np.pi*table.upf.z_valence/meta["omega"]*np.exp(-q[nz]**2/4)/q[nz]**2
        for axis in range(3):
            coefficients = -1j*p[:, axis]*phase*v
            ref = bg.array("vloc_bare_cartesian")[:, 3*atom+axis]
            result = field(coefficients)
            label = f"local_atom{atom+1}_axis{axis}"
            compare(report, label, result, ref)
            mutation(report, label+":wrong_derivative_sign", -result, ref)
            mutation(report, label+":omitted_Coulomb_tail", field(-1j*p[:, axis]*phase*(v+tail)), ref)
            mutation(report, label+":wrong_FFT_normalization", result/np.prod(shape), ref)


def run(args, report):
    elementary_checks(report)
    bundle, tau, upfs = load_inputs(args, report)
    bg = bundle.background
    g = bg.metadata["tpiba"]*bg.array("g").T
    qmax = float(np.linalg.norm(g, axis=1).max())
    for state in bundle.states:
        for side in ("source", "target"):
            p = g[state.array(f"igk_{side}")-1] + bg.metadata["tpiba"]*state.array(f"k_{side}")
            qmax = max(qmax, float(np.linalg.norm(p, axis=1).max()))
    tables = [build_radial_tables(upf, bg.metadata["omega"], qmax) for upf in upfs]
    report["table_build_count"] = len(tables)
    report["species"] = [{"element": t.upf.element, "nbeta": t.upf.nbeta, "nh": t.upf.nh,
                          "local_msh": t.upf.local_mesh, "beta_kkbeta": t.upf.beta_mesh,
                          "nqx": len(t.beta_values), "qmax_bohr_inv": qmax} for t in tables]
    for table in tables:
        quadrature_checks(report, table)
    compare_projectors(report, bundle, tau, tables)
    compare_local_fields(report, bundle, tau, tables)


def json_safe(value):
    """Keep nonfinite failure evidence writable under strict JSON."""
    if isinstance(value, float) and not math.isfinite(value):
        return {"nonfinite": str(value)}
    if isinstance(value, dict):
        return {key: json_safe(child) for key, child in value.items()}
    if isinstance(value, list):
        return [json_safe(child) for child in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--material-saved", type=Path, default=DEFAULT_MATERIAL)
    parser.add_argument("--inventory", type=Path, default=Path(__file__).with_name("source-access-inventory.json"))
    args = parser.parse_args()
    report = {"kind": "P254-0024-new-offline-NC-reconstruction-oracle", "checks": [], "mutations": [],
              "inputs": [], "quadrature_error_estimates": [], "rtol": RTOL, "atol": ATOL,
              "scope": "Actual0018 Gamma-q beta/D/local reconstruction, all existing k/spin states; no Native execution, no material action or derivative admission."}
    start = time.monotonic()
    try:
        run(args, report)
    except Exception as exc:
        report["exception"] = {"type": type(exc).__name__, "message": str(exc)}
    report["elapsed_seconds"] = time.monotonic()-start
    report["passed"] = ("exception" not in report and bool(report["checks"]) and
                        all(c["passed"] for c in report["checks"]) and bool(report["mutations"]) and
                        all(m["exposed"] for m in report["mutations"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(json_safe(report), indent=2, allow_nan=False)+"\n", encoding="utf-8")
    if not report["passed"]:
        raise SystemExit(f"offline NC oracle failed; metrics saved to {args.output}")
    print(json.dumps({"passed": True, "output": str(args.output), "elapsed_seconds": report["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
