#!/usr/bin/env python3
"""New saved-density reconstruction versus immutable original0018 exports.

Run only after Root source admission/build. No Native/SCF/PH invocation and no
old oracle replay. Raw fields/differences and prospective allowances are saved;
retained eigenpairs and finite literal H0 algebra are qualified diagnostics,
not Native unretained-complement action evidence or a convergence gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace
import zipfile

import numpy as np
# Preload lazy consumers before capturing the execution-source inventory.
import numpy.fft
import numpy.linalg
import scipy.special
from substrate_framework import qe_response_schema
from substrate_framework.qe_internal_pbe import NativePBE
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_material_reference import (
    COEFFICIENT_NAMES, load_material_inputs, reconstruct_reference,
)
from substrate_framework.qe_material_hamiltonian import (
    MaterialHamiltonian, PWFactors, SpeciesProjectors, build_pw_factors,
)
from substrate_framework.qe_response_bundle import load_bundle
from substrate_framework.qe_response_fft import NativeRhoFFT
from substrate_framework.qe_saved_native import load_wavefunctions, read_bound_bytes
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.spin_pbe_kernel import apply_spin_pbe_arrays, apply_spin_pbe_kernel

HERE = Path(__file__).resolve().parent
INVENTORY = dict(path=str(HERE.parent / "0024/source-access-inventory.json"), bytes="210197",
                 sha256="daa1d1509c50e01bff328f8ae215ad0bf4a71815da0aa1ab06295776f36824f1")


def binding(path, data=None):
    path = Path(path).resolve()
    data = path.read_bytes() if data is None else data
    return dict(path=str(path), bytes=str(len(data)), sha256=hashlib.sha256(data).hexdigest())


def records(value):
    if isinstance(value, dict):
        if {"path", "bytes", "sha256"} <= value.keys():
            yield value
        for child in value.values():
            yield from records(child)
    elif isinstance(value, list):
        for child in value:
            yield from records(child)


class Report:
    def __init__(self):
        self.checks = []
        self.fields = {}
        self.sources = {}
        self.sensitivity = []
        self.diagnostics = []

    def bound(self, record):
        data, receipt = read_bound_bytes(record)
        self.sources[receipt["path"]] = dict(receipt)
        return data

    def compare(self, name, actual, native, atol, rtol, *, gate=True):
        actual, native = np.asarray(actual), np.asarray(native)
        if actual.shape != native.shape:
            raise ValueError(f"{name}: shape mismatch {actual.shape} != {native.shape}")
        delta = actual - native
        maximum = float(np.max(np.abs(delta)))
        limit = float(atol + rtol*np.max(np.abs(native)))
        key = name.replace(":", "_").replace("/", "_")
        self.fields.update({key+"_actual": actual, key+"_native": native, key+"_difference": delta})
        entry = dict(name=name, max_abs_error=maximum, allowance=limit,
                     reference_max=float(np.max(np.abs(native))), passed=maximum <= limit,
                     rms_error=float(np.sqrt(np.mean(np.abs(delta)**2))),
                     raw_mean_difference_real=np.mean(delta, axis=0).real.tolist(),
                     raw_mean_difference_imag=np.mean(delta, axis=0).imag.tolist(),
                     no_fitted_offset=True)
        (self.checks if gate else self.diagnostics).append(entry)
        return limit

    def mutation(self, name, field, reference, atol, rtol):
        delta = np.asarray(field)-np.asarray(reference)
        limit = float(atol+rtol*np.max(np.abs(reference)))
        maximum = float(np.max(np.abs(delta)))
        self.fields[name+"_mutant"] = field
        self.fields[name+"_difference"] = delta
        self.sensitivity.append(dict(name=name, max_abs_error=maximum, allowance=limit,
                                     exposed=maximum > limit,
                                     interpretation="Only separation above the preregistered allowance exposes this convention error."))


def literal_fourier_h0(report):
    """Independent finite trigonometric sums, not a Native action comparator."""
    shape = (5, 5, 5)
    mill = np.array([[0, 1, -1, 0], [0, 0, 1, -1], [0, 1, 0, 1]], dtype=int)
    reciprocal = np.diag([.7, .9, 1.1])
    fft = PhysicalRhoFFT.create(shape, mill, reciprocal, 17.)
    coordinates = np.indices(shape).reshape(3, -1, order="F") / np.array(shape)[:, None]
    waves = np.exp(2j*np.pi*(coordinates.T @ mill))
    x = np.array([[1+.3j, -.2j], [.4-.2j, .7], [-.3j, .1+.6j], [.8, -.4j]])
    report.compare("literal_fourier_inverse", fft.from_g(x), waves @ x, 2e-13, 2e-13)
    report.compare("literal_fourier_forward", fft.to_g(waves @ x), x, 2e-13, 2e-13)
    local = .3 + .2*np.cos(2*np.pi*coordinates[0]) + .17*np.sin(2*np.pi*coordinates[1])
    beta = np.array([[.1+.2j, .3], [-.2j, .4], [.7, -.1j], [-.2, .3+.1j]])
    d = np.array([[.8, -.17], [-.17, -.3]])
    kinetic = np.sum((reciprocal @ mill + np.array([.1, -.2, .05])[:, None])**2, axis=0)
    species = SpeciesProjectors("finite-algebra", (1,), beta, d)
    factors = PWFactors(fft, kinetic, (species,), np.zeros(3), 1, ({"literal": True},))
    background = SimpleNamespace(fft=fft, inputs=SimpleNamespace(sources=()), vrs=np.column_stack((local, local)))
    actual = MaterialHamiltonian(background, factors, 1, 1).apply(x)
    expected = kinetic[:, None]*x + waves.conj().T @ (local[:, None]*(waves @ x))/fft.nnr
    # Explicit projector-channel sums, independent of production reshape/einsum.
    for row in range(len(mill.T)):
        for col in range(x.shape[1]):
            expected[row, col] += sum(beta[row, a]*d[a, b]*sum(beta[j, b].conjugate()*x[j, col]
                                          for j in range(x.shape[0])) for a in range(2) for b in range(2))
    report.compare("literal_H0_action", actual, expected, 5e-13, 5e-13)
    report.mutation("literal_H0_wrong_FFT_normalization", kinetic[:, None]*x +
                    (expected-kinetic[:, None]*x)/fft.nnr, expected, 5e-13, 5e-13)
    # Literal Coulomb zero/near-zero rule, including a nonzero q below the
    # Native dimensionless threshold. No XC coefficients obscure the witness.
    density = np.ones((fft.nnr, 2), dtype=np.complex128)
    density[:, 1] *= .3j
    zeros = {name: np.zeros((fft.nnr, 2, 2)) for name in COEFFICIENT_NAMES}
    for label, qx in (("zero", 0.), ("below", .5e-4), ("above", 2e-4)):
        shifted = fft.at_q((qx, 0., 0.))
        actual = apply_spin_pbe_arrays(shifted, density, grho=np.zeros((3, fft.nnr, 2)),
                                       coefficients=zeros, tpiba=1.)
        expected = np.zeros_like(density)
        if qx*qx > 1e-8:
            expected[:] = 8*np.pi*(1+.3j)/(qx*qx)
        report.compare("literal_Hartree_threshold_"+label, actual, expected, 1e-9, 5e-13)


def native_checks(report, bundle, saved, bindings, backend):
    bg = bundle.background
    native_fft = NativeRhoFFT.from_bundle(bundle)
    if native_fft.physical != native_fft.storage:
        raise ValueError("original0018 oracle requires its actual unpadded serial storage")
    xml_binding = bindings[str(saved / "data-file-schema.xml")]
    ref = load_saved_reference(xml_binding)
    upfs = {s.name: bindings[str(saved / s.pseudo_file)] for s in ref.species}
    inputs = load_material_inputs(xml_binding, bindings[str(saved / "charge-density.dat")], upfs)
    rebuilt = reconstruct_reference(inputs, backend)
    report.sources.update({r["path"]: dict(r) for r in rebuilt.sources})
    report.compare("rho_total_magnetization", rebuilt.rho_total_magnetization,
                   bg.array("rho_total_magnetization"), 1e-10, 1e-10)
    report.compare("grho", rebuilt.grho, bg.array("grho"), 1e-10, 1e-10)
    report.compare("vrs", rebuilt.vrs, bg.array("vrs"), 1e-8, 1e-8)
    for name in COEFFICIENT_NAMES:
        report.compare(name, rebuilt.coefficients[name], bg.array(name), 1e-6, 1e-6)
    for name in ("rho_updown", "ionic_ry", "hartree_ry", "local_xc_ry", "gradient_flux_ry_bohr", "xc_ry"):
        report.fields["reconstructed_"+name] = getattr(rebuilt, name)
    report.diagnostics.append(dict(name="imaginary_roundoff_before_real_projection", maxima=dict(rebuilt.imaginary_maxima)))
    zero = np.flatnonzero(np.all(inputs.density.miller == 0, axis=0))
    if zero.size != 1:
        raise ValueError("unique saved density G0 required")
    ionic_g0 = sum(inputs.tables[name].local_g0 for name in ref.atom_names)
    report.mutation("wrong_zero_ionic_G0", rebuilt.vrs-ionic_g0, bg.array("vrs"), 1e-8, 1e-8)
    wrong_hartree_g0 = 8*np.pi*inputs.density.components[zero[0], 0].real
    report.mutation("wrong_extra_tpiba_gradient", rebuilt.grho*bg.metadata["tpiba"],
                    bg.array("grho"), 1e-10, 1e-10)
    report.mutation("wrong_finite_Hartree_G0_denominator", rebuilt.vrs+wrong_hartree_g0, bg.array("vrs"), 1e-8, 1e-8)
    report.mutation("wrong_Ry_to_Ha_vrs", rebuilt.vrs/2, bg.array("vrs"), 1e-8, 1e-8)
    report.mutation("omitted_negative_divergence", rebuilt.ionic_ry[:, None]+rebuilt.hartree_ry[:, None]+rebuilt.local_xc_ry,
                    bg.array("vrs"), 1e-8, 1e-8)
    report.mutation("wrong_divergence_sign", 2*(rebuilt.ionic_ry[:, None]+rebuilt.hartree_ry[:, None]+rebuilt.local_xc_ry)-rebuilt.vrs,
                    bg.array("vrs"), 1e-8, 1e-8)
    # Actual wrong background spin interpretation, not just relabeling a result.
    wrong_rho = np.column_stack((rebuilt.rho_total_magnetization.sum(axis=1),
                                rebuilt.rho_total_magnetization[:, 0]-rebuilt.rho_total_magnetization[:, 1]))
    wrong_gradient = rebuilt.fft.gradient(rebuilt.rho_total_magnetization).real
    wrong = backend.evaluate(wrong_rho, wrong_gradient)
    wrong_vrs = rebuilt.ionic_ry[:, None]+rebuilt.hartree_ry[:, None]+wrong.local_xc_ry-rebuilt.fft.divergence(wrong.gradient_flux_ry_bohr).real
    report.mutation("wrong_background_spin_basis", wrong_vrs, bg.array("vrs"), 1e-8, 1e-8)
    kernel_checks(report, bundle, rebuilt, "gamma")
    # Two explicit spin sectors, first saved k. Ingredient agreement is separate
    # from saved eigenpair residuals; no new convergence threshold is imposed.
    factor_cache = None
    for spin in (1, 2):
        wbinding = bindings[str(saved / f"wfc{'up' if spin == 1 else 'dw'}1.dat")]
        wfc = load_wavefunctions(wbinding)
        report.sources[wfc.source["path"]] = dict(wfc.source)
        if factor_cache is None or not factor_cache.matches(wfc):
            factor_cache = build_pw_factors(inputs, wfc)
        state = bundle.documents[f"state-{spin}-{(spin-1)*ref.nks+1}.json"]
        target = bg.array("mill")[:, state.array("igk_source").astype(int)-1]
        rows = {tuple(m): i for i, m in enumerate(wfc.miller.T)}
        order = np.array([rows[tuple(m)] for m in target.T])
        report.compare(f"spin{spin}_kinetic", factor_cache.kinetic_ry[order],
                       state.array("g2kin_target")[:wfc.igwx], 2e-12, 2e-12)
        offset = 0
        for species in factor_cache.species:
            count = species.beta.shape[1]
            report.compare(f"spin{spin}_{species.name}_packed_beta", species.beta[order],
                           state.array("vkb_source")[:wfc.igwx, offset:offset+count], 2e-10, 2e-10)
            for atom in species.atom_indices:
                nh = species.active_d.shape[0]
                report.compare(f"spin{spin}_atom{atom}_active_D", species.active_d,
                               state.array("deeq")[:nh, :nh, atom-1, spin-1], 0., 0.)
            offset += count
        diagnostic = MaterialHamiltonian(rebuilt, factor_cache, spin).saved_eigenpair_diagnostics(wfc)
        for name in ("action_ry", "residual_ry", "gram", "projected_h0_ry", "residual_norm_ry", "eigenvalues_ry"):
            report.fields[f"spin{spin}_saved_{name}"] = diagnostic[name]
        report.diagnostics.append(dict(name=f"spin{spin}_saved_eigenpair_residual",
                                       per_band_norm_ry=diagnostic["residual_norm_ry"].tolist(),
                                       qualification=diagnostic["qualification"], convergence_gate=False))
    return rebuilt


def kernel_checks(report, bundle, rebuilt, label):
    """New coefficients and preserved fixture caller on original complex probes."""
    bg = bundle.background
    native_fft = NativeRhoFFT.from_bundle(bundle)
    if native_fft.physical != native_fft.storage or native_fft.physical != rebuilt.fft.physical:
        raise ValueError("native K comparator must use the reconstructed physical grid")
    q = bg.metadata["tpiba"]*bg.array("q_cart_2pi_alat")
    report.diagnostics.append(dict(name=label+"_native_q", q_bohr_inv=q.tolist(),
                                   bundle=str(bundle.root), material_q_admission=False))
    probe_in, probe_out = bg.array("kernel_probe_in"), bg.array("kernel_probe_out")
    for probe in range(probe_in.shape[2]):
        n = probe_in[:, :, probe]
        actual = rebuilt.apply_kernel(n, q_bohr_inv=q)
        report.compare(f"{label}_Kprobe_{probe}", actual, probe_out[:, :, probe], 1e-7, 1e-6)
        report.compare(f"{label}_fixture_caller_Kprobe_{probe}", apply_spin_pbe_kernel(bundle, n),
                       probe_out[:, :, probe], 1e-7, 1e-6)
        report.fields[f"{label}_Kprobe_{probe}_input"] = n
        wrong_n = np.column_stack(((n[:, 0]+n[:, 1])/2, (n[:, 0]-n[:, 1])/2))
        report.mutation(f"{label}_Kprobe_{probe}_wrong_spin", rebuilt.apply_kernel(wrong_n, q_bohr_inv=q),
                        probe_out[:, :, probe], 1e-7, 1e-6)


def extra_kernel_bundle(report, binding_file, rebuilt, label):
    """Consume only an explicitly Root-bound original manifest; no producer."""
    data = binding_file.read_bytes()
    report.sources[str(binding_file.resolve())] = binding(binding_file, data)
    manifest_binding = json.loads(data)
    original_manifest = report.bound(manifest_binding)
    bundle = load_bundle(Path(manifest_binding["path"]).parent)
    if bundle.buffers["manifest.json"] != original_manifest:
        raise ValueError("extra original manifest changed between binding and decode")
    for name, buffer in bundle.buffers.items():
        receipt = binding(bundle.root/name, buffer)
        report.sources[receipt["path"]] = receipt
    # Reuse the single reconstruction only after checking that the independent
    # export really has this same physical background; never create a fake bundle.
    bg = bundle.background
    report.compare(label+"_same_density", rebuilt.rho_total_magnetization,
                   bg.array("rho_total_magnetization"), 1e-10, 1e-10)
    report.compare(label+"_same_grho", rebuilt.grho, bg.array("grho"), 1e-10, 1e-10)
    report.compare(label+"_same_vrs", rebuilt.vrs, bg.array("vrs"), 1e-8, 1e-8)
    for name in COEFFICIENT_NAMES:
        report.compare(label+"_same_"+name, rebuilt.coefficients[name], bg.array(name), 1e-6, 1e-6)
    kernel_checks(report, bundle, rebuilt, label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library-binding", type=Path, required=True,
                        help="JSON path/bytes/sha256 binding of Root-built shared backend")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--extra-manifest-binding", type=Path, action="append", default=[],
                        help="Root-admitted JSON binding of an original nonzero-q Native manifest; repeat for plus/minus")
    args = parser.parse_args()
    archive = args.output.with_suffix(".npz")
    if args.output.exists() or archive.exists():
        raise ValueError("choose fresh immutable JSON and NPZ output paths")
    report = Report()
    started = time.monotonic()
    fatal = None
    # Capture preloaded Python/binary source consumers before numerical work.
    execution_sources = {}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, "__file__", None)
        if filename and Path(filename).is_file():
            receipt = binding(filename)
            execution_sources[receipt["path"]] = receipt
    execution_sources[str(Path(__file__).resolve())] = binding(__file__)
    try:
        inventory = json.loads(report.bound(INVENTORY))
        bindings = {r["path"]: r for r in records(inventory) if Path(r["path"]).is_absolute()}
        manifest_binding = inventory["existing_comparator_manifest"]
        report.bound(manifest_binding)
        bundle = load_bundle(Path(manifest_binding["path"]).parent)
        # Compare actual decoded immutable buffers to original inventory records.
        for name, data in bundle.buffers.items():
            path = str(bundle.root / name)
            receipt = binding(path, data)
            original = bindings[path]
            if receipt["bytes"] != str(original["bytes"]) or receipt["sha256"] != original["sha256"]:
                raise ValueError(f"original0018 source identity mismatch: {path}")
            report.sources[path] = receipt
        library_data = args.library_binding.read_bytes()
        report.sources[str(args.library_binding.resolve())] = binding(args.library_binding, library_data)
        backend = NativePBE(json.loads(library_data))
        saved = bundle.root.parent.parent / "scratch/scf.save"
        literal_fourier_h0(report)
        rebuilt = native_checks(report, bundle, saved, bindings, backend)
        for number, binding_file in enumerate(args.extra_manifest_binding, start=1):
            extra_kernel_bundle(report, binding_file, rebuilt, f"extra_q_{number}")
    except Exception as exc:
        fatal = dict(type=type(exc).__name__, message=str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(archive, **report.fields)
    passed = fatal is None and bool(report.checks) and all(c["passed"] for c in report.checks)
    result = dict(kind="P254-0025-independent-saved-reference-reconstruction", passed=passed,
                  fatal=fatal, checks=report.checks, sensitivity=report.sensitivity,
                  diagnostics=report.diagnostics, sources=list(report.sources.values()),
                  execution_sources=list(execution_sources.values()), full_fields=binding(archive),
                  elapsed_seconds=time.monotonic()-started,
                  allowances=dict(density_gradient="1e-10*(1+reference_max)",
                                  vrs="1e-8 Ry+1e-8*reference_max",
                                  coefficients="1e-6+1e-6*reference_max",
                                  Kprobe="1e-7 Ry+1e-6*reference_max"),
                  scope="Original0018 Gamma fixture reconstruction plus explicitly bound optional native-q kernel comparisons, and preserved fixture caller; literal finite Fourier/H0 algebra. No directly exported Native unretained H action, material-q admission, convergence gate, physical preparation or full BvK claim.")
    text = json.dumps(result, indent=2, allow_nan=False)+"\n"
    args.output.write_text(text)
    print(text, end="", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
