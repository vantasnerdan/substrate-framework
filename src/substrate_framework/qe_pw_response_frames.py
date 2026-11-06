"""Build complete-PW, common-mu FD response frames from verified Native bundles.

The reference is the spectrum of reconstructed H0, not the saved SCF density.
Every representably nonzero FD occupation is retained; target spectra are never
truncated. No response algorithm or retained-band completion lives here.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from .pw_response_data import PWResponseChannel, PWResponseModel
from .qe_ionic_operator import BareIonicOperator
from .qe_pw_hamiltonian import assemble_target_hamiltonian
from .qe_response_bundle import QEResponseBundle, ResponseDocument, load_bundle
from .qe_response_fft import NativeRhoFFT


_COORD_ATOL = 2e-12
_N_ATOL = 1e-10


@dataclass(frozen=True)
class _Spectrum:
    energies: np.ndarray
    vectors: np.ndarray
    slots: np.ndarray
    hermiticity_error: float
    hermiticity_allowance: float
    normalization_error: float


def _bundle(value: QEResponseBundle | str | Path) -> QEResponseBundle:
    return value if isinstance(value, QEResponseBundle) else load_bundle(value)


def _same(a: np.ndarray, b: np.ndarray, label: str) -> None:
    if a.shape != b.shape:
        raise ValueError(f"actual source backgrounds/bases disagree in shape: {label}")
    if a.dtype.kind in "iu" or b.dtype.kind in "iu":
        equal = np.array_equal(a, b)
    else:
        allowance = 64 * np.finfo(float).eps * max(1., float(np.max(np.abs(a))))
        equal = np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and np.max(np.abs(a - b)) <= allowance
    if not equal:
        raise ValueError(f"actual source backgrounds/bases disagree beyond roundoff: {label}")


def _miller(bundle: QEResponseBundle, state: ResponseDocument, side: str) -> np.ndarray:
    name = f"igk_{side}"
    indices = state.array(name, (state.metadata[f"ngw_{side}"],), "int32")
    if state.descriptors[name].index_base != 1:
        raise ValueError(f"{state.path}:{name}: explicit Native one-based map required")
    return bundle.background.array("mill")[:, indices.astype(np.int64) - 1]


def _permutation(reference: np.ndarray, wanted: np.ndarray) -> np.ndarray:
    """Rows of reference coefficients in wanted Miller order (no folding)."""
    index = {tuple(column): row for row, column in enumerate(reference.T)}
    if len(index) != reference.shape[1] or wanted.shape != reference.shape:
        raise ValueError("source PW bases have different dimensions or duplicate Miller indices")
    try:
        order = np.array([index[tuple(column)] for column in wanted.T], dtype=np.int64)
    except KeyError as error:
        raise ValueError("true missing source reciprocal-basis mapping") from error
    if np.unique(order).size != order.size:
        raise ValueError("source Miller map is not bijective")
    return order


def _coordinate_equal(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.allclose(a, b, rtol=0, atol=_COORD_ATOL))


def _source_shift(bundle: QEResponseBundle, source: ResponseDocument,
                  state: ResponseDocument) -> np.ndarray | None:
    shift = np.linalg.solve(bundle.background.array("bg"),
                            state.array("k_source") - source.array("k_source"))
    integer = np.rint(shift)
    return integer.astype(np.int64) if _coordinate_equal(shift, integer) else None


def _match(bundle: QEResponseBundle, source: ResponseDocument) -> ResponseDocument:
    matches = [state for state in bundle.states
               if state.metadata["spin"] == source.metadata["spin"]
               and _source_shift(bundle, source, state) is not None]
    if len(matches) != 1:
        raise ValueError(f"{bundle.root}: need one actual reciprocal-equivalent source k/spin match, found {len(matches)}")
    return matches[0]


def _common_slots(bundle: QEResponseBundle, state: ResponseDocument,
                  shift: np.ndarray, fft: NativeRhoFFT) -> np.ndarray:
    miller = _miller(bundle, state, "target") + shift[:, None]
    coordinates = tuple(miller[i] % fft.physical[i] for i in range(3))
    slots = np.ravel_multi_index(coordinates, fft.storage, order="F")
    if np.unique(slots).size != slots.size:
        raise ValueError("reciprocal-envelope transform aliases distinct target plane waves")
    return slots


def _check_envelope(bundle: QEResponseBundle, state: ResponseDocument) -> None:
    q = bundle.background.array("q_cart_2pi_alat")
    shift = state.array("k_target") - state.array("k_source") - q
    if not _coordinate_equal(shift, np.zeros(3)):
        reciprocal = np.linalg.solve(bundle.background.array("bg"), shift)
        raise ValueError(f"{bundle.root}/{state.path}: actual ktarget != ksource+q; "
                         f"true missing reciprocal-envelope mapping, shift={reciprocal.tolist()}")


def _check_inputs(bundles: dict[str, QEResponseBundle]) -> NativeRhoFFT:
    gamma = bundles["gamma"]
    fft = NativeRhoFFT.from_bundle(gamma)
    if len({bundle.root for bundle in bundles.values()}) != 3:
        raise ValueError("Gamma/+q/-q require three distinct actual input bundles")
    if len({hashlib.sha256(bundle.buffers["manifest.json"]).hexdigest()
            for bundle in bundles.values()}) != 3:
        raise ValueError("copied identical manifests do not supply explicit reverse data")
    q = {label: bundle.background.array("q_cart_2pi_alat") for label, bundle in bundles.items()}
    if (not _coordinate_equal(q["gamma"], np.zeros(3))
            or _coordinate_equal(q["plus"], np.zeros(3))
            or not _coordinate_equal(q["minus"], -q["plus"])):
        raise ValueError("require actual Gamma and distinct opposite nonzero q branches")
    half = np.linalg.solve(gamma.background.array("bg"), q["plus"])
    if not _coordinate_equal(half, np.array([0.5, 0., 0.])):
        raise ValueError("half model requires the actual +half-b1 branch")
    for label, bundle in bundles.items():
        other = NativeRhoFFT.from_bundle(bundle)
        if (other.physical, other.storage, other.omega) != (fft.physical, fft.storage, fft.omega):
            raise ValueError("actual source FFT geometry differs between branches")
        for key in ("nelec", "alat", "tpiba", "nspin", "nat", "functional"):
            if bundle.background.metadata[key] != gamma.background.metadata[key]:
                raise ValueError(f"actual source background scalar differs: {key}")
        for name in ("g", "mill", "nl", "at", "bg", "vrs", "rho_total_magnetization", "ityp"):
            _same(bundle.background.array(name), gamma.background.array(name), name)
        if len(bundle.states) != len(gamma.states):
            raise ValueError(f"{label}: unmatched actual source-state count")
        for state in bundle.states:
            _check_envelope(bundle, state)
    return fft


def _check_source(gamma: QEResponseBundle, source: ResponseDocument,
                  bundle: QEResponseBundle, state: ResponseDocument) -> np.ndarray:
    shift = _source_shift(bundle, source, state)
    if shift is None:
        raise ValueError("source reciprocal shift is not an integer lattice vector")
    order = _permutation(_miller(gamma, source, "source"),
                         _miller(bundle, state, "source") + shift[:, None])
    n = source.metadata["ngw_source"]
    _same(source.array("vkb_source")[:n][order], state.array("vkb_source")[:n], "source projectors")
    for name in ("deeq", "nh", "ityp", "source_kweights"):
        _same(source.array(name), state.array(name), name)
    for key in ("spin", "ngauss", "degauss", "nkb"):
        if source.metadata[key] != state.metadata[key]:
            raise ValueError(f"actual source metadata differs: {key}")
    return order


def _spectrum_key(bundle: QEResponseBundle, state: ResponseDocument) -> bytes:
    digest = hashlib.sha256()
    digest.update(str(state.metadata["spin"]).encode())
    digest.update(bundle.background.descriptors["vrs"].sha256.encode())
    for array in (_miller(bundle, state, "target"), state.array("k_target"),
                  state.array("g2kin_target"),
                  state.array("vkb_target")[:state.metadata["ngw_target"]],
                  state.array("deeq"), state.array("nh"), state.array("ityp")):
        digest.update(str(array.shape).encode())
        digest.update(array.tobytes())
    return digest.digest()


def _diagonalize(item: tuple[QEResponseBundle, ResponseDocument]) -> _Spectrum:
    bundle, state = item
    h = assemble_target_hamiltonian(bundle, state)
    scale = max(1., float(np.max(np.abs(h))))
    error = float(np.max(np.abs(h - h.conj().T)))
    allowance = 64 * np.finfo(float).eps * h.shape[0] * scale
    if error > allowance:
        raise ValueError(f"{bundle.root}/{state.path}: H0 is not Hermitian ({error} > {allowance}); no repair")
    energies, vectors = eigh(h, overwrite_a=True, check_finite=False, driver="evd")
    norm_error = float(np.max(np.abs(np.sum(np.abs(vectors) ** 2, axis=0) - 1)))
    if norm_error > 64 * np.finfo(float).eps * vectors.shape[0]:
        raise ValueError("complete H0 eigenvectors failed column normalization")
    fft = NativeRhoFFT.from_bundle(bundle)
    indices = state.array("igk_target").astype(np.int64) - 1
    return _Spectrum(energies, vectors, fft.nl[indices], error, allowance, norm_error)


def _fermi(x: np.ndarray) -> np.ndarray:
    """Stable logistic retaining representable subnormal FD tails."""
    tail = np.exp(-np.abs(x))
    denominator = 1 + tail
    return np.where(x >= 0, 1 / denominator, tail / denominator)


def _chemical_potential(sources: list[tuple[ResponseDocument, _Spectrum]],
                        degauss: float, nelec: float) -> tuple[float, list[np.ndarray]]:
    weights = [float(state.array("source_kweights")[0]) for state, _ in sources]
    capacity = sum(weight * spectrum.energies.size for weight, (_, spectrum) in zip(weights, sources))
    if not 0 < nelec < capacity:
        raise ValueError("Native total electron count must lie inside full-PW FD capacity")

    def count(mu: float) -> float:
        return sum(weight * float(np.sum(_fermi((mu - spectrum.energies) / degauss)))
                   for weight, (_, spectrum) in zip(weights, sources))

    lo = min(float(spectrum.energies[0]) for _, spectrum in sources) - 64 * degauss
    hi = max(float(spectrum.energies[-1]) for _, spectrum in sources) + 64 * degauss
    mu = brentq(lambda value: count(value) - nelec, lo, hi,
                xtol=np.finfo(float).tiny, rtol=4 * np.finfo(float).eps)
    occupations = [_fermi((mu - spectrum.energies) / degauss) for _, spectrum in sources]
    if abs(count(mu) - nelec) > _N_ATOL:
        raise ValueError("common-mu full-FD total N failed absolute 1e-10 electron root accuracy")
    return float(mu), occupations


def _periodic_grid(fft: NativeRhoFFT, slots: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    padded = np.zeros((fft.nnr, vectors.shape[1]), dtype=np.complex128)
    padded[slots] = vectors
    grid = fft._embed(np.fft.ifftn(fft.physical_view(padded), axes=(0, 1, 2), norm="forward"))
    norms = np.mean(np.abs(fft.physical_view(grid)) ** 2, axis=(0, 1, 2))
    if np.max(np.abs(norms - 1), initial=0) > 1e-10:
        raise ValueError("source PW/real-grid Parseval normalization failed")
    return grid


def _ionic_actions(bundle: QEResponseBundle, state: ResponseDocument,
                   vectors: np.ndarray) -> np.ndarray:
    operator = BareIonicOperator.from_bundle(bundle, state)
    actions = np.empty((6, operator.target_slots.size, vectors.shape[1]), dtype=np.complex128)
    for axis in range(6):
        displacement = np.zeros(6)
        displacement[axis] = 1
        actions[axis] = operator.apply(displacement, vectors)
    return actions


def _density_diagnostics(fft: NativeRhoFFT, density: np.ndarray,
                         saved: np.ndarray) -> dict[str, Any]:
    reference = fft.physical_view(density)
    actual = fft.physical_view(saved)
    difference = reference - actual
    integral = fft.omega * np.mean(reference, axis=(0, 1, 2))
    return {"integral_up_down_electrons": integral.tolist(), "integral_total_electrons": float(integral.sum()),
            "saved_integral_up_down_electrons": (fft.omega * np.mean(actual, axis=(0, 1, 2))).tolist(),
            "difference_maxabs_up_down_bohr_minus3": np.max(np.abs(difference), axis=(0, 1, 2)).tolist(),
            "difference_rms_up_down_bohr_minus3": np.sqrt(np.mean(difference ** 2, axis=(0, 1, 2))).tolist(),
            "difference_L1_up_down_electrons": (fft.omega * np.mean(np.abs(difference), axis=(0, 1, 2))).tolist(),
            "difference_relative_L2": float(np.linalg.norm(difference.ravel()) / np.linalg.norm(actual.ravel())),
            "interpretation": "Measured full-FD H0-reference versus saved SCF-background applicability discrepancy; no convergence or material-accuracy claim and no discrepancy rejection threshold."}


def build_pw_response_models(
    gamma_bundle: QEResponseBundle | str | Path,
    plus_bundle: QEResponseBundle | str | Path,
    minus_bundle: QEResponseBundle | str | Path,
) -> dict[str, PWResponseModel]:
    """Return Gamma and explicit +/-half-b1 frames, sharing all source arrays.

    Uses at most six CPU workers; the caller must set BLAS/OpenMP threads to one.
    Inputs may be byte-verified typed bundles or paths loaded by the typed loader.
    No dense ionic operators, H0 copies, or on-disk frame caches are retained.
    """
    bundles = dict(zip(("gamma", "plus", "minus"), map(_bundle, (gamma_bundle, plus_bundle, minus_bundle))))
    fft = _check_inputs(bundles)
    gamma = bundles["gamma"]
    matches = [{label: _match(bundle, source) for label, bundle in bundles.items()}
               for source in gamma.states]
    permutations = [{label: _check_source(gamma, source, bundles[label], state)
                     for label, state in states.items()}
                    for source, states in zip(gamma.states, matches)]
    widths = {state.metadata["degauss"] for bundle in bundles.values() for state in bundle.states}
    if len(widths) != 1:
        raise ValueError("one actual common FD degauss is required across both spin channels")
    degauss = float(next(iter(widths)))
    unique: dict[bytes, tuple[QEResponseBundle, ResponseDocument]] = {}
    keys = []
    for states in matches:
        channel_keys = {}
        for label, state in states.items():
            key = _spectrum_key(bundles[label], state)
            unique.setdefault(key, (bundles[label], state))
            channel_keys[label] = key
        keys.append(channel_keys)
    with ThreadPoolExecutor(max_workers=6) as pool:
        spectra = dict(zip(unique, pool.map(_diagonalize, unique.values())))
    sources = [(source, spectra[key["gamma"]]) for source, key in zip(gamma.states, keys)]
    weight_sums = [sum(float(state.array("source_kweights")[0]) for state in gamma.states
                       if state.metadata["spin"] == spin) for spin in (1, 2)]
    if not np.allclose(weight_sums, [1., 1.], rtol=0, atol=1e-12):
        raise ValueError(f"actual explicit-spin source quadrature must sum to one per spin, got {weight_sums}")
    nelec = float(gamma.background.metadata["nelec"])
    mu, occupations = _chemical_potential(sources, degauss, nelec)
    channels: dict[str, list[PWResponseChannel]] = {"gamma": [], "half": []}
    density = np.zeros((fft.nnr, 2))
    source_records = []

    def channel_frames(index: int) -> tuple[PWResponseChannel, PWResponseChannel, dict[str, Any]]:
        source, spectrum = sources[index]
        states, orders, key = matches[index], permutations[index], keys[index]
        shifts = {label: _source_shift(bundles[label], source, state)
                  for label, state in states.items()}
        # Gamma target/source bases are verified, not assumed to have the same ordering.
        to_source = _permutation(_miller(gamma, source, "target"), _miller(gamma, source, "source"))
        _same(source.array("vkb_target")[:spectrum.energies.size][to_source],
              source.array("vkb_source")[:source.metadata["ngw_source"]],
              "Gamma target/source projector equivalence in mapped PW order")
        active = occupations[index] != 0
        vectors = spectrum.vectors[:, active]
        if not np.array_equal(to_source, np.arange(to_source.size)):
            vectors = vectors[to_source]
        energies = spectrum.energies[active]
        weighted = float(source.array("source_kweights")[0]) * occupations[index][active]
        slots = fft.nl[source.array("igk_source").astype(np.int64) - 1]
        grid = _periodic_grid(fft, slots, vectors)
        actions = {label: _ionic_actions(bundles[label], state, vectors if np.array_equal(
                    orders[label], np.arange(orders[label].size)) else vectors[orders[label]])
                   for label, state in states.items()}
        identity = f"spin{source.metadata['spin']}:ik{source.metadata['ik']}"
        common = dict(identity=identity, spin=source.metadata["spin"] - 1,
                      source_energies=energies, weighted_occupations=weighted, source_grid=grid)

        def frame(plus: str, minus: str) -> PWResponseChannel:
            p, m = spectra[key[plus]], spectra[key[minus]]
            return PWResponseChannel(**common, plus_energies=p.energies, plus_vectors=p.vectors,
                                     minus_energies=m.energies, minus_vectors=m.vectors,
                                     plus_slots=_common_slots(bundles[plus], states[plus], shifts[plus], fft),
                                     minus_slots=_common_slots(bundles[minus], states[minus], shifts[minus], fft),
                                     g_plus=actions[plus], g_minus=actions[minus])

        record = {"identity": identity, "native_spin": str(source.metadata["spin"]),
                  "native_ik": str(source.metadata["ik"]), "k_source_2pi_alat": source.array("k_source").tolist(),
                  "quadrature_weight": float(source.array("source_kweights")[0]),
                  "complete_source_dimension": str(spectrum.energies.size), "active_source_dimension": str(energies.size),
                  "active_rule": "f_FD != 0 in float64; no occupancy threshold", "occupancy_electrons": float(weighted.sum()),
                  "smallest_retained_FD": float(np.min(occupations[index][active], initial=1)),
                  "branches": {label: {"state_document": state.path, "native_ik": str(state.metadata["ik"]),
                                          "k_target_2pi_alat": state.array("k_target").tolist(),
                                          "source_reciprocal_shift": shifts[label].tolist(),
                                          "complete_target_dimension": str(spectra[key[label]].energies.size),
                                          "H0_hermiticity_maxabs_Ry": spectra[key[label]].hermiticity_error,
                                          "H0_hermiticity_roundoff_allowance_Ry": spectra[key[label]].hermiticity_allowance,
                                          "eigenvector_column_norm_maxabs_error": spectra[key[label]].normalization_error}
                               for label, state in states.items()}}
        return frame("gamma", "gamma"), frame("plus", "minus"), record

    with ThreadPoolExecutor(max_workers=6) as pool:
        for gamma_channel, half_channel, record in pool.map(channel_frames, range(len(sources))):
            channels["gamma"].append(gamma_channel)
            channels["half"].append(half_channel)
            source_records.append(record)
            density[:, gamma_channel.spin] += (np.abs(gamma_channel.source_grid) ** 2 @ gamma_channel.weighted_occupations) / fft.omega
    saved_total_m = gamma.background.array("rho_total_magnetization", (fft.nnr, 2), "float64")
    saved = np.column_stack(((saved_total_m[:, 0] + saved_total_m[:, 1]) / 2,
                             (saved_total_m[:, 0] - saved_total_m[:, 1]) / 2))
    diagnostics = _density_diagnostics(fft, density, saved)
    spin_counts = [sum(float(channel.weighted_occupations.sum()) for channel in channels["gamma"]
                       if channel.spin == spin) for spin in (0, 1)]
    if (abs(sum(spin_counts) - nelec) > _N_ATOL
            or np.max(np.abs(np.array(diagnostics["integral_up_down_electrons"]) - spin_counts)) > _N_ATOL):
        raise ValueError("full-FD density integral/explicit-spin occupancies fail total-N normalization")
    metadata = {"reference": "full-FD spectral reference of actual reconstructed H0; not exact saved SCF rho",
                "common_mu_Ry": mu, "actual_degauss_Ry": degauss, "native_nelec": nelec,
                "source_total_occupancy_electrons": sum(spin_counts), "source_spin_occupancies_electrons": spin_counts,
                "source_quadrature_sums_up_down": weight_sums, "density_diagnostics": diagnostics,
                "source_channels": source_records, "unique_complete_H0_spectra": str(len(spectra)),
                "assertions": {"verified_byte_closed_inputs": True, "explicit_distinct_actual_plus_minus": True,
                               "source_k_spin_and_Miller_maps_verified": True,
                               "source_reciprocal_gauge_and_target_envelope_transformed": True,
                               "common_mu_fixed_total_N": True, "all_nonzero_FD_sources": True,
                               "complete_target_spectra": True, "FFT_Parseval_and_density_integral": True},
                "input_manifests": {label: {"path": str(bundle.root / "manifest.json"),
                                            "sha256": hashlib.sha256(bundle.buffers["manifest.json"]).hexdigest()}
                                    for label, bundle in bundles.items()}}
    return {name: PWResponseModel(fft if name == "gamma" else NativeRhoFFT.from_bundle(bundles["plus"]),
                                 bundles["gamma" if name == "gamma" else "plus"].background.array("q_cart_2pi_alat"),
                                 tuple(values), {**metadata, "branch": name,
                                                 "reverse_sector": "same actual q=0 operator" if name == "gamma"
                                                 else "distinct actual minus-half-b1 bundle"})
            for name, values in channels.items()}
