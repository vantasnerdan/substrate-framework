"""Source-rank validation for the ordinary-NC response export, without a solver."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .qe_response_fft import NativeRhoFFT

if TYPE_CHECKING:
    from .qe_response_bundle import QEResponseBundle, ResponseDocument


def _positive(doc: ResponseDocument, *keys: str) -> tuple[int, ...]:
    values = tuple(doc.metadata[key] for key in keys)
    if any(type(value) is not int or value < 1 for value in values):
        raise ValueError(f"{doc.path}: positive native dimensions/identities required: {keys}")
    return values


def _indices(doc: ResponseDocument, name: str, size: int | None = None) -> np.ndarray:
    array = doc.array(name, None if size is None else (size,), "int32")
    if array.ndim != 1 or doc.descriptors[name].index_base != 1 or np.any(array < 1):
        raise ValueError(f"{doc.path}:{name}: positive one-based native index vector required")
    return array


def validate_bundle(bundle: QEResponseBundle) -> None:
    """Retain padded PW matrices and distinguish quadrature weights from f_FD.

    This validates ingredients, not completeness of the electronic complement,
    irrep completion, a Hamiltonian solve, or a finite-frequency response.
    """
    fft = NativeRhoFFT.from_bundle(bundle)
    bg = bundle.background
    ngm, nat, nksq, npwx = _positive(bg, "ngm", "nat", "nksq", "npwx")
    if nat != 2:
        raise ValueError("registered response background requires two atoms")
    if bg.metadata["okvan"] is not False or bg.metadata["doublegrid"] is not False:
        raise ValueError("registered response background requires ordinary NC and one FFT grid")
    if bg.metadata["density_basis"] != "up,down" or bg.metadata["energy_units"] != "Ry":
        raise ValueError("native response spin basis / energy units disagree")
    for name in ("dmuxc", "dvxc_rr", "dvxc_sr", "dvxc_ss", "dvxc_s"):
        bg.array(name, (fft.nnr, 2, 2), "float64")
    bg.array("grho", (3, fft.nnr, 2), "float64")
    bg.array("vrs", (fft.nnr, 2), "float64")
    bg.array("mill", (3, ngm), "int32")
    bg.array("at", (3, 3), "float64")
    bg.array("bg", (3, 3), "float64")
    bg.array("kernel_probe_in", (fft.nnr, 2, 3), "complex128")
    bg.array("kernel_probe_out", (fft.nnr, 2, 3), "complex128")
    bg.array("kernel_probe_modes", (3, 3), "int32")
    bg.array("kernel_probe_coefficients", (2, 3), "complex128")
    states = {}
    for state in bundle.states:
        spin, ik, source_n, target_n, state_npwx, nbnd, nkb = _positive(
            state, "spin", "ik", "ngw_source", "ngw_target", "npwx", "nbnd", "nkb")
        if spin > 2 or ik > nksq or state_npwx != npwx or max(source_n, target_n) > npwx:
            raise ValueError(f"{state.path}: invalid native spin/k/PW dimensions")
        identity = (spin, ik)
        if identity in states:
            raise ValueError("duplicate native state identity")
        states[identity] = state
        for name in ("evc", "evq"):
            state.array(name, (npwx, nbnd), "complex128")
        for name, size in (("igk_source", source_n), ("igk_target", target_n)):
            indices = _indices(state, name, size)
            if np.any(indices > ngm) or np.unique(indices).size != size:
                raise ValueError(f"{state.path}:{name}: invalid global G mapping")
        state.array("g2kin_target", (target_n,), "float64")
        state.array("vkb_target", (npwx, nkb), "complex128")
        nh = state.array("nh", dtype="int32")
        ityp = _indices(state, "ityp", nat)
        deeq = state.array("deeq", dtype="float64")
        if nh.ndim != 1 or np.any(nh < 0) or np.any(ityp > nh.size):
            raise ValueError("invalid NC species/projector counts")
        if (deeq.ndim != 4 or deeq.shape[0] != deeq.shape[1] or
                deeq.shape[2:] != (nat, 2) or deeq.shape[0] < np.max(nh)):
            raise ValueError("deeq must retain actual (nhm,nhm,nat,nspin) rank")
        for name in ("eigenvalues_source", "eigenvalues_target", "weighted_native_wg",
                     "weighted_native_wg_target"):
            state.array(name, (nbnd,), "float64")
        for name in ("source_kweights", "target_kweights"):
            weights = state.array(name, (1,), "float64")
            if np.any(weights < 0):
                raise ValueError("native quadrature weights must be nonnegative")
        for name in ("k_source", "k_target"):
            state.array(name, (3,), "float64")
        for name in ("ikks", "ikqs"):
            _indices(state, name, nksq)
        # isk is a native-capacity k-to-spin label array, not a vector of
        # k indices. Only the referenced physical k slots identify states.
        isk = state.array("isk", dtype="int32")
        if isk.ndim != 1 or state.descriptors["isk"].index_base != 1:
            raise ValueError("isk must retain the one-based native k-to-spin map")
        ikk, ikq = _positive(state, "ikk", "ikq")
        if max(ikk, ikq) > isk.size or isk[ikk - 1] != spin or isk[ikq - 1] != spin:
            raise ValueError("native k/partner/spin identity mismatch")
        if state.metadata["ngauss"] != -99 or state.metadata["degauss"] <= 0:
            raise ValueError("registered electronic reference requires Fermi-Dirac smearing")
        if not np.isfinite(state.metadata["ef"]) or not np.isfinite(state.metadata["degauss"]):
            raise ValueError("nonfinite reference Fermi energy / smearing width")
        if state.metadata["energy_units"] != "Ry" or not isinstance(state.metadata["smearing"], str):
            raise ValueError("electronic energy/smearing convention missing")
    bare_identities = set()
    for bare in bundle.bare:
        spin, ik, target_n, bare_npwx, nbnd, mode = _positive(
            bare, "spin", "ik", "ngw_target", "npwx", "nbnd", "mode")
        _positive(bare, "irrep", "perturbation_index")
        state = states.get((spin, ik))
        if state is None or bare_npwx != npwx or mode > 3 * nat:
            raise ValueError("bare RHS does not identify an exported native electronic state")
        if target_n != state.metadata["ngw_target"] or nbnd != state.metadata["nbnd"]:
            raise ValueError("bare RHS/state physical dimensions disagree")
        if (spin, ik, mode) in bare_identities:
            raise ValueError("duplicate bare perturbation/state identity")
        bare_identities.add((spin, ik, mode))
        bare.array("fullPWbareRHS_dvpsi_before_screening", (npwx, nbnd), "complex128")
        bare.array("Cartesian_complex_pattern", (3 * nat,), "complex128")
