"""Full target-PW Hamiltonian of the registered 0016 generation-003 NC fixture.

In the *native target basis*, in Ry,

    H_ij = g2kin_i delta_ij + FFT_forward(vrs_spin)[G_i-G_j]
           + (B D B^dagger)_ij.

This is QE 7.5 PW/src/h_psi.f90 (kinetic initialization),
PW/src/vloc_psi_acc.f90:vloc_psi_k_acc (local multiplication), and
PW/src/add_vuspsi.f90:add_vuspsi_k (NC projector action).  D contains only
active deeq blocks, ordered species then atom as upflib/init_us_1.f90 sets
ofsbeta.  The forward FFT is normalized by the number of physical grid points;
there is no cell-volume or spin-multiplicity factor in this PW matrix.

No diagonalization, retained-band completion, eigenvalue fit, hermitization or
nonfinite-data repair is performed.  This is not a general electronic solver.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .qe_response_bundle import QEResponseBundle, ResponseDocument
from .qe_response_fft import NativeRhoFFT


def _ingredients(
    bundle: QEResponseBundle, state: ResponseDocument,
) -> tuple[NativeRhoFFT, NDArray[Any], NDArray[Any], NDArray[Any], NDArray[Any]]:
    """Return geometry, target Miller indices, kinetic diagonal, B and active D."""
    if not any(state is document for document in bundle.states):
        raise ValueError("state must be a target-state document of this verified bundle")
    bg = bundle.background
    if (bg.metadata["nat"] != 2 or bg.metadata["okvan"] is not False
            or bg.metadata["doublegrid"] is not False):
        raise ValueError("Hamiltonian scope is the two-atom ordinary-NC single-grid fixture")
    fft = NativeRhoFFT.from_bundle(bundle)
    n = state.metadata["ngw_target"]
    npwx, nkb, spin = (state.metadata[key] for key in ("npwx", "nkb", "spin"))
    indices = state.array("igk_target", (n,), "int32").astype(np.int64) - 1
    if state.descriptors["igk_target"].index_base != 1:
        raise ValueError("target PW map must retain Native one-based indices")
    mill = bg.array("mill", (3, bg.metadata["ngm"]), "int32").astype(np.int64)
    # Check the actual Native map, rather than guessing an FFT axis/order from g.
    slots = np.ravel_multi_index(tuple(mill[axis] % fft.physical[axis]
                                      for axis in range(3)), fft.storage, order="F")
    if not np.array_equal(slots, fft.nl):
        raise ValueError("Native mill and nl disagree on serial FFT placement")
    target_mill = mill[:, indices]
    if np.unique(slots[indices]).size != n:
        raise ValueError("target PW basis aliases on the actual FFT grid")
    kinetic = state.array("g2kin_target", (n,), "float64")
    beta = state.array("vkb_target", (npwx, nkb), "complex128")[:n]
    nh = state.array("nh", dtype="int32")
    ityp = state.array("ityp", (bg.metadata["nat"],), "int32")
    deeq = state.array("deeq", dtype="float64")
    d = np.zeros((nkb, nkb), dtype=np.float64)
    offset = 0
    for species, count_native in enumerate(nh, start=1):
        count = int(count_native)
        if count == 0:
            continue
        for atom, atom_species in enumerate(ityp):
            if atom_species != species:
                continue
            end = offset + count
            if end > nkb:
                raise ValueError("active species/atom projector blocks exceed Native nkb")
            d[offset:end, offset:end] = deeq[:count, :count, atom, spin - 1]
            offset = end
    if offset != nkb:
        raise ValueError("active species/atom projector counts do not exhaust Native nkb")
    return fft, target_mill, kinetic, beta, d


def _local_matrix(
    bundle: QEResponseBundle, state: ResponseDocument,
    fft: NativeRhoFFT, mill: NDArray[Any],
) -> NDArray[np.complex128]:
    potential = bundle.background.array("vrs", (fft.nnr, 2), "float64")
    spectrum = np.fft.fftn(fft.physical_view(potential[:, state.metadata["spin"] - 1]),
                           norm="forward")
    # A single small FFT followed by G_i-G_j lookup; never ngw full-grid FFTs.
    differences = tuple((mill[axis, :, None] - mill[axis, None, :]) % fft.physical[axis]
                        for axis in range(3))
    return spectrum[differences]


def assemble_target_hamiltonian(
    bundle: QEResponseBundle, state: ResponseDocument,
) -> NDArray[np.complex128]:
    """Return complex (ngw_target,ngw_target) H in exact Native target PW order.

    The complement of the exported eigenvectors remains present.  Only padded
    npwx rows and inactive deeq storage are excluded; no physical PW is dropped.
    Both input objects come from :func:`qe_response_bundle.load_bundle`.
    """
    fft, mill, kinetic, beta, d = _ingredients(bundle, state)
    h = _local_matrix(bundle, state, fft, mill)
    h[np.diag_indices(kinetic.size)] += kinetic
    h += (beta @ d) @ beta.conj().T
    if not np.isfinite(h).all():
        raise ValueError("assembled Hamiltonian is nonfinite; no repair is permitted")
    return h


def hamiltonian_sensitivity_actions(
    bundle: QEResponseBundle, state: ResponseDocument, vectors: NDArray[Any],
) -> dict[str, NDArray[np.complex128]]:
    """Apply three named scientific mutations without allocating mutated H copies.

    ``wrong_fft_normalization`` omits the 1/Ngrid forward normalization, not a
    fitted perturbation.  Exposure must be measured against Native eigenpairs;
    a mutation not exposed by those retained vectors supplies no sensitivity.
    """
    fft, mill, kinetic, beta, d = _ingredients(bundle, state)
    if vectors.ndim != 2 or vectors.shape[0] != kinetic.size:
        raise ValueError("sensitivity vectors must have ngw_target rows")
    local_action = _local_matrix(bundle, state, fft, mill) @ vectors
    nonlocal_action = beta @ (d @ (beta.conj().T @ vectors))
    kinetic_action = kinetic[:, None] * vectors
    return {
        "omitted_local": kinetic_action + nonlocal_action,
        "omitted_nonlocal": kinetic_action + local_action,
        "wrong_fft_normalization": (kinetic_action + np.prod(fft.physical) * local_action
                                    + nonlocal_action),
    }
