"""Factorized bare Cartesian ionic action in the actual Native source/target PWs.

Scope: selected ordinary-NC/collinear-PBE/no-NLCC/no-Hubbard/two-atom,
serial physical 24^3 dense=smooth backgrounds. Source authority is QE 7.5
compute_dvloc, dvqpsi_us_only, init_us_2 and init_us_2_acc. Native beta already
contains exp[-i(k+G).tau]; displacement amplitudes are linear, not conjugated.
There is no retained-band closure, fitted vertex, extra atom-q phase,
hermitization, nonfinite-data repair, or independent Native-complement claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .qe_response_bundle import QEResponseBundle, ResponseDocument
from .qe_response_fft import NativeRhoFFT


@dataclass(frozen=True)
class _AtomProjectors:
    atom: int
    source: NDArray[Any]
    target: NDArray[Any]
    d: NDArray[Any]


def _pw_map(bundle: QEResponseBundle, state: ResponseDocument,
            fft: NativeRhoFFT, side: str) -> NDArray[np.int64]:
    n = state.metadata[f"ngw_{side}"]
    name = f"igk_{side}"
    indices = state.array(name, (n,), "int32").astype(np.int64) - 1
    if state.descriptors[name].index_base != 1:
        raise ValueError(f"{name} must retain explicit Native one-based indices")
    if (np.any(indices < 0) or np.any(indices >= fft.nl.size)
            or np.unique(fft.nl[indices]).size != n):
        raise ValueError(f"{side} physical PW map is invalid or aliases on the FFT grid")
    return indices


def _atom_blocks(bundle: QEResponseBundle, state: ResponseDocument,
                 ns: int, nt: int) -> tuple[_AtomProjectors, ...]:
    bg = bundle.background
    npwx, nkb, spin = (state.metadata[key] for key in ("npwx", "nkb", "spin"))
    source = state.array("vkb_source", (npwx, nkb), "complex128")[:ns]
    target = state.array("vkb_target", (npwx, nkb), "complex128")[:nt]
    nh = state.array("nh", dtype="int32")
    ityp = bg.array("ityp", (2,), "int32")
    if (bg.descriptors["ityp"].index_base != 1
            or not np.array_equal(ityp, state.array("ityp", (2,), "int32"))):
        raise ValueError("background/state Native species identities disagree")
    deeq = state.array("deeq", dtype="float64")
    blocks = []
    offset = 0
    # init_us_2_acc and dvqpsi_us_only both order species, then actual atom.
    for species, native_count in enumerate(nh, start=1):
        count = int(native_count)
        for atom, atom_species in enumerate(ityp):
            if atom_species != species or count == 0:
                continue
            end = offset + count
            if end > nkb:
                raise ValueError("active atom projector blocks exceed Native nkb")
            blocks.append(_AtomProjectors(atom, source[:, offset:end],
                                          target[:, offset:end],
                                          deeq[:count, :count, atom, spin - 1]))
            offset = end
    if offset != nkb:
        raise ValueError("active species/atom blocks do not exhaust Native nkb")
    return tuple(blocks)


@dataclass(frozen=True)
class BareIonicOperator:
    """Apply six atom-major Cartesian amplitudes to any number of source PWs.

    Construct with ``from_bundle`` after byte-closed loading. ``apply(u, psi)``
    accepts complex (6,) u in displacement amplitudes and complex
    (ngw_source,nvectors) psi, returning (ngw_target,nvectors) in Ry/bohr.
    Only physical PW rows are accepted or returned; npwx padding is not a basis.
    """

    fft: NativeRhoFFT
    source_slots: NDArray[Any]
    target_slots: NDArray[Any]
    local_cartesian: NDArray[Any]
    p_source: NDArray[Any]
    p_target: NDArray[Any]
    atoms: tuple[_AtomProjectors, ...]

    @classmethod
    def from_bundle(cls, bundle: QEResponseBundle,
                    state: ResponseDocument) -> BareIonicOperator:
        if not any(state is document for document in bundle.states):
            raise ValueError("state must belong to this verified bundle")
        bg = bundle.background
        if (bg.metadata["nat"] != 2 or bg.metadata["okvan"] is not False
                or bg.metadata["doublegrid"] is not False
                or bg.metadata["lda_plus_u"] is not False):
            raise ValueError("ionic scope requires two-atom ordinary NC/no Hubbard/single grid")
        fft = NativeRhoFFT.from_bundle(bundle)
        source = _pw_map(bundle, state, fft, "source")
        target = _pw_map(bundle, state, fft, "target")
        mill = bg.array("mill", (3, bg.metadata["ngm"]), "int32")
        slots = np.ravel_multi_index(tuple(mill[a] % fft.physical[a] for a in range(3)),
                                    fft.storage, order="F")
        if not np.array_equal(slots, fft.nl):
            raise ValueError("actual Native mill/nl FFT storage maps disagree")
        local = bg.array("vloc_bare_cartesian", (fft.nnr, 6), "complex128")
        # Native unit descriptors may append a semicolon-delimited basis note.
        if bg.descriptors["vloc_bare_cartesian"].units.partition(";")[0].strip() != "Ry/bohr":
            raise ValueError("bare Cartesian local fields must be in Ry/bohr")
        g = bg.array("g", (3, bg.metadata["ngm"]), "float64")
        tpiba = bg.metadata["tpiba"]
        ps = tpiba * (g[:, source].T + state.array("k_source", (3,), "float64"))
        pt = tpiba * (g[:, target].T + state.array("k_target", (3,), "float64"))
        return cls(fft, fft.nl[source], fft.nl[target], local, ps, pt,
                   _atom_blocks(bundle, state, source.size, target.size))

    def _inputs(self, displacement: NDArray[Any],
                vectors: NDArray[Any]) -> tuple[NDArray[Any], NDArray[Any]]:
        u = np.asarray(displacement, dtype=np.complex128)
        psi = np.asarray(vectors, dtype=np.complex128)
        if u.shape != (6,) or psi.ndim != 2 or psi.shape[0] != self.source_slots.size:
            raise ValueError("expected displacement (6,) and source vectors (ngw_source,nvectors)")
        if not np.isfinite(u).all() or not np.isfinite(psi).all():
            raise ValueError("nonfinite ionic input; no repair is permitted")
        return u, psi

    def _local(self, u: NDArray[Any], psi: NDArray[Any]) -> NDArray[np.complex128]:
        # Work on physical serial FFT cells while retaining actual storage maps.
        padded = np.zeros((self.fft.nnr, psi.shape[1]), dtype=np.complex128)
        padded[self.source_slots] = psi
        envelope = np.fft.ifftn(self.fft.physical_view(padded), axes=(0, 1, 2),
                               norm="forward")
        potential = self.fft.physical_view(self.local_cartesian @ u)
        envelope *= potential[..., None]
        spectrum = np.fft.fftn(envelope, axes=(0, 1, 2), norm="forward")
        coordinates = np.unravel_index(self.target_slots, self.fft.storage, order="F")
        return spectrum[coordinates]

    def _nonlocal(self, u: NDArray[Any], psi: NDArray[Any], *,
                  wrong_target_sign: bool = False,
                  conjugated_projectors: bool = False) -> NDArray[np.complex128]:
        out = np.zeros((self.target_slots.size, psi.shape[1]), dtype=np.complex128)
        for block in self.atoms:
            ua = u[3 * block.atom:3 * block.atom + 3]
            source, target = block.source, block.target
            if conjugated_projectors:
                source, target = source.conj(), target.conj()
            # Contract momenta with complex u, never with u.conj(). This is
            # i beta_t D beta_s^dagger (p_s.u psi)
            # - i (p_t.u) beta_t D beta_s^dagger psi.
            derivative_source = source.conj().T @ ((self.p_source @ ua)[:, None] * psi)
            overlap = source.conj().T @ psi
            out += 1j * (target @ (block.d @ derivative_source))
            sign = 1j if wrong_target_sign else -1j
            out += sign * (self.p_target @ ua)[:, None] * (target @ (block.d @ overlap))
        return out

    def apply(self, displacement: NDArray[Any],
              vectors: NDArray[Any]) -> NDArray[np.complex128]:
        """Return the full physical target-PW action, without a dense PW matrix."""
        u, psi = self._inputs(displacement, vectors)
        out = self._local(u, psi) + self._nonlocal(u, psi)
        if not np.isfinite(out).all():
            raise ValueError("nonfinite bare ionic action; no repair is permitted")
        return out

    def sensitivity_actions(self, displacement: NDArray[Any],
                            vectors: NDArray[Any]) -> dict[str, NDArray[Any]]:
        """Scientific mutations; exposure is decided against actual Native RHS.

        Conjugated projectors reverse their actual complex phases, including the
        k+G atomic structure factors and angular (-i)^l phases. It is not claimed
        to isolate only the atomic factor. Unexposed mutations earn no credit.
        """
        u, psi = self._inputs(displacement, vectors)
        local, nonlocal_action = self._local(u, psi), self._nonlocal(u, psi)
        return {
            "dropped_local": nonlocal_action,
            "dropped_nonlocal": local,
            "wrong_target_derivative_sign": local + self._nonlocal(u, psi, wrong_target_sign=True),
            "conjugated_projector_phases": local + self._nonlocal(u, psi, conjugated_projectors=True),
            "conjugated_displacement": self.apply(u.conj(), psi),
        }
