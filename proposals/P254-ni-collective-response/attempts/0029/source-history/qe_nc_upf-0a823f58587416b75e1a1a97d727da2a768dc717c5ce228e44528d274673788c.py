"""Bound ordinary-NC UPF data and QE radial-to-angular channel indexing.

Source: QE 7.5 read_upf_new, read_pseudo, init_us_1 and newd_nous.
PP_LOCAL stays a radial Ry potential; PP_BETA stores r*beta. PP_DIJ is
radial-projector data, not an already angular-expanded solid-state matrix.
Only scalar Ni/H NC, no NLCC/SOC/PAW/USPP, is supported here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class NCUPF:
    path: str
    original_bytes: bytes
    sha256: str
    element: str
    z_valence: float
    r: NDArray[np.float64]
    rab: NDArray[np.float64]
    local_ry: NDArray[np.float64]
    beta_r: NDArray[np.float64]  # (mesh, nbeta), UPF r*beta convention
    angular_l: NDArray[np.int64]
    cutoff_indices: NDArray[np.int64]  # original one-based kbeta
    raw_dij: NDArray[np.float64]
    local_mesh: int
    beta_mesh: int
    channel_radial: NDArray[np.int64]  # indv - 1
    channel_l: NDArray[np.int64]
    channel_lm: NDArray[np.int64]  # nhtolm - 1
    expanded_d: NDArray[np.float64]

    @property
    def nbeta(self) -> int:
        return self.angular_l.size

    @property
    def nh(self) -> int:
        return self.channel_radial.size


def _array(node: ET.Element | None, count: int, label: str) -> NDArray[np.float64]:
    if node is None or node.text is None:
        raise ValueError(f"missing UPF {label}")
    values = np.array([float(x.replace("D", "E").replace("d", "e"))
                       for x in node.text.split()], dtype=np.float64)
    if values.size != count or not np.isfinite(values).all():
        raise ValueError(f"UPF {label}: wrong count or nonfinite values")
    if "size" in node.attrib and int(node.attrib["size"]) != count:
        raise ValueError(f"UPF {label}: inconsistent declared size")
    return values


def expand_channels(angular_l: NDArray[np.int64], raw_dij: NDArray[np.float64]
                    ) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """QE init_us_1 order: radial index, then m=1..2*l+1.

    Returns zero-based (indv, nhtol, nhtolm, active_D). For ordinary NC,
    newd_nous copies this D unchanged into every atom and explicit spin.
    """
    radial = np.repeat(np.arange(len(angular_l)), 2 * angular_l + 1)
    ell = angular_l[radial]
    lm = np.concatenate([np.arange(l*l, (l+1)*(l+1)) for l in angular_l])
    d = raw_dij[np.ix_(radial, radial)] * ((ell[:, None] == ell[None, :]) &
                                          (lm[:, None] == lm[None, :]))
    return radial, ell, lm, d


def read_nc_upf(path: str | Path, *, expected_sha256: str,
                expected_bytes: int | str) -> NCUPF:
    """Verify the very bytes decoded, then retain their original binding.

    r/positions are bohr, local and D follow the original UPF Ry convention.
    Radial integration limits are *different*: local msh follows the first
    r>10bohr then odd rule; beta kkbeta is max(kbeta), even values retained.
    """
    path = Path(path).resolve()
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if len(data) != int(expected_bytes) or digest != expected_sha256:
        raise ValueError(f"UPF original-byte binding mismatch: {path}")
    root = ET.fromstring(data)
    header = root.find("PP_HEADER")
    if header is None:
        raise ValueError("UPF v2 PP_HEADER required")
    h = header.attrib
    if h.get("pseudo_type") != "NC" or h.get("element", "").strip() not in ("Ni", "H"):
        raise ValueError("selected Ni/H ordinary NC UPF required")
    for key in ("is_ultrasoft", "is_paw", "is_coulomb", "has_so", "core_correction"):
        if h.get(key, "").strip().lower() not in ("f", "false", ".false."):
            raise ValueError(f"unsupported or unspecified UPF flag {key}")
    if h.get("relativistic") != "scalar" or h.get("functional", "").strip() != "PBE":
        raise ValueError("selected scalar-relativistic PBE UPF required")
    mesh, nbeta = int(h["mesh_size"]), int(h["number_of_proj"])
    if mesh < 3 or nbeta < 1:
        raise ValueError("radial mesh and nonlocal projectors required")
    r = _array(root.find("PP_MESH/PP_R"), mesh, "PP_R")
    rab = _array(root.find("PP_MESH/PP_RAB"), mesh, "PP_RAB")
    local = _array(root.find("PP_LOCAL"), mesh, "PP_LOCAL")
    if r[0] < 0 or np.any(np.diff(r) <= 0) or np.any(rab <= 0):
        raise ValueError("invalid UPF radial grid")
    nodes = [root.find(f"PP_NONLOCAL/PP_BETA.{i+1}") for i in range(nbeta)]
    if any(node is None for node in nodes):
        raise ValueError("missing indexed radial beta")
    ell = np.array([int(node.attrib["angular_momentum"]) for node in nodes], dtype=np.int64)
    cutoff = np.array([int(node.attrib["cutoff_radius_index"]) for node in nodes], dtype=np.int64)
    if np.any((ell < 0) | (ell > 2)) or np.any((cutoff < 3) | (cutoff > mesh)):
        raise ValueError("selected s/p/d channels and valid kbeta required")
    beta = np.column_stack([_array(node, mesh, "PP_BETA") for node in nodes])
    dij = _array(root.find("PP_NONLOCAL/PP_DIJ"), nbeta*nbeta, "PP_DIJ").reshape((nbeta, nbeta), order="F")
    if not np.array_equal(dij, dij.T):
        raise ValueError("real symmetric NC PP_DIJ required")
    past = np.flatnonzero(r > 10.0)
    msh = int(past[0] + 1) if past.size else mesh
    msh = 2 * ((msh + 1) // 2) - 1
    radial, channel_l, lm, expanded = expand_channels(ell, dij)
    for array in (r, rab, local, beta, ell, cutoff, dij, radial, channel_l, lm, expanded):
        array.setflags(write=False)
    return NCUPF(str(path), data, digest, h["element"].strip(), float(h["z_valence"]),
                 r, rab, local, beta, ell, cutoff, dij, msh, int(cutoff.max()),
                 radial, channel_l, lm, expanded)
