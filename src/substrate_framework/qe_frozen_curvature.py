"""Import QE irrep-zero frozen curvature from content-bound original XML bytes.

This is ionic plus fixed-wavefunction explicit curvature, not dyn00 alone and
not the screened dynamical matrix. The existing decimal XML is authoritative;
its import is not a bit-exact serialization of Native in-memory DP values.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from types import MappingProxyType
from typing import Any, Mapping
import xml.etree.ElementTree as ET

import numpy as np
from numpy.typing import NDArray

# Reuse the producer project's established complex-pair / Fortran-order reader.
# In particular, do not use its path-reading read_xml entry point here.
from catdft.qe_phonon_io import complex_record, complex_values, element, integer

SOURCE_AUTHORITY = {
    "total_frozen_hook": "PHonon/PH/dynmat0.f90:41-94",
    "restart_write": "PHonon/PH/ph_restart.f90:269-295",
    "pattern_columns": "PHonon/PH/ph_restart.f90:224-267",
    "cartesian_transform": "PHonon/PH/rotate_pattern_add.f90:50-85",
    "complex_serialization": "upflib/xmltools.f90:103-105,548-566",
    "mass_conversion_later": "PHonon/PH/dyndia.f90:61-75",
}


def verify_buffer(data: bytes, source: Mapping[str, Any]) -> dict[str, str]:
    """Bind the very buffer to be decoded, without filesystem access or rereads."""
    if not isinstance(data, bytes):
        raise TypeError("original immutable byte buffer required")
    receipt = {"path": str(source["path"]), "bytes": str(len(data)),
               "sha256": hashlib.sha256(data).hexdigest()}
    if str(source["bytes"]) != receipt["bytes"] or source["sha256"] != receipt["sha256"]:
        raise ValueError(f"source byte count/SHA256 mismatch: {receipt['path']}")
    return receipt


def _root(data: bytes) -> ET.Element:
    try:
        root = ET.fromstring(data)
    except ET.ParseError as error:
        raise ValueError("incomplete native XML byte buffer") from error
    if root.tag != "Root":
        raise ValueError("expected native XML Root")
    return root


def _boolean(parent: ET.Element, tag: str) -> bool:
    text = (element(parent, tag).text or "").strip().lower()
    if text in ("true", ".true.", "t"):
        return True
    if text in ("false", ".false.", "f"):
        return False
    raise ValueError(f"invalid native Boolean {tag}")


def _patterns(root: ET.Element, nat: int, iq: int) -> tuple[NDArray[np.complex128], tuple[tuple[int, int], ...]]:
    info = element(root, "IRREPS_INFO")
    if integer(info, "QPOINT_NUMBER") != iq:
        raise ValueError("pattern q-point number differs from source identity")
    nirr = integer(info, "NUMBER_IRR_REP")
    representations = [node for node in info if node.tag.startswith("REPRESENTION.")]
    if nirr < 1 or len(representations) != nirr:
        raise ValueError("incomplete native irreducible representation list")
    columns, ordering = [], []
    for irr in range(1, nirr + 1):
        representation = element(info, f"REPRESENTION.{irr}")
        count = integer(representation, "NUMBER_OF_PERTURBATIONS")
        perturbations = [node for node in representation if node.tag.startswith("PERTURBATION.")]
        if count < 1 or len(perturbations) != count:
            raise ValueError("incomplete native perturbation list")
        for ipert in range(1, count + 1):
            columns.append(complex_values(element(representation, f"PERTURBATION.{ipert}"),
                                          "DISPLACEMENT_PATTERN", (3 * nat,)))
            ordering.append((irr, ipert))
    if len(columns) != 3 * nat:
        raise ValueError("native patterns must span every atomic Cartesian displacement")
    u = np.column_stack(columns)
    if not np.allclose(u.conj().T @ u, np.eye(3 * nat), rtol=0, atol=1e-12):
        raise ValueError("XML patterns are not a complete unitary basis within decimal XML tolerance")
    return u, tuple(ordering)


def pattern_to_cartesian(curvature: NDArray[np.complex128],
                         patterns: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Apply the actual complex U C_pattern U†, without projection or repair."""
    if curvature.ndim != 2 or curvature.shape[0] != curvature.shape[1] or patterns.shape != curvature.shape:
        raise ValueError("curvature and full pattern basis must have equal square shapes")
    return patterns @ curvature @ patterns.conj().T


@dataclass(frozen=True)
class FrozenCurvature:
    """Original source buffers and unsymmetrized decimal-XML curvature in two bases."""

    pattern_matrix: NDArray[np.complex128]
    cartesian_matrix: NDArray[np.complex128]
    patterns: NDArray[np.complex128]
    pattern_order: tuple[tuple[int, int], ...]
    sources: tuple[Mapping[str, str], Mapping[str, str]]
    identity: Mapping[str, Any]
    buffers: tuple[bytes, bytes]
    minus_q_symmetry: bool
    q_group_rank: int

    def record(self) -> dict[str, Any]:
        """Return JSON-ready physical data, not frequencies or interpretations."""
        return {
            "source_identity": dict(self.identity), "sources": [dict(source) for source in self.sources],
            "done_irr": True, "irrep": 0, "shape": list(self.pattern_matrix.shape),
            "dtype": "complex128", "units": "Ry/bohr2", "mass_weighted": False,
            "induced_density_screening_included": False,
            "contribution": "irrep-zero: ionic plus fixed-wavefunction explicit terms",
            "cartesian_order": "atom-major: 3*(atom-1)+Cartesian_direction",
            "pattern_order": [{"mode": i, "irrep": irr, "perturbation": ipert}
                              for i, (irr, ipert) in enumerate(self.pattern_order, 1)],
            "serialization": {"complex": "consecutive real/imaginary pairs", "order": "F",
                              "matrix_writer_format": "(1p2es24.15)",
                              "patterns_writer": "native list-directed complex-vector real/imaginary pairs",
                              "precision_scope": "original decimal XML, not bit-exact in-memory Native DP"},
            "cartesian_transform": "U C_pattern Udagger",
            "pattern_matrix": complex_record(self.pattern_matrix),
            "cartesian_matrix": complex_record(self.cartesian_matrix),
            "patterns_U": complex_record(self.patterns),
            "native_pattern_flags": {"minus_q_sym": self.minus_q_symmetry,
                                     "qpoint_group_rank": self.q_group_rank},
            "parser_projection": "none; no symmetrization, mass division or acoustic-sum-rule repair",
            "source_authority": SOURCE_AUTHORITY,
        }


def parse_frozen_curvature(
    dynmat_bytes: bytes, patterns_bytes: bytes, *,
    dynmat_source: Mapping[str, Any], patterns_source: Mapping[str, Any],
    identity: Mapping[str, Any],
) -> FrozenCurvature:
    """Parse bound buffers; caller supplies independently bound q/cell/bundle identity.

    Native restart dynmat XML does not encode q or cell itself. ``identity`` must
    therefore name nat, iq and the corresponding bundle/cell/q provenance; those
    identities must not be inferred from numerical curvature entries.
    """
    receipts = (verify_buffer(dynmat_bytes, dynmat_source),
                verify_buffer(patterns_bytes, patterns_source))
    nat, iq = identity["nat"], identity["iq"]
    if type(nat) is not int or nat < 1 or type(iq) is not int or iq < 1:
        raise ValueError("positive integer nat and iq source identities required")
    for key in ("q_cart_2pi_alat", "cell", "bundle_manifest", "applicability_flags"):
        if key not in identity:
            raise ValueError(f"missing bound source identity: {key}")
    root = _root(dynmat_bytes)
    if not _boolean(element(root, "PM_HEADER"), "DONE_IRR"):
        raise ValueError("irrep-zero curvature is not DONE_IRR")
    cpattern = complex_values(element(root, "PARTIAL_MATRIX"), "PARTIAL_DYN", (3 * nat, 3 * nat))
    pattern_root = _root(patterns_bytes)
    u, ordering = _patterns(pattern_root, nat, iq)
    info = element(pattern_root, "IRREPS_INFO")
    ccart = pattern_to_cartesian(cpattern, u)
    if not np.isfinite(ccart).all():
        raise ValueError("nonfinite Cartesian curvature")
    for matrix in (cpattern, ccart, u):
        matrix.flags.writeable = False
    return FrozenCurvature(cpattern, ccart, u, ordering,
                           tuple(MappingProxyType(receipt) for receipt in receipts),
                           MappingProxyType(dict(identity)), (dynmat_bytes, patterns_bytes),
                           _boolean(info, "MINUS_Q_SYM"), integer(info, "QPOINT_GROUP_RANK"))
