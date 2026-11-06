"""QE 7.5 saved XML geometry and reference metadata, not a quantum preparation.

XML energies are Hartree; converted Ry fields are explicitly separate. Raw
Marzari-Vanderbilt occupations are retained unchanged, never clipped or called
probabilities. XML k and reciprocal columns use 2pi/alat; binary saved records
already use bohr^-1. Geometry uses bohr, with cell vectors stored as columns.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Mapping
import xml.etree.ElementTree as ET

import numpy as np

from .qe_saved_native import read_bound_bytes
from .qe_saved_symmetry import SavedSymmetries, parse_saved_symmetries


def _node(parent, name):
    result = parent.find(name)
    if result is None:
        raise ValueError(f"missing QE XML element {name}")
    return result


def _text(parent, name):
    return (_node(parent, name).text or "").strip()


def _bool(parent, name):
    text = _text(parent, name)
    if text not in ("true", "false"):
        raise ValueError(f"invalid QE XML Boolean {name}")
    return text == "true"


def _values(node, count):
    values = np.array([float(x) for x in (node.text or "").split()], dtype=np.float64)
    if values.size != count or not np.isfinite(values).all():
        raise ValueError(f"invalid QE XML vector {node.tag}")
    if "size" in node.attrib and int(node.attrib["size"]) != count:
        raise ValueError(f"QE XML declared size mismatch: {node.tag}")
    values.setflags(write=False)
    return values


def _readonly(values):
    values.setflags(write=False)
    return values


def _grid(basis, name):
    node = _node(basis, name)
    result = tuple(int(node.attrib[f"nr{i}"]) for i in (1, 2, 3))
    if min(result) < 1:
        raise ValueError("invalid physical FFT dimensions")
    return result


@dataclass(frozen=True)
class SavedSpecies:
    name: str
    mass_amu: float
    pseudo_file: str


@dataclass(frozen=True)
class SavedKPoint:
    ik: int                             # spin-local index, original XML order
    npw: int
    weight: float                       # native quadrature, not mode probability
    k_cart_2pi_alat: np.ndarray          # (3,), original XML values
    k_cart_bohr_inv: np.ndarray          # (3,), explicitly converted
    eigenvalues_ha: tuple[np.ndarray, ...]  # spin-order arrays (nbnd_spin,)
    eigenvalues_ry: tuple[np.ndarray, ...]
    raw_occupations: tuple[np.ndarray, ...] # original values, no filtering


@dataclass(frozen=True)
class SavedReference:
    source: Mapping[str, str]
    creator_version: str
    source_energy_units: str
    nat: int
    species: tuple[SavedSpecies, ...]
    atom_names: tuple[str, ...]
    atom_indices: tuple[int, ...]         # original one-based indices
    positions_bohr: np.ndarray           # (3,nat), original atom order
    cell_bohr: np.ndarray                # (3,3), columns
    alat_bohr: float
    reciprocal_2pi_alat: np.ndarray      # (3,3), original XML columns
    reciprocal_bohr_inv: np.ndarray      # (3,3), converted columns
    fft_grid: tuple[int, int, int]
    fft_smooth: tuple[int, int, int]
    fft_box: tuple[int, int, int]
    ngm: int
    ngms: int
    npwx: int
    gamma_only: bool
    lsda: bool
    noncolin: bool
    spinorbit: bool
    nbnd_by_spin: tuple[int, ...]
    nks: int                             # XML spin-local k count
    kpoints: tuple[SavedKPoint, ...]
    nelec: float
    occupations_kind: str
    smearing: str | None
    degauss_ha: float | None
    degauss_ry: float | None
    fermi_energies_ha: np.ndarray
    fermi_energies_ry: np.ndarray
    total_energies_ha: Mapping[str, float]
    total_energies_ry: Mapping[str, float]
    ecutwfc_ha: float
    ecutrho_ha: float
    monkhorst_pack: Mapping[str, int]
    saved_symmetries: SavedSymmetries | None = None

    def select(self, ik: int, ispin: int) -> tuple[SavedKPoint, np.ndarray, np.ndarray]:
        """Return k record, Ry energies, raw occupations for exact local k/spin."""
        if not 1 <= ik <= self.nks or not 1 <= ispin <= len(self.nbnd_by_spin):
            raise ValueError("saved reference requires one-based spin-local ik and explicit ispin")
        point = self.kpoints[ik - 1]
        return point, point.eigenvalues_ry[ispin - 1], point.raw_occupations[ispin - 1]


def _kpoint(node, ik, band_counts, tpiba):
    k_node = _node(node, "k_point")
    k = _values(k_node, 3)
    energies = _values(_node(node, "eigenvalues"), sum(band_counts))
    occupations = _values(_node(node, "occupations"), sum(band_counts))
    offsets = np.cumsum((0,) + band_counts)
    slices = tuple(slice(int(a), int(b)) for a, b in zip(offsets[:-1], offsets[1:]))
    ha = tuple(energies[s] for s in slices)
    ry = tuple(_readonly(2.0 * e) for e in ha)
    raw = tuple(occupations[s] for s in slices)
    npw, weight = int(_text(node, "npw")), float(k_node.attrib["weight"])
    if npw < 1 or not np.isfinite(weight):
        raise ValueError("invalid saved k-point metadata")
    return SavedKPoint(ik, npw, weight, k, _readonly(tpiba * k), ha, ry, raw)


def load_saved_reference(binding: Mapping[str, Any]) -> SavedReference:
    data, source = read_bound_bytes(binding, max_bytes=16 << 20)
    root = ET.fromstring(data)
    if root.tag != "{http://www.quantum-espresso.org/ns/qes/qes-1.0}espresso" or root.attrib.get("Units") != "Hartree atomic units":
        raise ValueError("unsupported saved XML envelope/units")
    creator = _node(root, "general_info/creator")
    if creator.attrib.get("NAME") != "PWSCF" or creator.attrib.get("VERSION") != "7.5":
        raise ValueError("saved reference supports QE PWSCF 7.5 only")
    out = _node(root, "output")
    structure, basis, bands = (_node(out, n) for n in ("atomic_structure", "basis_set", "band_structure"))
    nat, alat = int(structure.attrib["nat"]), float(structure.attrib["alat"])
    if nat < 1 or not np.isfinite(alat) or alat <= 0:
        raise ValueError("invalid saved geometry nat/alat")
    atoms = _node(structure, "atomic_positions").findall("atom")
    indices = tuple(int(a.attrib["index"]) for a in atoms)
    if len(atoms) != nat or indices != tuple(range(1, nat + 1)):
        raise ValueError("incomplete or reordered original atom index sequence")
    names = tuple(a.attrib["name"] for a in atoms)
    positions = _readonly(np.column_stack([_values(a, 3) for a in atoms]))
    cell = _readonly(np.column_stack([_values(_node(structure, f"cell/a{i}"), 3) for i in (1, 2, 3)]))
    species_node = _node(out, "atomic_species")
    species = tuple(SavedSpecies(s.attrib["name"], float(_text(s, "mass")), _text(s, "pseudo_file")) for s in species_node.findall("species"))
    if len(species) != int(species_node.attrib["ntyp"]) or any(n not in {s.name for s in species} for n in names):
        raise ValueError("saved species/atom identity mismatch")
    lsda, noncolin, spinorbit = (_bool(bands, n) for n in ("lsda", "noncolin", "spinorbit"))
    if noncolin or spinorbit:
        raise ValueError("spinor saved reference unsupported")
    counts = tuple(int(_text(bands, n)) for n in (("nbnd_up", "nbnd_dw") if lsda else ("nbnd",)))
    if min(counts) < 1:
        raise ValueError("invalid saved band count")
    reciprocal = _readonly(np.column_stack([_values(_node(basis, f"reciprocal_lattice/b{i}"), 3) for i in (1, 2, 3)]))
    tpiba = 2.0 * np.pi / alat
    kpoints = tuple(_kpoint(n, i, counts, tpiba) for i, n in enumerate(bands.findall("ks_energies"), 1))
    nks = int(_text(bands, "nks"))
    if nks != len(kpoints) or nks < 1:
        raise ValueError("saved XML k count mismatch")
    smearing = bands.find("smearing")
    degauss = float(smearing.attrib["degauss"]) if smearing is not None else None
    fermi = bands.find("two_fermi_energies")
    if fermi is None:
        fermi = bands.find("fermi_energy")
    ef = _values(fermi, 2 if fermi.tag == "two_fermi_energies" else 1) if fermi is not None else _readonly(np.empty(0))
    totals = {n.tag: float(n.text) for n in _node(out, "total_energy")}
    mesh = bands.find("starting_k_points/monkhorst_pack")
    return SavedReference(
        source, creator.attrib["VERSION"], "Ha", nat, species, names, indices,
        positions, cell, alat, reciprocal, _readonly(tpiba * reciprocal),
        _grid(basis, "fft_grid"), _grid(basis, "fft_smooth"), _grid(basis, "fft_box"),
        int(_text(basis, "ngm")), int(_text(basis, "ngms")), int(_text(basis, "npwx")),
        _bool(basis, "gamma_only"), lsda, noncolin, spinorbit, counts, nks, kpoints,
        float(_text(bands, "nelec")), _text(bands, "occupations_kind"),
        (smearing.text or "").strip() if smearing is not None else None,
        degauss, 2 * degauss if degauss is not None else None, ef, _readonly(2 * ef),
        MappingProxyType(totals), MappingProxyType({k: 2 * v for k, v in totals.items()}),
        float(_text(basis, "ecutwfc")), float(_text(basis, "ecutrho")),
        MappingProxyType({k: int(v) for k, v in mesh.attrib.items()} if mesh is not None else {}),
        parse_saved_symmetries(root, source, nat, noncolin=noncolin),
    )
