"""Audit complete QE 7.5 SCF XML against the declared spin-PBE ensemble."""

from itertools import product
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from ase import Atoms
from ase.units import create_units

from .qe_io import SCF_MIXING_BETA, SCF_MIXING_MODE, sha256_file, validate_geometry
from .qe_scf_pseudo import pseudo_header

# Match ASE espresso-out's explicit PWSCF CODATA-2006 conversion.
Bohr, Hartree = (create_units("2006")[key] for key in ("Bohr", "Hartree"))


def node(parent, path):
    found = parent.findall(path)
    if len(found) != 1:
        raise ValueError(f"Expected exactly one SCF XML {path}")
    return found[0]


def text(parent, path):
    return (node(parent, path).text or "").strip()


def numbers(value, count):
    result = np.array([float(x.replace("D", "E")) for x in value.split()])
    if result.size != count or not np.isfinite(result).all():
        raise ValueError("Incomplete/nonfinite SCF XML array")
    return result


def scalar(parent, path):
    return float(numbers(text(parent, path), 1)[0])


def equal(actual, expected, name, atol=1e-10):
    if not np.allclose(actual, expected, rtol=0, atol=atol):
        raise ValueError(f"SCF XML {name} differs from declaration")


def _geometry(parent, expected):
    structure = node(parent, "atomic_structure")
    atoms = node(structure, "atomic_positions").findall("atom")
    if int(structure.attrib["nat"]) != len(expected) or len(atoms) != len(expected):
        raise ValueError("SCF XML atom count differs")
    positions = np.array([numbers(a.text or "", 3) for a in atoms]) * Bohr
    cell = np.array([numbers(text(structure, f"cell/a{i}"), 3) for i in (1, 2, 3)]) * Bohr
    actual = Atoms([a.attrib["name"] for a in atoms], positions=positions, cell=cell, pbc=True)
    validate_geometry(actual, expected)
    return {"symbols": actual.get_chemical_symbols(), "positions_angstrom": positions.tolist(),
            "cell_angstrom": cell.tolist()}


def _species(parent, spec, expected):
    species = node(parent, "atomic_species").findall("species")
    names = [s.attrib["name"] for s in species]
    if len(set(names)) != len(names) or set(names) != set(expected.get_chemical_symbols()):
        raise ValueError("SCF XML species differ")
    records = {}
    for entry in species:
        symbol = entry.attrib["name"]
        filename = text(entry, "pseudo_file")
        if filename != spec["execution"]["pseudopotentials"][symbol]:
            raise ValueError("SCF XML pseudopotential identity differs")
        path = Path(spec["execution"]["pseudo_dir"]) / filename
        header = pseudo_header(path)
        if header["element"].strip() != symbol:
            raise ValueError("Invalid UPF species identity")
        records[symbol] = {"filename": filename, "sha256": sha256_file(path),
                           "valence": float(header["z_valence"]), "header": header}
    return records


def _bands(root, spec, expected, species):
    s = spec["settings"]
    bands = node(root, "output/band_structure")
    for path in ("input/spin", "output/magnetization", "output/band_structure"):
        if [text(root, path + "/" + key) for key in ("lsda", "noncolin", "spinorbit")] != ["true", "false", "false"]:
            raise ValueError("SCF XML spin contract differs")
    for path in ("input/bands", "output/band_structure"):
        smear = node(root, path + "/smearing")
        if (smear.text or "").strip() != s["smearing"]:
            raise ValueError("SCF XML occupation smearing differs")
        equal(float(smear.attrib["degauss"]) * 2, s["degauss_ry"], "degauss")
    if text(root, "input/bands/occupations") != "smearing" or text(bands, "occupations_kind") != "smearing":
        raise ValueError("SCF XML occupations kind differs")
    for path in ("input/bands/nbnd", "output/band_structure/nbnd_up", "output/band_structure/nbnd_dw"):
        equal(scalar(root, path), s["nbnd"], "band count", atol=0)
    nelec = sum(species[symbol]["valence"] for symbol in expected.get_chemical_symbols())
    equal(scalar(root, "input/bands/tot_charge"), 0, "neutral charge")
    equal(scalar(bands, "nelec"), nelec, "valence electrons")
    mesh_size = int(np.prod(s["kpts"]))
    nks_value = scalar(bands, "nks")
    nks = int(nks_value)
    if nks_value != nks or not 1 <= nks <= mesh_size:
        raise ValueError("Invalid SCF XML k mesh size")
    full_mesh = s["nosym"] and s["noinv"]
    if full_mesh and nks != mesh_size:
        raise ValueError("SCF XML full k mesh differs")
    entries = bands.findall("ks_energies")
    if len(entries) != nks:
        raise ValueError("Incomplete SCF k-point records")
    reciprocal = np.linalg.inv(expected.cell.array).T
    alat = float(node(root, "output/atomic_structure").attrib["alat"]) * Bohr
    mesh, records, population = [], [], 0.0
    ef = scalar(bands, "fermi_energy")
    for entry in entries:
        kp = node(entry, "k_point")
        cart = numbers(kp.text or "", 3)
        fractional = cart @ np.linalg.inv(reciprocal * alat)
        mesh.append(np.mod(fractional, 1))
        weight = float(kp.attrib["weight"])
        if not np.isfinite(weight) or weight <= 0:
            raise ValueError("Invalid SCF XML k weight")
        equal(weight * mesh_size, round(weight * mesh_size), "mesh weight multiplicity")
        if full_mesh:
            equal(weight, 1 / nks, "k weight")
        arrays = {}
        for key in ("eigenvalues", "occupations"):
            array = node(entry, key)
            equal(int(array.attrib["size"]), 2 * s["nbnd"], key + " size", atol=0)
            arrays[key] = numbers(array.text or "", 2 * s["nbnd"])
        occupation = arrays["occupations"]
        if s["smearing"] == "fd":
            argument = np.clip((arrays["eigenvalues"] - ef) / (s["degauss_ry"] / 2), -700, 700)
            equal(occupation, 1 / (1 + np.exp(argument)), "FD occupations", atol=2e-9)
        else:
            # QE7.5 Modules/wgauss.f90:59-65, x=(Ef-E)/degauss.
            x = (ef - arrays["eigenvalues"]) / (s["degauss_ry"] / 2) - 1 / math.sqrt(2)
            mv = np.array([0.5 * math.erf(float(v)) + math.exp(-min(200, float(v * v))) /
                           math.sqrt(2 * math.pi) + 0.5 for v in x])
            equal(occupation, mv, "MV occupations", atol=2e-9)
        population += weight * occupation.sum()
        records.append({"k_cartesian_2pi_alat": cart.tolist(), "k_fractional": fractional.tolist(),
                        "weight": weight, "eigenvalues_hartree": arrays["eigenvalues"].reshape(2, -1).tolist(),
                        "occupations": occupation.reshape(2, -1).tolist()})
    # Compare modulo reciprocal lattice vectors; native ordering is retained above.
    def canonical(rows):
        return sorted(tuple(np.round(np.mod(r, 1) % 1, 7)) for r in rows)
    target = [np.array(row) / s["kpts"] for row in product(*(range(n) for n in s["kpts"]))]
    # Remove floating negative-zero / modulo-one roundoff before sorting.
    cleaned = [np.where(np.isclose(row, 1, atol=1e-7), 0, row) for row in mesh]
    observed, declared = canonical(cleaned), canonical(target)
    if (len(set(observed)) != nks or not set(observed) <= set(declared)
            or (full_mesh and observed != declared)):
        raise ValueError("SCF XML k coordinates differ")
    equal(sum(record["weight"] for record in records), 1, "k weight normalization")
    equal(population, nelec, "occupation electron sum", atol=2e-7)
    return {"nbnd_per_spin": s["nbnd"], "nelec": nelec, "occupation_electron_sum": population,
            "nks": nks, "fermi_energy_hartree": ef, "ks_energies": records}


def read_scf_xml(path, version):
    """Read a complete native XML envelope for the declared QE version."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise ValueError(f"Incomplete SCF XML: {exc}") from exc
    if root.tag != "{http://www.quantum-espresso.org/ns/qes/qes-1.0}espresso" or root.attrib.get("Units") != "Hartree atomic units":
        raise ValueError("Unexpected SCF XML format/units")
    if text(root, "exit_status") != "0" or root.find("closed") is None:
        raise ValueError("Unfinished SCF XML")
    if node(root, "general_info/creator").attrib != {"NAME": "PWSCF", "VERSION": version}:
        raise ValueError("SCF XML creator differs")
    return root


def validate_scf_metadata(root, spec, expected):
    """Audit shared geometry, model, solver and occupied-band metadata."""
    s = spec["settings"]
    for path in ("input/dft/functional", "output/dft/functional"):
        if text(root, path) != "PBE":
            raise ValueError("SCF XML functional differs")
    for flag in ("nosym", "noinv"):
        if text(root, "input/symmetry_flags/" + flag) != str(s[flag]).lower():
            raise ValueError("SCF XML symmetry contract differs")
    for path in ("input/k_points_IBZ/monkhorst_pack", "output/band_structure/starting_k_points/monkhorst_pack"):
        mp = node(root, path).attrib
        if [int(mp[k]) for k in ("nk1", "nk2", "nk3", "k1", "k2", "k3")] != s["kpts"] + [0, 0, 0]:
            raise ValueError("SCF XML requested k mesh differs")
    for key in ("ecutwfc", "ecutrho"):
        equal(scalar(root, "input/basis/" + key) * 2, s[key + "_ry"], key)
        equal(scalar(root, "output/basis_set/" + key) * 2, s[key + "_ry"], key)
    equal(scalar(root, "input/electron_control/conv_thr") * 2, s["conv_thr_ry"], "conv_thr", atol=1e-20)
    equal(scalar(root, "input/electron_control/max_nstep"), s["electron_maxstep"], "electron_maxstep", atol=0)
    equal(scalar(root, "input/electron_control/mixing_beta"), SCF_MIXING_BETA, "mixing_beta")
    if text(root, "input/electron_control/mixing_mode") != SCF_MIXING_MODE:
        raise ValueError("SCF XML mixing mode differs")
    if text(root, "output/convergence_info/scf_conv/convergence_achieved") != "true":
        raise ValueError("SCF XML lacks convergence")
    for path in ("input", "output"):
        geometry = _geometry(node(root, path), expected)
    species = _species(node(root, "input"), spec, expected)
    output_species = _species(node(root, "output"), spec, expected)
    if species != output_species:
        raise ValueError("SCF XML input/output pseudo identity differs")
    for entry in node(root, "input/atomic_species").findall("species"):
        requested = s["starting_magnetization_ni"] if entry.attrib["name"] == "Ni" else 0
        equal(scalar(entry, "starting_magnetization"), requested, "starting magnetization")
    return {"species": species, "geometry": geometry,
            "bands": _bands(root, spec, expected, species)}


def parse_scf_xml(path, spec, expected, physical):
    """Parse public XML independently of stdout and refuse mismatched science."""
    root = read_scf_xml(path, physical["qe_version"])
    metadata = validate_scf_metadata(root, spec, expected)
    equal(scalar(root, "output/total_energy/etot") * Hartree, physical["energy_ev"], "energy", atol=2e-6)
    forces = numbers(text(root, "output/forces"), 3 * len(expected)).reshape(-1, 3) * Hartree / Bohr
    equal(forces, physical["forces_ev_per_angstrom"], "forces", atol=3e-6)
    s = spec["settings"]
    return {"model_contract": {"generation": "explicit-ensemble-symmetry-v3", "functional": "PBE",
             "nspin": 2, "nosym": s["nosym"], "noinv": s["noinv"], "smearing": s["smearing"],
             "degauss_ry": s["degauss_ry"], "occupation_interpretation":
             "declared one-electron Fermi-Dirac ensemble; not measured temperature" if s["smearing"] == "fd"
             else "numerical Marzari-Vanderbilt smearing; not a physical ensemble or loss"},
            **metadata,
            "native_timing_seconds": {k: scalar(root, "timing_info/total/" + k) for k in ("cpu", "wall")},
            "converged": True, "xml_units": "Hartree atomic units"}
