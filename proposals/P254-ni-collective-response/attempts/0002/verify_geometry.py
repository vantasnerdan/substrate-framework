"""Independent parity-lattice oracle for preregistered W1 seed inputs only."""

import argparse
from copy import deepcopy
import hashlib
from itertools import product
import json
import os
from pathlib import Path

import sympy as sp
from ase.io import read as read_atoms

from substrate_framework.verification import CheckLedger

TOL = sp.Rational(1, 10**12)  # Angstrom representation tolerance, not a force tolerance.


def rational(value):
    return sp.Rational(str(value))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def delta_periodic(left, right, lengths):
    return tuple((x - y + sp.Rational(length) / 2) % length - sp.Rational(length) / 2
                 for x, y, length in zip(left, right, lengths, strict=True))


def distance2(left, right, lengths):
    return sum(value**2 for value in delta_periodic(left, right, lengths))


def half_grid(position, half_a, period):
    integers = tuple(int(sp.floor(rational(x) / half_a + sp.Rational(1, 2)))
                     for x in position)
    if any(abs(rational(x) - n * half_a) > TOL
           for x, n in zip(position, integers, strict=True)):
        raise ValueError("Position is not on the half-lattice within declared representation error")
    return tuple(n % length for n, length in zip(integers, period, strict=True))


def symbol_sites(record, symbol, half_a, period):
    return [half_grid(position, half_a, period)
            for name, position in zip(record["symbols"], record["positions_angstrom"], strict=True)
            if name == symbol]


def valid_octahedra(hydrogen, nickel, period):
    if len(set(hydrogen)) != len(hydrogen):
        return False
    for point in hydrogen:
        distances = [distance2(point, site, period) for site in nickel]
        if min(distances) != 1 or distances.count(1) != 6:
            return False
    return True


def verify_provenance(ledger, payload):
    for name, expected in payload["input_files_sha256"].items():
        ledger.check(f"bound producer input {name}", digest(name) == expected)
    for name, record in payload["structures"].items():
        physical = read_atoms(record["path"], format="extxyz")
        ledger.check(f"analyzed geometry is actual SCF-parsed seed {name}",
                     physical.get_chemical_symbols() == record["symbols"]
                     and physical.positions.tolist() == record["positions_angstrom"]
                     and physical.cell.array.tolist() == record["cell_angstrom"]
                     and physical.get_masses().tolist() == record["masses_amu"]
                     and not physical.constraints)
        expected_h = 6 if name in {"loaded-ni31h6", "dispersed-ni32h6"} else 0
        ledger.check(f"declared chemical inventory {name}",
                     record["symbols"].count("H") == expected_h
                     and set(record["symbols"]) <= {"Ni", "H"})
        ledger.check(f"actual saved seed hash {name}", digest(record["path"]) == record["sha256"])
        ledger.check(f"unconstrained periodic seed {name}",
                     record["constraint_count"] == 0 and record["pbc"] == [True] * 3)
        ledger.check(f"declared isotope/mean masses {name}", all(
            rational(mass) == rational(payload["config"][
                "ni_mass_amu" if symbol == "Ni" else "h1_mass_amu"]) and rational(mass) > 0
            for symbol, mass in zip(record["symbols"], record["masses_amu"], strict=True)
        ))


def verify_parent_and_vacancy(ledger, payload):
    config = payload["config"]
    half_a = rational(config["a_angstrom"]) / 2
    period = tuple(2 * n for n in config["conventional_repeat"])
    # Independent representation: the even-parity integer half-grid, not the
    # constructor's four-site FCC basis or its ASE primitive/conventional choice.
    nickel = {point for point in product(*(range(n) for n in period)) if sum(point) % 2 == 0}
    for name, record in payload["structures"].items():
        expected = nickel - {(0, 0, 0)} if name in {"empty-vacancy-ni31", "loaded-ni31h6"} else nickel
        sites = symbol_sites(record, "Ni", half_a, period)
        ledger.check(f"full FCC parent/vacancy identity {name}", set(sites) == expected and len(sites) == len(expected))
        expected_cell = sp.diag(*(half_a * n for n in period))
        actual_cell = sp.Matrix(record["cell_angstrom"]).applyfunc(rational)
        ledger.check(f"fixed conventional supercell {name}",
                     max(map(abs, actual_cell - expected_cell)) <= TOL)
    return half_a, period, nickel


def verify_disperse_selection(ledger, payload, half_a, period, nickel):
    record = payload["structures"]["dispersed-ni32h6"]
    hydrogen = symbol_sites(record, "H", half_a, period)
    ledger.check("six distinct six-coordinated octahedra", len(hydrogen) == 6 and valid_octahedra(hydrogen, nickel, period))
    octahedra = {point for point in product(*(range(n) for n in period)) if sum(point) % 2 == 1}
    ledger.check("registered first body centre", hydrogen[0] == (1, 1, 1))
    selected = [hydrogen[0]]
    for point in hydrogen[1:]:
        score = min(distance2(point, other, period) for other in selected)
        possible = [min(distance2(candidate, other, period) for other in selected)
                    for candidate in octahedra - set(selected)]
        ledger.check(f"farthest-point score at placement {len(selected) + 1}", score == max(possible))
        selected.append(point)
    wrong_scale = [tuple((2 * coordinate) % length for coordinate, length in zip(point, period, strict=True))
                   for point in hydrogen]
    ledger.check("supercell-fraction scale mutation is exposed", not valid_octahedra(wrong_scale, nickel, period))
    wrong_vacancy = deepcopy(payload["structures"]["empty-vacancy-ni31"])
    wrong_vacancy["positions_angstrom"][0] = [0, 0, 0]
    ledger.check("wrong vacancy with unchanged atom count is exposed",
                 set(symbol_sites(wrong_vacancy, "Ni", half_a, period)) != nickel - {(0, 0, 0)})


def inspect_loaded_shell(ledger, payload, half_a, period):
    record = payload["structures"]["loaded-ni31h6"]
    lengths = tuple(half_a * n for n in period)
    hydrogen = [tuple(map(rational, position))
                for name, position in zip(record["symbols"], record["positions_angstrom"], strict=True) if name == "H"]
    nickel = [tuple(map(rational, position))
              for name, position in zip(record["symbols"], record["positions_angstrom"], strict=True) if name == "Ni"]
    centered = [delta_periodic(point, (0, 0, 0), lengths) for point in hydrogen]
    radius = max(map(abs, centered[0]))
    axes = [tuple(sign * radius if index == axis else 0 for index in range(3))
            for axis in range(3) for sign in (-1, 1)]
    ledger.check("six distinct axial archived shell positions",
                 len(centered) == 6 and all(any(max(abs(x - y) for x, y in zip(point, target, strict=True)) <= TOL
                                                   for target in axes) for point in centered)
                 and len(set(centered)) == 6 and radius > 0)
    hh2 = [distance2(left, right, lengths) for i, left in enumerate(hydrogen) for right in hydrogen[i + 1:]]
    ledger.check("shell has no periodic H collisions", min(hh2) > 0)
    return {
        "scope": "Initial archived seed only; not a stationary structure or bonding assignment",
        "radius_angstrom": float(radius),
        "h_h_distances_angstrom": sorted(float(sp.sqrt(value)) for value in hh2),
        "h_ni_sorted_distances_angstrom": [sorted(float(sp.sqrt(distance2(h, ni, lengths))) for ni in nickel) for h in hydrogen],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    ledger = CheckLedger("P254-W1-seed-geometry")
    verify_provenance(ledger, payload)
    half_a, period, nickel = verify_parent_and_vacancy(ledger, payload)
    verify_disperse_selection(ledger, payload, half_a, period, nickel)
    inspection = inspect_loaded_shell(ledger, payload, half_a, period)
    ledger.finish()
    receipt = {
        "schema_version": 1, "status": "complete", "checks": ledger.passed,
        "evidence_scope": "Exact FCC/octahedral combinatorics with explicit floating-geometry representation tolerance and raw seed inspection",
        "does_not_license": ["stationarity", "mode stability", "coherence", "nuclear enhancement"],
        "source_sha256": {str(path.resolve()): digest(path) for path in (args.input, Path(__file__))},
        "producer_input_sha256": payload["input_sha256"], "inspection": inspection,
        "sympy_version": sp.__version__,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp")
    with temporary.open("w") as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, args.output)
    descriptor = os.open(args.output.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    print(f"P254_GEOMETRY_COMPLETE {args.output}", flush=True)


if __name__ == "__main__":
    main()
